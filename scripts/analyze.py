"""Landmark survival associations and censoring-adjusted prediction validation."""
from pathlib import Path
import json
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.exceptions import ConvergenceError, ConvergenceWarning
from lifelines.statistics import proportional_hazard_test
from scipy.optimize import minimize, brentq
from scipy.special import logit, expit
from scipy.stats import chi2
from sklearn.ensemble import (ExtraTreesClassifier, RandomForestClassifier,
                              GradientBoostingClassifier, HistGradientBoostingClassifier)
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.tree import DecisionTreeClassifier

from prepare_data import eligible, load_all

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
SEED = 11202026
COHORTS = ("TCGA", "AOCS", "UKE", "ICON7")
BASE = ["age10", "iv", "rd"]
STATUS = BASE + ["early"]
TIMING = STATUS + ["timing6"]
FEATURE_LABELS = {"age10": "Age", "iv": "Stage IV", "rd": "Residual >1 cm",
                  "early": "Progression status", "timing6": "Time since progression"}


def cohort_data(raw, landmark=12, include_age=True, include_rd=True,
                exact=True, documented=False, arms=None):
    required = [k for k, use in (("age", include_age), ("rd", include_rd)) if use]
    cohorts = COHORTS if include_age else (*COHORTS, "Japan")
    frames, flow = {}, []
    for cohort in cohorts:
        frame = raw[cohort]
        if arms and cohort in arms:
            frame = frame.loc[frame.treatment.eq(arms[cohort])]
        data, steps = eligible(frame, landmark, require=required, inclusive=exact)
        if documented and cohort == "TCGA":
            data = data.loc[data.early.eq(0) | data.documented_event].reset_index(drop=True)
        data["age10"] = (data.age - 60) / 10
        data["timing6"] = np.where(data.early.eq(1), (data.duration - 6) / 6, 0)
        data["group"] = np.select([data.early.eq(0), data.duration.lt(6)], [0, 1], default=2)
        frames[cohort] = data
        flow.extend({"cohort": cohort, "step": k, "n": v} for k, v in steps.items())
    return frames, flow


def censor_curve(time, event):
    unique, idx, counts = np.unique(time, return_inverse=True, return_counts=True)
    cens = np.bincount(idx, weights=1 - event, minlength=len(unique))
    risk = len(time) - np.r_[0, counts.cumsum()[:-1]]
    return unique, np.r_[1.0, np.cumprod(1 - cens / risk)]


def horizon_data(data, horizon):
    time = data.time.to_numpy()
    event = data.death.to_numpy()
    return horizon_arrays(time, event, horizon)


def horizon_arrays(time, event, horizon):
    t, g = censor_curve(time, event)
    cases = (event == 1) & (time <= horizon)
    controls = time > horizon
    weights = np.zeros(len(time))
    weights[cases] = 1 / g[np.searchsorted(t, time[cases], side="left")]
    gh = g[np.searchsorted(t, horizon, side="right")]
    if controls.any():
        weights[controls] = 1 / gh
    return cases.astype(int), weights, gh, int(controls.sum())


def log_model(c=100.0):
    return LogisticRegression(C=c, solver="lbfgs", max_iter=2000)


def algorithms():
    return {
        "Logistic": log_model(),
        "LASSO logistic": LogisticRegression(C=1, l1_ratio=1, solver="saga", max_iter=5000, random_state=SEED),
        "Elastic net": LogisticRegression(C=1, l1_ratio=.5, solver="saga", max_iter=5000, random_state=SEED),
        "Decision tree": DecisionTreeClassifier(max_depth=3, min_samples_leaf=15, random_state=SEED),
        "Random forest": RandomForestClassifier(n_estimators=300, min_samples_leaf=12, max_features=None, random_state=SEED),
        "Extra trees": ExtraTreesClassifier(n_estimators=300, min_samples_leaf=12, random_state=SEED),
        "Gradient boosting": GradientBoostingClassifier(n_estimators=80, max_depth=2, min_samples_leaf=15,
                                                       learning_rate=.04, random_state=SEED),
        "Hist. gradient boosting": HistGradientBoostingClassifier(max_iter=80, max_leaf_nodes=7,
                                                                 min_samples_leaf=20, learning_rate=.04,
                                                                 random_state=SEED),
    }


def fit_predict(model, training, testing, columns, horizon):
    y, weights, _, _ = horizon_data(training, horizon)
    observed = weights > 0
    model.fit(training.loc[observed, columns], y[observed], sample_weight=weights[observed])
    return {key: model.predict_proba(data[columns])[:, 1] for key, data in testing.items()}


def calibration(y, w, p):
    lp = logit(np.clip(p, 1e-6, 1 - 1e-6))
    if np.ptp(lp[w > 0]) < 1e-9:
        return np.nan, np.nan
    intercept = brentq(lambda a: np.dot(w, expit(lp + a) - y), -30, 30)
    design = np.column_stack((np.ones(len(lp)), lp))

    def loss(beta):
        fitted = design @ beta
        value = np.dot(w, np.logaddexp(0, fitted) - y * fitted)
        gradient = design.T @ (w * (expit(fitted) - y))
        return value, gradient

    slope = minimize(loss, [intercept, 1], jac=True, method="BFGS").x[1]
    return intercept, slope


def metrics(data, p, horizon):
    y, w, g, controls = horizon_data(data, horizon)
    intercept, slope = calibration(y, w, p) if y.sum() and controls else (np.nan, np.nan)
    return {
        "n": len(y), "deaths": int(y.sum()), "controls": controls, "G": g,
        "auc": roc_auc_score(y, p, sample_weight=w) if y.sum() and controls else np.nan,
        "brier": np.dot(w, (y - p) ** 2) / len(y),
        "observed": np.dot(w, y) / len(y), "expected": p.mean(),
        "oe": (np.dot(w, y) / len(y)) / p.mean(),
        "intercept": intercept, "slope": slope,
    }


def compare_bootstrap(frame, preds, horizon, repetitions=1000):
    rng = np.random.default_rng(SEED + len(frame))
    time, event = frame.time.to_numpy(), frame.death.to_numpy()
    keys = list(preds)
    rng_values = []
    for _ in range(repetitions):
        idx = rng.integers(0, len(time), len(time))
        y, w, _, controls = horizon_arrays(time[idx], event[idx], horizon)
        if not y.sum() or not controls:
            continue
        p = {k: preds[k][idx] for k in keys}
        auc = [roc_auc_score(y, p[k], sample_weight=w) for k in keys]
        brier = [np.dot(w, (y - p[k]) ** 2) / len(y) for k in keys]
        rng_values.append([*auc, auc[1] - auc[0], auc[2] - auc[1],
                           brier[1] - brier[0], brier[2] - brier[1]])
    columns = [*keys, "Status vs baseline", "Timing vs status", "Brier status vs baseline", "Brier timing vs status"]
    first, last = np.percentile(rng_values, [2.5, 97.5], axis=0)
    return {key: (lo, hi) for key, lo, hi in zip(columns, first, last, strict=True)}


def saved_predictions(results, frames, variant, horizon):
    saved = []
    for cohort, values in results.items():
        for model, predictions in values.items():
            saved.append(pd.DataFrame({"cohort": cohort, "variant": variant, "model": model,
                                       "id": frames[cohort].id, "time": frames[cohort].time,
                                       "death": frames[cohort].death, "risk": predictions,
                                       "horizon": horizon}))
    return pd.concat(saved, ignore_index=True)


def run_models(frames, horizon, features=BASE, extra=False, oof=False):
    sets = {"Baseline": features, "Status": features + ["early"], "Timing": features + ["early", "timing6"]}
    models = {k: (log_model(), cols) for k, cols in sets.items()}
    if extra:
        models.update({k: (model, features + ["early", "timing6"]) for k, model in algorithms().items()})
    train = frames["TCGA"]
    result = {k: {} for k in frames}
    for name, (model, cols) in models.items():
        fitted = fit_predict(model, train, frames, cols, horizon)
        for cohort, predictions in fitted.items():
            result[cohort][name] = predictions
    if oof:
        y, w, _, _ = horizon_data(train, horizon)
        strata = y + 2 * (w == 0)
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
        for name, (model, cols) in models.items():
            pred = np.empty(len(train))
            for learn, test in cv.split(train, strata):
                pred[test] = fit_predict(model, train.iloc[learn], {"test": train.iloc[test]}, cols, horizon)["test"]
            result["TCGA"][name] = pred
    return result


def cox(data, cols):
    f = data[cols + ["time", "death"]].copy()
    f["death"] = f.death.eq(1) & f.time.le(24)
    f["time"] = f.time.clip(upper=24)
    return CoxPHFitter().fit(f, duration_col="time", event_col="death"), f


def cox_stats(frames):
    rows, ph = [], []
    for cohort, d in frames.items():
        for termset, cols in (("Status", STATUS), ("Timing", TIMING)):
            mod, f = cox(d, cols)
            for term, entry in mod.summary.iterrows():
                rows.append({"cohort": cohort, "model": termset, "term": term,
                             "hr": entry["exp(coef)"], "lo": entry["exp(coef) lower 95%"],
                             "hi": entry["exp(coef) upper 95%"], "p": entry["p"],
                             "events": int(f.death.sum())})
            for term, p in proportional_hazard_test(mod, f, time_transform="rank").summary.p.items():
                ph.append({"cohort": cohort, "model": termset, "term": term, "p": p})
    return pd.DataFrame(rows), pd.DataFrame(ph)


def subgroup_stats(frames):
    groups = {"<60 years": ("age10", lambda d: d.age.lt(60)), "≥60 years": ("age10", lambda d: d.age.ge(60)),
              "Stage III": ("iv", lambda d: d.iv.eq(0)), "Stage IV": ("iv", lambda d: d.iv.eq(1)),
              "Residual ≤1 cm": ("rd", lambda d: d.rd.eq(0)), "Residual >1 cm": ("rd", lambda d: d.rd.eq(1))}
    rows = []
    for cohort, data in frames.items():
        main, _ = cox(data, STATUS)
        for label, (omit, condition) in groups.items():
            d = data.loc[condition(data)].copy()
            d_int = data.assign(interaction=data.early * (data.age10.ge(0) if omit == "age10" else data[omit]))
            with warnings.catch_warnings():
                warnings.simplefilter("error", ConvergenceWarning)
                try:
                    if omit == "age10":
                        d_int["age60"] = data.age10.ge(0).astype(int)
                        m0, _ = cox(d_int, STATUS + ["age60"])
                        m1, _ = cox(d_int, STATUS + ["age60", "interaction"])
                    else:
                        m0 = main
                        m1, _ = cox(d_int, STATUS + ["interaction"])
                    pint = chi2.sf(max(0, 2 * (m1.log_likelihood_ - m0.log_likelihood_)), 1)
                except (ConvergenceError, ConvergenceWarning):
                    pint = np.nan
                try:
                    cols = [c for c in STATUS if c != omit and d[c].nunique() > 1]
                    model, _ = cox(d, cols)
                    est = model.summary.loc["early"]
                    hr, lo, hi = est[["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%"]]
                except (ConvergenceError, ConvergenceWarning, KeyError):
                    hr, lo, hi = (np.nan,) * 3
            deaths = (d.death.eq(1) & d.time.le(24)).sum()
            rows.append({"cohort": cohort, "group": label, "n": len(d), "events": int(deaths),
                         "hr": hr, "lo": lo, "hi": hi, "interaction_p": pint})
    return pd.DataFrame(rows)


def describe(frames):
    rows, survival = [], []
    for cohort, f in frames.items():
        for group, d in [("All", f), *list(f.groupby("group"))]:
            if not len(d):
                continue
            km = KaplanMeierFitter().fit(d.time, d.death)
            y, _, g, controls = horizon_data(d, 12)
            row = {"cohort": cohort, "group": group, "n": len(d), "age_median": d.age.median(),
                   "age_q1": d.age.quantile(.25), "age_q3": d.age.quantile(.75),
                   "stage_iv": int(d.iv.sum()), "suboptimal": int(d.rd.sum()), "early": int(d.early.sum()),
                   "bevacizumab": int(d.treatment.eq("Bevacizumab").sum()), "dead12": int(y.sum()),
                   "known_alive12": controls, "dep_followup12": g,
                   "followup_median": KaplanMeierFitter().fit(d.time, 1 - d.death).median_survival_time_}
            for h in (12, 24):
                last = km.confidence_interval_.loc[:h].iloc[-1]
                row.update({f"risk{h}": 1 - km.predict(h), f"risk{h}_lo": 1 - last.iloc[1],
                            f"risk{h}_hi": 1 - last.iloc[0]})
            rows.append(row)
            curve = pd.concat([km.survival_function_, km.confidence_interval_], axis=1).reset_index()
            curve.columns = ["months", "survival", "lo", "hi"]
            survival.append(curve.assign(cohort=cohort, group=group))
    return pd.DataFrame(rows), pd.concat(survival, ignore_index=True)


def run():
    OUT.mkdir(exist_ok=True)
    raw = load_all()
    primary, flow = cohort_data(raw)
    pd.DataFrame(flow).to_csv(OUT / "flow.csv", index=False)
    # Only the internal validation is out of fold. All external predictions
    # are from the untouched full development fit.
    result = run_models(primary, 12, extra=True, oof=True)
    report, contrasts, stored = [], [], []
    contrast_keys = ["Status vs baseline", "Timing vs status", "Brier status vs baseline", "Brier timing vs status"]

    def report_variant(key, frames, horizon, results, boot=False):
        for cohort, predictions in results.items():
            estimates = {name: metrics(frames[cohort], pred, horizon) for name, pred in predictions.items()}
            ci = compare_bootstrap(frames[cohort], {k: predictions[k] for k in ("Baseline", "Status", "Timing")},
                                   horizon, 1000 if boot else 350) if cohort != "TCGA" else {}
            for model, estimate in estimates.items():
                low, high = ci.get(model, [np.nan, np.nan])
                report.append({"variant": key, "cohort": cohort, "model": model, **estimate,
                               "auc_lo": low, "auc_hi": high})
            for contrast in contrast_keys:
                field = "brier" if contrast.startswith("Brier") else "auc"
                new, old = ("Status", "Baseline") if "baseline" in contrast else ("Timing", "Status")
                low, high = ci.get(contrast, [np.nan, np.nan])
                contrasts.append({"variant": key, "cohort": cohort, "contrast": contrast,
                                  "delta": estimates[new][field] - estimates[old][field],
                                  "lo": low, "hi": high})
        stored.append(saved_predictions(results, frames, key, horizon))

    report_variant("Primary", primary, 12, result, boot=True)
    variants = [
        ("Landmark 6", 12, {"landmark": 6}, BASE),
        ("Landmark 18", 12, {"landmark": 18}, BASE),
        ("Horizon 24", 24, {}, BASE),
        ("Progression strictly before landmark", 12, {"exact": False}, BASE),
        ("Documented events in TCGA", 12, {"documented": True}, BASE),
        ("No age", 12, {"include_age": False}, ["iv", "rd"]),
        ("No residual disease variable", 12, {"include_rd": False}, ["age10", "iv"]),
        ("Standard treatment", 12, {"arms": {"AOCS": "Platinum", "UKE": "Platinum", "ICON7": "Platinum"}}, BASE),
    ]
    for name, horizon, opts, features in variants:
        data, _ = cohort_data(raw, **opts)
        report_variant(name, data, horizon, run_models(data, horizon, features, oof=True))
        print("Completed", name, flush=True)

    pd.DataFrame(report).to_csv(OUT / "model_metrics.csv", index=False)
    pd.DataFrame(contrasts).to_csv(OUT / "model_contrasts.csv", index=False)
    pd.concat(stored, ignore_index=True).to_csv(OUT / "predictions.csv", index=False)
    description, curves = describe(primary)
    description.to_csv(OUT / "cohort_characteristics.csv", index=False)
    curves.to_csv(OUT / "survival_curves.csv", index=False)
    associations, ph = cox_stats(primary)
    associations.to_csv(OUT / "cox_associations.csv", index=False)
    ph.to_csv(OUT / "cox_ph_tests.csv", index=False)
    subgroup_stats(primary).to_csv(OUT / "subgroups.csv", index=False)
    dev = primary["TCGA"]
    eq = {}
    for key, cols in (("Baseline", BASE), ("Status", STATUS), ("Timing", TIMING)):
        model = log_model()
        fit_predict(model, dev, {}, cols, 12)
        eq[key] = {"intercept": model.intercept_[0],
                   **dict(zip(cols, model.coef_[0], strict=True))}
    (OUT / "equations.json").write_text(json.dumps(eq, indent=2), encoding="utf-8")


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*penalty.*deprecated.*", category=FutureWarning)
        run()

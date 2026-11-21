"""Second-revision clinical analyses.

A  Post-progression survival among all patients with a recorded progression.
B  Multi-landmark supermodel and an illness-death (clock-forward) model.
C  Internal-external cross-validation, recalibration, multiple imputation, precision.
E  Smooth calibration, period-specific hazards and prognostic-enrichment arithmetic.
"""
import sys
import warnings

from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import proportional_hazard_test
from lifelines.utils import concordance_index
import numpy as np
import pandas as pd
from patsy import build_design_matrices, dmatrix
from scipy import stats
from scipy.optimize import brentq
from scipy.special import expit
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

from v2_common import (BASE, OUT, PRIMARY, SEED, SETS, STATUS, TIMING, all_sources, base_population,
                       clip_logit, cohort_data, fit_pooled, horizon_arrays, horizon_data, meta, metrics,
                       weighted_logistic)

warnings.filterwarnings("ignore")
LANDMARKS = (6, 9, 12, 15, 18, 21, 24)


def hr_row(fitted, term, /, **extra):
    s = fitted.summary.loc[term]
    return {**extra, "term": term, "coef": s["coef"], "se": s["se(coef)"], "hr": s["exp(coef)"],
            "lo": s["exp(coef) lower 95%"], "hi": s["exp(coef) upper 95%"], "p": s["p"]}


# ---------------------------------------------------------------- A. post-progression survival
def pps_frames(raw):
    frames, flow = {}, []
    for cohort, frame in raw.items():
        d = base_population(frame)
        prog = d.loc[d.progression_event.eq(1) & d.pfs.gt(0)].copy()
        prog["pps"] = prog.os - prog.pfs
        coincident = prog.pps.le(1e-6)
        flow.append({"cohort": cohort, "stage_serous_eligible": len(d), "progressed": len(prog),
                     "death_or_last_contact_at_progression": int(coincident.sum())})
        prog = prog.loc[~coincident].copy()
        prog["ttp6"] = prog.pfs / 6
        prog["cohort"] = cohort
        flow[-1].update(analysed=len(prog), deaths=int(prog.death.sum()),
                        ttp_median=prog.pfs.median(), ttp_q1=prog.pfs.quantile(.25), ttp_q3=prog.pfs.quantile(.75),
                        ttp_le12=int(prog.pfs.le(12).sum()))
        frames[cohort] = prog
    pd.DataFrame(flow).to_csv(OUT / "pps_flow.csv", index=False)
    return frames


def adjusters(d):
    return [c for c in ("age10", "iv", "rd") if d[c].notna().all() and d[c].nunique() > 1]


def pps_analysis(raw):
    frames = pps_frames(raw)
    rows, cohort_rows = [], []
    for cohort, d in frames.items():
        d = d.dropna(subset=[c for c in ("age10", "rd") if d[c].notna().mean() > .5])
        cols = adjusters(d)
        m = CoxPHFitter().fit(d[cols + ["ttp6", "pps", "death"]], "pps", "death")
        cohort_rows.append(hr_row(m, "ttp6", cohort=cohort, n=len(d), deaths=int(d.death.sum()),
                                  adjusted_for="+".join(cols)))
        short = d.loc[d.pfs.le(12)]
        if short.death.sum() >= 10:
            cols = adjusters(short)
            m = CoxPHFitter().fit(short[cols + ["ttp6", "pps", "death"]], "pps", "death")
            cohort_rows.append(hr_row(m, "ttp6", cohort=cohort, n=len(short), deaths=int(short.death.sum()),
                                      adjusted_for="+".join(cols) + " | progression within 12 months"))
    cohort_table = pd.DataFrame(cohort_rows)
    cohort_table.to_csv(OUT / "pps_cohort_hr.csv", index=False)
    for label, part in (("All progressions", cohort_table.loc[~cohort_table.adjusted_for.str.contains("within")]),
                        ("Progression within 12 months", cohort_table.loc[cohort_table.adjusted_for.str.contains("within")])):
        for scope, sub in (("All cohorts", part), ("Primary cohorts", part.loc[part.cohort.isin(PRIMARY)])):
            r = meta(sub.coef, sub.se)
            rows.append({"analysis": f"Two-stage random effects, {label}", "scope": scope, "k": r["k"],
                         "n": int(sub.n.sum()), "deaths": int(sub.deaths.sum()), "hr": np.exp(r["estimate"]),
                         "lo": np.exp(r["lo"]), "hi": np.exp(r["hi"]), "i2": r["i2"],
                         "pi_lo": np.exp(r.get("pi_lo", np.nan)), "pi_hi": np.exp(r.get("pi_hi", np.nan))})

    primary = pd.concat([frames[c] for c in PRIMARY]).dropna(subset=["age10", "rd"]).reset_index(drop=True)
    everyone = pd.concat(frames.values()).reset_index(drop=True)
    base = CoxPHFitter().fit(primary[BASE + ["pps", "death", "cohort"]], "pps", "death", strata=["cohort"])
    full = CoxPHFitter().fit(primary[BASE + ["ttp6", "pps", "death", "cohort"]], "pps", "death", strata=["cohort"])
    lr = 2 * (full.log_likelihood_ - base.log_likelihood_)
    rows.append({"analysis": "One-stage stratified Cox, age, stage and residual adjusted", "scope": "Primary cohorts",
                 "k": 4, "n": len(primary), "deaths": int(primary.death.sum()), **{k: v for k, v in hr_row(full, "ttp6").items()
                                                                                    if k in ("hr", "lo", "hi", "p")},
                 "lr_chi2": lr, "lr_p": stats.chi2.sf(lr, 1)})
    wide = CoxPHFitter().fit(everyone[["iv", "ttp6", "pps", "death", "cohort"]], "pps", "death", strata=["cohort"])
    rows.append({"analysis": "One-stage stratified Cox, stage adjusted", "scope": "All cohorts", "k": everyone.cohort.nunique(),
                 "n": len(everyone), "deaths": int(everyone.death.sum()),
                 **{k: v for k, v in hr_row(wide, "ttp6").items() if k in ("hr", "lo", "hi", "p")}})
    ph = proportional_hazard_test(full, primary[BASE + ["ttp6", "pps", "death", "cohort"]], time_transform="rank")
    rows.append({"analysis": "Schoenfeld test for the progression-time term", "scope": "Primary cohorts",
                 "p": float(ph.summary.loc["ttp6", "p"])})
    # Early and late post-progression follow-up, because the timing effect need not be proportional.
    for name, lo, hi in (("0-12 months after progression", 0, 12), (">12 months after progression", 12, np.inf)):
        part = primary.loc[primary.pps.gt(lo)].copy()
        part["event"] = (part.death.eq(1) & part.pps.le(hi)).astype(int)
        part["stop"] = part.pps.clip(upper=hi if np.isfinite(hi) else None)
        part["start"] = lo
        m = CoxPHFitter().fit(part[BASE + ["ttp6", "start", "stop", "event", "cohort"]], "stop", "event",
                              entry_col="start", strata=["cohort"])
        rows.append({"analysis": f"Period-specific hazard ratio, {name}", "scope": "Primary cohorts", "n": len(part),
                     "deaths": int(part.event.sum()),
                     **{k: v for k, v in hr_row(m, "ttp6").items() if k in ("hr", "lo", "hi", "p")}})
    pd.DataFrame(rows).to_csv(OUT / "pps_pooled.csv", index=False)

    # Categories and Kaplan-Meier curves.
    bins = [0, 6, 12, 18, 24, np.inf]
    labels = ["<6", "6 to <12", "12 to <18", "18 to <24", "≥24"]
    everyone["band"] = pd.cut(everyone.pfs, bins, right=False, labels=labels)
    primary["band"] = pd.cut(primary.pfs, bins, right=False, labels=labels)
    cat = pd.get_dummies(primary.band, prefix="band", dtype=float).drop(columns="band_6 to <12")
    m = CoxPHFitter().fit(pd.concat([primary[BASE + ["pps", "death", "cohort"]], cat], axis=1), "pps", "death",
                          strata=["cohort"])
    band_rows, curves = [], []
    for label in labels:
        d = everyone.loc[everyone.band.eq(label)]
        km = KaplanMeierFitter().fit(d.pps, d.death)
        ci = km.confidence_interval_.loc[:12].iloc[-1]
        row = {"band": label, "n_all": len(d), "deaths_all": int(d.death.sum()), "median_pps": km.median_survival_time_,
               "mortality12": 1 - km.predict(12), "mortality12_lo": 1 - ci.iloc[1], "mortality12_hi": 1 - ci.iloc[0],
               "n_primary": int(primary.band.eq(label).sum())}
        if f"band_{label}" in m.summary.index:
            row.update({k: v for k, v in hr_row(m, f"band_{label}").items() if k in ("hr", "lo", "hi", "p")})
        else:
            row.update(hr=1.0)
        band_rows.append(row)
        curve = km.survival_function_.reset_index()
        curve.columns = ["months", "survival"]
        curves.append(curve.assign(band=label))
    pd.DataFrame(band_rows).to_csv(OUT / "pps_bands.csv", index=False)
    pd.concat(curves).to_csv(OUT / "pps_km.csv", index=False)

    # Restricted cubic spline of time to progression, reference 12 months.
    x = primary.pfs.clip(upper=primary.pfs.quantile(.99))
    basis = dmatrix("cr(x, df=4) - 1", {"x": x}, return_type="dataframe")
    basis.columns = [f"s{i}" for i in range(basis.shape[1])]
    # Drop one column because the natural spline basis sums to one within strata.
    sp = CoxPHFitter().fit(pd.concat([primary[BASE + ["pps", "death", "cohort"]], basis.iloc[:, 1:]], axis=1),
                           "pps", "death", strata=["cohort"])
    grid = np.linspace(2, 48, 93)
    design = np.asarray(build_design_matrices([basis.design_info], {"x": np.r_[12.0, grid]})[0])[:, 1:]
    contrast = design[1:] - design[0]
    names = list(basis.columns[1:])
    beta = sp.params_[names].to_numpy()
    cov = sp.variance_matrix_.loc[names, names].to_numpy()
    est = contrast @ beta
    se = np.sqrt(np.einsum("ij,jk,ik->i", contrast, cov, contrast))
    lr_s = 2 * (sp.log_likelihood_ - full.log_likelihood_)
    pd.DataFrame({"months": grid, "hr": np.exp(est), "lo": np.exp(est - 1.96 * se), "hi": np.exp(est + 1.96 * se),
                  "nonlinearity_p": stats.chi2.sf(lr_s, len(names) - 1)}).to_csv(OUT / "pps_spline.csv", index=False)

    # Leave-one-cohort-out concordance for post-progression survival with and without the timing term.
    conc = []
    rng = np.random.default_rng(SEED + 71)
    for held in PRIMARY:
        train = primary.loc[primary.cohort.ne(held)]
        test = primary.loc[primary.cohort.eq(held)].reset_index(drop=True)
        scores = []
        for cols in (BASE, BASE + ["ttp6"]):
            mod = CoxPHFitter().fit(train[cols + ["pps", "death", "cohort"]], "pps", "death", strata=["cohort"])
            scores.append(-(test[cols].to_numpy() @ mod.params_[cols].to_numpy()))
        point = [concordance_index(test.pps, s, test.death) for s in scores]
        boot = []
        for _ in range(500):
            i = rng.integers(0, len(test), len(test))
            c = [concordance_index(test.pps.to_numpy()[i], s[i], test.death.to_numpy()[i]) for s in scores]
            boot.append(c[1] - c[0])
        conc.append({"cohort": held, "n": len(test), "deaths": int(test.death.sum()), "c_baseline": point[0],
                     "c_with_timing": point[1], "delta": point[1] - point[0], "se": np.std(boot, ddof=1),
                     "lo": np.percentile(boot, 2.5), "hi": np.percentile(boot, 97.5)})
    conc = pd.DataFrame(conc)
    r = meta(conc.delta, conc.se)
    conc = pd.concat([conc, pd.DataFrame([{"cohort": "Random-effects summary", "n": conc.n.sum(),
                                           "deaths": conc.deaths.sum(), "delta": r["estimate"], "lo": r["lo"],
                                           "hi": r["hi"]}])])
    conc.to_csv(OUT / "pps_concordance.csv", index=False)
    print("A done", flush=True)


# ---------------------------------------------------------------- B. landmark supermodel, illness-death
def stacked(raw):
    parts = []
    for lm in LANDMARKS:
        frames, _ = cohort_data(raw, landmark=lm)
        for cohort, d in frames.items():
            parts.append(d.assign(lm=lm, cohort=cohort, pid=cohort + "_" + d.id.astype(str),
                                  event=(d.death.eq(1) & d.time.le(12)).astype(int), stop=d.time.clip(upper=12)))
    return pd.concat(parts, ignore_index=True)


def supermodel(raw):
    s = stacked(raw)
    s["early_lm"] = s.early * (s.lm - 12) / 12
    rows, dist = [], []
    cols = STATUS + ["early_lm", "timing6"]
    m = CoxPHFitter().fit(s[cols + ["stop", "event", "cohort", "lm", "pid"]], "stop", "event",
                          strata=["cohort", "lm"], cluster_col="pid")
    for term in cols:
        rows.append(hr_row(m, term, model="Supermodel, landmarks 6 to 24 months", landmark="all",
                           n=s.pid.nunique(), rows=len(s), deaths=int(s.event.sum())))
    late = s.loc[s.lm.ge(12)]
    m = CoxPHFitter().fit(late[cols + ["stop", "event", "cohort", "lm", "pid"]], "stop", "event",
                          strata=["cohort", "lm"], cluster_col="pid")
    rows.append(hr_row(m, "timing6", model="Supermodel, landmarks 12 to 24 months", landmark="12-24",
                       n=late.pid.nunique(), rows=len(late), deaths=int(late.event.sum())))
    for lm, d in s.groupby("lm"):
        m = CoxPHFitter().fit(d[TIMING + ["stop", "event", "cohort"]], "stop", "event", strata=["cohort"])
        for term in ("early", "timing6"):
            rows.append(hr_row(m, term, model="Single landmark", landmark=lm, n=len(d), rows=len(d),
                               deaths=int(d.event.sum())))
        for cohort, c in d.groupby("cohort"):
            p = c.loc[c.early.eq(1)]
            dist.append({"landmark": lm, "cohort": cohort, "n": len(c), "deaths12": int(c.event.sum()),
                         "progressed": len(p), "progressed_ge6": int(p.duration.ge(6).sum()),
                         "duration_median": p.duration.median(), "duration_q1": p.duration.quantile(.25),
                         "duration_q3": p.duration.quantile(.75), "duration_sd": p.duration.std()})
    pd.DataFrame(rows).to_csv(OUT / "supermodel_hr.csv", index=False)
    pd.DataFrame(dist).to_csv(OUT / "landmark_timing_distribution.csv", index=False)
    s.loc[s.early.eq(1), ["lm", "cohort", "duration"]].to_csv(OUT / "landmark_durations.csv", index=False)

    # Internal-external discrimination contrasts at every landmark.
    rng = np.random.default_rng(SEED + 19)
    out = []
    for lm in LANDMARKS:
        frames, _ = cohort_data(raw, landmark=lm)
        for held in PRIMARY:
            train = [d for c, d in frames.items() if c != held]
            test = frames[held]
            pred = {k: fit_pooled(train, cols).predict_proba(test[cols])[:, 1] for k, cols in SETS.items()}
            point = metrics(test, pred["Status"], 12)
            auc = {k: metrics(test, p, 12)["auc"] for k, p in pred.items()}
            t, e = test.time.to_numpy(), test.death.to_numpy()
            boot = []
            for _ in range(300):
                i = rng.integers(0, len(t), len(t))
                y, w, _, controls = horizon_arrays(t[i], e[i], 12)
                if not y.sum() or not controls:
                    continue
                a = {k: roc_auc_score(y, p[i], sample_weight=w) for k, p in pred.items()}
                boot.append([a["Status"] - a["Baseline"], a["Timing"] - a["Status"]])
            se = np.std(boot, axis=0, ddof=1)
            out.append({"landmark": lm, "cohort": held, "n": point["n"], "deaths": point["deaths"],
                        **{f"auc_{k}": v for k, v in auc.items()},
                        "status_delta": auc["Status"] - auc["Baseline"], "status_se": se[0],
                        "timing_delta": auc["Timing"] - auc["Status"], "timing_se": se[1]})
    out = pd.DataFrame(out)
    out.to_csv(OUT / "landmark_iecv.csv", index=False)
    pooled = []
    for lm, d in out.groupby("landmark"):
        for name in ("status", "timing"):
            r = meta(d[f"{name}_delta"], d[f"{name}_se"].clip(lower=1e-4))
            pooled.append({"landmark": lm, "contrast": name, "n": int(d.n.sum()), "deaths": int(d.deaths.sum()),
                           **{k: r[k] for k in ("estimate", "lo", "hi", "fixed", "fixed_lo", "fixed_hi", "i2")}})
    pd.DataFrame(pooled).to_csv(OUT / "landmark_iecv_pooled.csv", index=False)
    print("B supermodel done", flush=True)


def illness_death(raw):
    rows = []
    for cohort in PRIMARY:
        d = base_population(raw[cohort]).dropna(subset=["age10", "rd"])
        for r in d.itertuples():
            pid = f"{cohort}_{r.id}"
            base = {"pid": pid, "cohort": cohort, "age10": r.age10, "iv": r.iv, "rd": r.rd}
            progressed = r.progression_event == 1 and r.pfs < r.os - 1e-6
            first_stop = r.pfs if (progressed or r.progression_event == 0) else r.os
            first_stop = min(first_stop, r.os)
            died_without = (not progressed) and r.death == 1 and r.progression_event == 1
            if first_stop > 0:
                rows.append({**base, "start": 0.0, "stop": first_stop, "event": int(died_without), "prog": 0,
                             "prog_ttp": 0.0})
            if progressed:
                rows.append({**base, "start": r.pfs, "stop": r.os, "event": int(r.death), "prog": 1,
                             "prog_ttp": (r.pfs - 12) / 6})
    s = pd.DataFrame(rows)
    s = s.loc[s.stop.gt(s.start)]
    out = []
    for name, cols in (("Markov", BASE + ["prog"]), ("With time of progression", BASE + ["prog", "prog_ttp"])):
        m = CoxPHFitter().fit(s[cols + ["start", "stop", "event", "cohort", "pid"]], "stop", "event",
                              entry_col="start", strata=["cohort"], cluster_col="pid")
        for term in cols:
            out.append(hr_row(m, term, model=name, patients=s.pid.nunique(), deaths=int(s.event.sum()),
                              deaths_after_progression=int(s.loc[s.prog.eq(1), "event"].sum()),
                              loglik=m.log_likelihood_))
    pd.DataFrame(out).to_csv(OUT / "illness_death.csv", index=False)
    print("B illness-death done", flush=True)


# ---------------------------------------------------------------- C. IECV, recalibration, imputation, precision
def boot_eval(d, preds, reps, rng):
    t, e = d.time.to_numpy(), d.death.to_numpy()
    draws = []
    for _ in range(reps):
        i = rng.integers(0, len(t), len(t))
        y, w, _, controls = horizon_arrays(t[i], e[i], 12)
        if not y.sum() or not controls:
            continue
        auc = {k: roc_auc_score(y, p[i], sample_weight=w) for k, p in preds.items()}
        brier = {k: np.dot(w, (y - p[i]) ** 2) / len(i) for k, p in preds.items()}
        p = preds["Status"][i]
        lp = clip_logit(p)
        b = weighted_logistic(np.column_stack((np.ones(len(i)), lp)), y, w)
        draws.append({**{f"auc_{k}": v for k, v in auc.items()}, **{f"brier_{k}": v for k, v in brier.items()},
                      "oe": np.dot(w, y) / len(i) / p.mean(), "slope": b[1],
                      "d_status_baseline": auc["Status"] - auc["Baseline"],
                      "d_status_progonly": auc["Status"] - auc["Progression only"],
                      "d_timing_status": auc["Timing"] - auc["Status"],
                      "db_status_progonly": brier["Status"] - brier["Progression only"]})
    return pd.DataFrame(draws)


def smooth_curve(d, p, label, cohort):
    y, w, _, _ = horizon_data(d, 12)
    lp = clip_logit(p)
    basis = dmatrix("cr(x, df=3)", {"x": lp}, return_type="dataframe")
    beta = weighted_logistic(np.asarray(basis), y, w)
    fitted = expit(np.asarray(basis) @ beta)
    grid = np.linspace(np.quantile(lp, .02), np.quantile(lp, .98), 60)
    curve = expit(np.asarray(build_design_matrices([basis.design_info], {"x": grid})[0]) @ beta)
    ici = np.average(np.abs(fitted - p))
    return (pd.DataFrame({"model": label, "cohort": cohort, "predicted": expit(grid), "observed": curve}),
            {"model": label, "cohort": cohort, "ici": ici, "e90": np.quantile(np.abs(fitted - p), .9)})


def iecv(raw):
    frames, _ = cohort_data(raw)
    rng = np.random.default_rng(SEED + 5)
    cohort_rows, curves, icis, preds_all = [], [], [], {}
    designs = {
        "IECV": lambda held: [d for c, d in frames.items() if c != held],
        "TCGA only (locked)": lambda held: [frames["TCGA"]],
        "IECV without German training cohorts": lambda held: [d for c, d in frames.items()
                                                              if c != held and c not in ("UKE", "ICON7")],
    }
    for design, trainer in designs.items():
        for held in PRIMARY:
            if design != "IECV" and held == "TCGA":
                continue
            if design.endswith("cohorts") and held not in ("UKE", "ICON7"):
                continue
            train, test = trainer(held), frames[held]
            pred = {k: fit_pooled(train, cols).predict_proba(test[cols])[:, 1] for k, cols in SETS.items()}
            if design == "IECV":
                preds_all[held] = pred
            m = {k: metrics(test, p, 12) for k, p in pred.items()}
            b = boot_eval(test, pred, 1000, rng)
            row = {"design": design, "cohort": held, "n": m["Status"]["n"], "deaths": m["Status"]["deaths"]}
            point = {**{f"auc_{k}": v["auc"] for k, v in m.items()}, **{f"brier_{k}": v["brier"] for k, v in m.items()},
                     "oe": m["Status"]["oe"], "slope": m["Status"]["slope"],
                     "d_status_baseline": m["Status"]["auc"] - m["Baseline"]["auc"],
                     "d_status_progonly": m["Status"]["auc"] - m["Progression only"]["auc"],
                     "d_timing_status": m["Timing"]["auc"] - m["Status"]["auc"],
                     "db_status_progonly": m["Status"]["brier"] - m["Progression only"]["brier"]}
            for k, v in point.items():
                row.update({k: v, f"{k}_lo": b[k].quantile(.025), f"{k}_hi": b[k].quantile(.975),
                            f"{k}_se": b[k].std()})
            row["logit_auc_Status_se"] = clip_logit(b["auc_Status"]).std()
            row["logit_auc_Baseline_se"] = clip_logit(b["auc_Baseline"]).std()
            row["log_oe_se"] = np.log(b["oe"]).std()
            cohort_rows.append(row)
            curve, ici = smooth_curve(test, pred["Status"], design, held)
            curves.append(curve)
            icis.append(ici)
    table = pd.DataFrame(cohort_rows)
    table.to_csv(OUT / "iecv_cohort.csv", index=False)
    pd.concat(curves).to_csv(OUT / "calibration_curves.csv", index=False)
    pd.DataFrame(icis).to_csv(OUT / "calibration_ici.csv", index=False)

    pooled = []
    for design in ("IECV", "TCGA only (locked)"):
        d = table.loc[table.design.eq(design)]
        specs = {
            "AUC Baseline": (clip_logit(d.auc_Baseline), d.logit_auc_Baseline_se, expit),
            "AUC Status": (clip_logit(d.auc_Status), d.logit_auc_Status_se, expit),
            "O/E Status": (np.log(d.oe), d.log_oe_se, np.exp),
            "Calibration slope Status": (d.slope, d.slope_se, lambda v: v),
            "ΔAUC Status − Baseline": (d.d_status_baseline, d.d_status_baseline_se, lambda v: v),
            "ΔAUC Status − Progression only": (d.d_status_progonly, d.d_status_progonly_se, lambda v: v),
            "ΔAUC Timing − Status": (d.d_timing_status, d.d_timing_status_se, lambda v: v),
            "ΔBrier Status − Progression only": (d.db_status_progonly, d.db_status_progonly_se, lambda v: v),
        }
        for name, (est, se, back) in specs.items():
            r = meta(est, se)
            pooled.append({"design": design, "measure": name, "k": r["k"], "estimate": back(r["estimate"]),
                           "lo": back(r["lo"]), "hi": back(r["hi"]), "pi_lo": back(r["pi_lo"]),
                           "pi_hi": back(r["pi_hi"]), "i2": r["i2"]})
    pd.DataFrame(pooled).to_csv(OUT / "iecv_pooled.csv", index=False)

    # Coefficients: locked TCGA equation, each IECV fit and the final four-cohort equation with bootstrap CIs.
    coef = []
    fits = {"TCGA only": [frames["TCGA"]], "All four cohorts": list(frames.values()),
            **{f"Without {c}": [d for k, d in frames.items() if k != c] for c in PRIMARY}}
    for name, train in fits.items():
        for label in ("Status", "Timing"):
            m = fit_pooled(train, SETS[label])
            coef.append({"fit": name, "model": label, "intercept": m.intercept_[0],
                         **dict(zip(SETS[label], m.coef_[0], strict=True))})
    boot = []
    for _ in range(500):
        sample = [d.iloc[rng.integers(0, len(d), len(d))] for d in frames.values()]
        for label in ("Status", "Timing"):
            m = fit_pooled(sample, SETS[label])
            boot.append({"model": label, "intercept": m.intercept_[0], **dict(zip(SETS[label], m.coef_[0], strict=True))})
    boot = pd.DataFrame(boot)
    for label, b in boot.groupby("model"):
        coef.append({"fit": "All four cohorts, 2.5%", "model": label, **b.drop(columns="model").quantile(.025).to_dict()})
        coef.append({"fit": "All four cohorts, 97.5%", "model": label, **b.drop(columns="model").quantile(.975).to_dict()})
    pd.DataFrame(coef).to_csv(OUT / "iecv_coefficients.csv", index=False)
    print("C IECV done", flush=True)
    return frames, preds_all


def recalibration(raw, frames):
    """Cross-validated updating of the locked TCGA Status equation inside each external cohort."""
    locked = fit_pooled([frames["TCGA"]], STATUS)
    rows = []
    for cohort in PRIMARY[1:]:
        d = frames[cohort].reset_index(drop=True)
        p = locked.predict_proba(d[STATUS])[:, 1]
        lp = clip_logit(p)
        y, w, _, _ = horizon_data(d, 12)
        updated = {"Intercept update": [], "Logistic recalibration": []}
        for rep in range(20):
            cv = StratifiedKFold(5, shuffle=True, random_state=SEED + rep)
            pa, pb = np.empty(len(d)), np.empty(len(d))
            for learn, test in cv.split(d, y + 2 * (w == 0)):
                yl, wl, _, _ = horizon_data(d.iloc[learn], 12)
                # Intercept with the slope held at one.
                a = brentq(lambda v: np.dot(wl, expit(lp[learn] + v) - yl), -10, 10)
                b = weighted_logistic(np.column_stack((np.ones(len(learn)), lp[learn])), yl, wl)
                pa[test] = expit(lp[test] + a)
                pb[test] = expit(b[0] + b[1] * lp[test])
            updated["Intercept update"].append(metrics(d, pa, 12))
            updated["Logistic recalibration"].append(metrics(d, pb, 12))
        rows.append({"cohort": cohort, "method": "No update", **{k: metrics(d, p, 12)[k]
                                                                    for k in ("oe", "slope", "brier", "auc", "intercept")}})
        for method, values in updated.items():
            v = pd.DataFrame(values)
            rows.append({"cohort": cohort, "method": method, **v[["oe", "slope", "brier", "auc", "intercept"]].mean().to_dict()})
    pd.DataFrame(rows).to_csv(OUT / "recalibration.csv", index=False)
    print("C recalibration done", flush=True)


def imputation(raw):
    """Multiple imputation of missing residual-disease status inside each landmark cohort."""
    frames, _ = cohort_data(raw, include_rd=False)
    rng = np.random.default_rng(SEED + 77)
    missing = {c: int(d.rd.isna().sum()) for c, d in frames.items()}
    draws, coefs = [], []
    for m in range(30):
        filled = {}
        for cohort, d in frames.items():
            d = d.copy()
            na = d.rd.isna()
            if na.any():
                y, _, _, _ = horizon_data(d, 12)
                x = np.column_stack((d.age10, d.iv, d.early, y, np.log(d.time.clip(lower=.1))))
                obs = np.flatnonzero(~na.to_numpy())
                i = rng.choice(obs, len(obs))
                model = LogisticRegression(C=10, max_iter=1000).fit(x[i], d.rd.to_numpy()[i])
                d.loc[na, "rd"] = rng.binomial(1, model.predict_proba(x[na.to_numpy()])[:, 1])
            filled[cohort] = d
        fits = {k: fit_pooled([filled["TCGA"]], cols) for k, cols in SETS.items()}
        coefs.append({"intercept": fits["Status"].intercept_[0], **dict(zip(STATUS, fits["Status"].coef_[0], strict=True))})
        for cohort in PRIMARY[1:]:
            d = filled[cohort]
            pred = {k: f.predict_proba(d[SETS[k]])[:, 1] for k, f in fits.items()}
            ms = {k: metrics(d, p, 12) for k, p in pred.items()}
            draws.append({"cohort": cohort, "n": len(d), "deaths": ms["Status"]["deaths"],
                          "auc_baseline": ms["Baseline"]["auc"], "auc_status": ms["Status"]["auc"],
                          "d_status_baseline": ms["Status"]["auc"] - ms["Baseline"]["auc"],
                          "d_timing_status": ms["Timing"]["auc"] - ms["Status"]["auc"], "oe": ms["Status"]["oe"]})
    draws = pd.DataFrame(draws)
    summary = draws.groupby("cohort").agg(["mean", "min", "max"])
    summary.columns = ["_".join(c) for c in summary.columns]
    summary = summary.reset_index()
    summary["missing_rd"] = summary.cohort.map(missing)
    summary["n_tcga"] = len(frames["TCGA"])
    summary["missing_rd_tcga"] = missing["TCGA"]
    summary.to_csv(OUT / "imputation.csv", index=False)
    pd.DataFrame(coefs).agg(["mean", "min", "max"]).to_csv(OUT / "imputation_coefficients.csv")
    print("C imputation done", flush=True)


def precision(frames, preds):
    """Riley-type precision targets for an external validation of a binary 12-month outcome."""
    rows = []
    lp = np.concatenate([clip_logit(preds[c]["Status"]) for c in PRIMARY])
    p = expit(lp)
    ia, ib, iab = np.mean(p * (1 - p)), np.mean(lp ** 2 * p * (1 - p)), np.mean(lp * p * (1 - p))
    phi = float(np.mean([metrics(frames[c], preds[c]["Status"], 12)["observed"] for c in PRIMARY]))
    n_oe = (1 - phi) / (phi * (0.2 / 1.96) ** 2)          # 95% CI for O/E of about 0.8 to 1.2
    n_slope = ia / ((0.1 / 1.96) ** 2 * (ia * ib - iab ** 2))  # 95% CI half-width of 0.1 for the slope

    def se_c(c, n):
        q1, q2 = c / (2 - c), 2 * c ** 2 / (1 + c)
        n1, n0 = n * phi, n * (1 - phi)
        return np.sqrt((c * (1 - c) + (n1 - 1) * (q1 - c ** 2) + (n0 - 1) * (q2 - c ** 2)) / (n1 * n0))

    n_c = next(n for n in range(50, 20000) if se_c(.8, n) <= 0.05 / 1.96)
    for name, n in (("O/E within about 0.8 to 1.2", n_oe), ("Calibration slope ±0.1", n_slope),
                    ("C statistic 0.80 ± 0.05", n_c)):
        rows.append({"target": name, "outcome_proportion": phi, "required_n": int(np.ceil(n)),
                     "required_events": int(np.ceil(n * phi))})
    for c in PRIMARY[1:]:
        m = metrics(frames[c], preds[c]["Status"], 12)
        rows.append({"target": f"Available in {c}", "outcome_proportion": m["observed"], "required_n": m["n"],
                     "required_events": m["deaths"]})
    pd.DataFrame(rows).to_csv(OUT / "precision.csv", index=False)


# ---------------------------------------------------------------- E. period hazards and enrichment arithmetic
def period_hazards(frames):
    rows = []
    for cohort, d in frames.items():
        for name, lo, hi in (("0-12", 0, 12), ("12-24", 12, 24)):
            part = d.loc[d.time.gt(lo)].copy()
            part["event"] = (part.death.eq(1) & part.time.le(hi)).astype(int)
            part["stop"] = part.time.clip(upper=hi)
            part["start"] = float(lo)
            try:
                m = CoxPHFitter().fit(part[STATUS + ["start", "stop", "event"]], "stop", "event", entry_col="start")
                rows.append(hr_row(m, "early", cohort=cohort, period=name, n=len(part), deaths=int(part.event.sum())))
            except Exception as exc:  # sparse strata
                rows.append({"cohort": cohort, "period": name, "n": len(part), "deaths": int(part.event.sum()),
                             "note": type(exc).__name__})
    pd.DataFrame(rows).to_csv(OUT / "period_hazards.csv", index=False)


def enrichment(frames):
    events = 4 * (stats.norm.ppf(.975) + stats.norm.ppf(.8)) ** 2 / np.log(.7) ** 2
    rows = []
    pooled = pd.concat([d.assign(cohort=c) for c, d in frames.items()])
    for label, d in [("All landmark survivors", pooled), ("Progressed by month 12", pooled.loc[pooled.early.eq(1)]),
                     ("Not progressed by month 12", pooled.loc[pooled.early.eq(0)])]:
        risks = []
        for _, c in d.groupby("cohort"):
            risks.append(1 - KaplanMeierFitter().fit(c.time, c.death).predict(12))
        risk = float(np.average(risks, weights=d.groupby("cohort").size()))
        rows.append({"population": label, "n": len(d), "mortality12": risk, "mortality12_min": min(risks),
                     "mortality12_max": max(risks), "events_required": int(np.ceil(events)),
                     "patients_required": int(np.ceil(events / risk))})
    pd.DataFrame(rows).to_csv(OUT / "enrichment.csv", index=False)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    sources = all_sources()
    if "--skip-pps" not in sys.argv:
        pps_analysis(sources)
    supermodel(sources)
    illness_death(sources)
    landmark_frames, iecv_predictions = iecv(sources)
    recalibration(sources, landmark_frames)
    imputation(sources)
    precision(landmark_frames, iecv_predictions)
    period_hazards(landmark_frames)
    enrichment(landmark_frames)
    print("clinical analyses complete")

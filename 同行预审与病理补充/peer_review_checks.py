"""Run exploratory checks added after the original manuscript was written."""
from pathlib import Path
import json
import sys

from lifelines import KaplanMeierFitter
import numpy as np
import pandas as pd
from scipy.special import expit
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "scripts"))

from analyze import BASE, SEED, STATUS, cohort_data, cox, fit_predict, horizon_arrays, horizon_data, log_model, metrics
from prepare_data import MONTH, common, eligible, load_all

OUT = HERE / "results"
DATA = HERE / "data"
SOURCE_MODEL = json.loads((ROOT / "results" / "equations.json").read_text(encoding="utf-8"))


def discrimination(y, w, p):
    return roc_auc_score(y, p, sample_weight=w)


def bootstrap_metrics(d, scores, repetitions=1000, seed=SEED):
    """Paired evaluation-set bootstrap; all predictions remain fixed."""
    rng = np.random.default_rng(seed)
    time, event = d.time.to_numpy(), d.death.to_numpy()

    def values(idx):
        y, w, _, n_controls = horizon_arrays(time[idx], event[idx], 12)
        if not y.sum() or not n_controls:
            return None
        first, second = (np.asarray(p)[idx] for p in scores)
        prev, new = discrimination(y, w, first), discrimination(y, w, second)
        brier_delta = np.dot(w, (y - second) ** 2 - (y - first) ** 2) / len(idx)
        oe = np.dot(w, y) / len(idx) / second.mean()
        return [prev, new, new - prev, brier_delta, oe]

    point = values(np.arange(len(d)))
    results = [v for _ in range(repetitions) if (v := values(rng.integers(0, len(d), len(d)))) is not None]
    lo, hi = np.percentile(results, [2.5, 97.5], axis=0)
    labels = ["auc_reference", "auc_added", "delta_auc", "delta_brier", "oe_added"]
    estimates = {}
    for name, est, low, high in zip(labels, point, lo, hi, strict=True):
        estimates.update({name: est, f"{name}_lo": low, f"{name}_hi": high})
    return {"n": len(d), "deaths12": int(((d.time <= 12) & d.death.eq(1)).sum()),
            "bootstrap_valid": len(results), **estimates}


def predict_locked(d, name):
    terms = SOURCE_MODEL[name]
    return expit(terms["intercept"] + sum(d[k].to_numpy() * b for k, b in terms.items() if k != "intercept"))


def original_checks(frames):
    primary = []
    progressed = []
    conditional = []
    highgrade = []
    status_only = fit_predict(log_model(), frames["TCGA"], frames, ["early"], 12)
    for order, (name, d) in enumerate(frames.items()):
        if name == "TCGA":
            continue
        early = d.early.to_numpy().astype(bool)
        status = predict_locked(d, "Status")
        timing = predict_locked(d, "Timing")
        primary.append({"cohort": name, **bootstrap_metrics(d, (status_only[name], status), seed=SEED + order)})
        subset = d.loc[early].reset_index(drop=True)
        progressed.append({"cohort": name, **bootstrap_metrics(subset, (status[early], timing[early]), seed=SEED + order)})
        fields = [np.zeros(len(d)), np.zeros(len(d))]
        for group in (False, True):
            rows = early == group
            y, weights, _, _ = horizon_data(d.loc[rows], 12)
            fields[0][rows], fields[1][rows] = y, weights
        conditional.append({"cohort": name, "status_auc": discrimination(*fields, status),
                            "timing_auc": discrimination(*fields, timing)})
        if name in ("AOCS", "ICON7"):
            grade = d.grade.astype(str).str.lower()
            known = grade.isin(["3", "3.0"]) if name == "AOCS" else grade.eq("high.grade")
            sub = d.loc[known].reset_index(drop=True)
            highgrade += [
                {"cohort": name, "contrast": "Status vs Baseline",
                 **bootstrap_metrics(sub, (predict_locked(sub, "Baseline"), status[known]), seed=SEED+order)},
                {"cohort": name, "contrast": "Timing vs Status",
                 **bootstrap_metrics(sub, (status[known], timing[known]), seed=SEED + order)},
            ]
    for name, rows in (("clinical_benchmark", primary), ("progressed_only", progressed),
                       ("conditional_ipcw", conditional), ("source_high_grade", highgrade)):
        pd.DataFrame(rows).to_csv(OUT / f"{name}.csv", index=False)


def tcga_stroma():
    cases = json.loads((DATA / "TCGA_OV_GDC_pathology.json").read_text(encoding="utf-8"))["data"]["hits"]
    records = []
    for case in cases:
        for sample in case.get("samples", []):
            if sample.get("sample_type") != "Primary Tumor":
                continue
            for portion in sample.get("portions", []):
                for slide in portion.get("slides", []):
                    if slide.get("percent_stromal_cells") is not None:
                        records.append({"id": case["submitter_id"], "slide": slide["submitter_id"],
                                        "stroma": slide["percent_stromal_cells"]})
    s = pd.DataFrame(records).dropna(subset=["stroma"]).drop_duplicates()
    assert s.slide.is_unique and s.stroma.between(0, 100).all()
    return s.groupby("id", as_index=False).agg(stroma=("stroma", "median"), n_slides=("slide", "size"))


def apollo():
    source = pd.read_excel(DATA / "APOLLO2_MOESM5.xlsx", sheet_name="2A")
    assert source.apollo_id.is_unique
    d = common("APOLLO-2", source.apollo_id, source.age_at_diagnosis,
               source.tumor_stage.replace({"IVA": "IV", "IVB": "IV"}),
               source.residualdzcat.map({"R0 (none or microscopic)": 0, "R1 (macroscopic >1 cm)": 1}),
               source.survivaltimeindays / MONTH, source.survival_status.map({"Dead": 1, "Alive": 0}),
               source.pfstimeindays / MONTH, source.pfs_status.map({"Event": 1, "Censored": 0}),
               pd.Series(True, index=source.index), pd.Series("Unavailable", index=source.index))
    d["stroma"] = source.PATH_stroma
    d["primary_sample"] = source.classification_of_tumor.eq("Primary")
    d["primary_debulking"] = source.typeofsurgery.eq("Primary debulking surgery")
    d["parp"] = source.ACT_PARP.eq("YES") | source.MACT_PARP.eq("YES")
    assert d.stroma.between(0, 100).all()
    assert d.pfs.le(d.os).all()
    d, _ = eligible(d)
    d["age10"] = (d.age - 60) / 10
    d["timing6"] = np.where(d.early == 1, (d.duration - 6) / 6, 0)
    return d


def oof(d, columns):
    y, w, _, _ = horizon_data(d, 12)
    stratify = y + 2 * (w == 0)
    split = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    p = np.empty(len(d))
    for train, test in split.split(d, stratify):
        p[test] = fit_predict(log_model(), d.iloc[train], {"test": d.iloc[test]}, columns, 12)["test"]
    return p


def pathology_checks(frames):
    raw_path = tcga_stroma()
    data = frames["TCGA"].merge(raw_path, how="left", on="id", validate="one_to_one")
    selected = data.loc[data.stroma.notna()].copy()
    selected["stroma10"] = selected.stroma / 10
    p0 = oof(selected, STATUS)
    p1 = oof(selected, STATUS + ["stroma10"])
    pd.DataFrame([bootstrap_metrics(selected, (p0, p1))]).to_csv(OUT / "tcga_stroma_oof.csv", index=False)

    external = apollo()
    report, tests = [], []
    for label, d in (("All sites", external), ("Primary-site specimens", external.loc[external.primary_sample]),
                     ("Primary debulking", external.loc[external.primary_debulking])):
        d = d.reset_index(drop=True)
        predictions = {m: predict_locked(d, m) for m in ("Baseline", "Status", "Timing")}
        for model, score in predictions.items():
            m = metrics(d, score, 12)
            report.append({"subset": label, "model": model, "n_progressed": int(d.early.sum()), **m})
        a, b, c = (predictions[k] for k in ("Baseline", "Status", "Timing"))
        for contrast, scores in (("Status vs Baseline", (a, b)), ("Timing vs Status", (b, c))):
            tests.append({"subset": label, "contrast": contrast, **bootstrap_metrics(d, scores)})
    pd.DataFrame(report).to_csv(OUT / "apollo_clinical_transport.csv", index=False)
    pd.DataFrame(tests).to_csv(OUT / "apollo_clinical_contrasts.csv", index=False)

    counts = []
    for name, d in (("TCGA all", data), ("TCGA with stroma", selected),
                    ("TCGA without stroma", data.loc[data.stroma.isna()]), ("APOLLO all sites", external),
                    ("APOLLO primary sites", external.loc[external.primary_sample])):
        y, _, _, _ = horizon_data(d, 12)
        counts.append({"subset": name, "n": len(d), "deaths12": y.sum(), "early": d.early.sum(),
                       "stroma_median": d.stroma.median() if d.stroma.notna().any() else np.nan,
                       "stroma_q1": d.stroma.quantile(.25), "stroma_q3": d.stroma.quantile(.75),
                       "median_slides": d.n_slides.median() if "n_slides" in d and d.n_slides.notna().any() else np.nan})
    pd.DataFrame(counts).to_csv(OUT / "pathology_coverage.csv", index=False)

    # Fit the pathology-adjusted association only where the event count supports it.
    model, model_data = cox(selected, STATUS + ["stroma10"])
    model.summary.assign(events24=int(model_data.death.sum())).to_csv(OUT / "tcga_stroma_cox.csv")

    rows = []
    for cohort, d in [("TCGA", selected), ("APOLLO primary sites", external.loc[external.primary_sample])]:
        for event, p in d.groupby("early"):
            km = KaplanMeierFitter().fit(p.time, p.death)
            ci = km.confidence_interval_.loc[:12].iloc[-1]
            rows.append({"cohort": cohort, "early": event, "n": len(p), "risk12": 1 - km.predict(12),
                         "lower": 1 - ci.iloc[1], "upper": 1 - ci.iloc[0]})
    pd.DataFrame(rows).to_csv(OUT / "pathology_subset_landmark_risk.csv", index=False)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    cohorts, _ = cohort_data(load_all())
    original_checks(cohorts)
    pathology_checks(cohorts)

"""Post-review checks using the original study definitions and locked models."""
from pathlib import Path
import json
import sys
import warnings

import numpy as np
import pandas as pd
import rdata
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA, OUT = HERE / "data", HERE / "results"
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "同行预审与病理补充"))

from analyze import (BASE, SEED, STATUS, TIMING, calibration, cohort_data,
                     fit_predict, horizon_arrays, horizon_data, log_model, metrics)
from prepare_data import MONTH, common, eligible, load_all
from peer_review_checks import bootstrap_metrics, predict_locked


def safe_phenotype(name):
    p = DATA / (name + "_eset.rda")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        s = next(iter(rdata.read_rda(p).values())).phenoData.data
    s.columns = s.columns.map(str)
    s = s.reset_index(names="accession")
    # Retain analytical fields only. Original sample metadata include calendar dates.
    return s[["accession", "histological_type", "tumorstage", "grade",
              "age_at_initial_pathologic_diagnosis", "debulking",
              "days_to_death", "vital_status", "days_to_tumor_recurrence",
              "recurrence_status", "sample_type"]]


def curated_cohort(name, label, required):
    s = safe_phenotype(name)
    st = s.tumorstage.map({1: "I", 2: "II", 3: "III", 4: "IV"})
    rd = s.debulking.map({"optimal": 0, "suboptimal": 1})
    if "rd" not in required:
        rd = pd.Series(np.nan, index=s.index)
    d = common(
        label, s.accession, s.age_at_initial_pathologic_diagnosis, st, rd,
        s.days_to_death / MONTH, s.vital_status.map({"deceased": 1, "living": 0}),
        s.days_to_tumor_recurrence / MONTH, s.recurrence_status.map({"recurrence": 1, "norecurrence": 0}),
        s.histological_type.eq("ser") & s.sample_type.eq("tumor"),
        pd.Series("Unavailable", index=s.index), s.grade,
    )
    d["pf_gt_os"] = d.pfs.gt(d.os + 2 / MONTH)
    invalid = int(d.loc[d.advanced & d.serous, "pf_gt_os"].sum())
    d = d.loc[~d.pf_gt_os].copy()
    data, steps = eligible(d, require=required)
    data["age10"] = (data.age - 60) / 10
    data["timing6"] = np.where(data.early == 1, (data.duration - 6) / 6, 0)
    export_columns = ["id", "stage", "age", "rd", "os", "death", "pfs", "progression_event", "early", "time"]
    data[export_columns].to_csv(OUT / (name + "_landmark.csv"), index=False)
    return data, steps, invalid, len(s)


def persist_extended_cohorts(dev):
    rows, contrasts, flow = [], [], []
    specs = [
        ("GSE30161", "UVA", ("age", "rd"), BASE, "Full"),
        ("GSE26193", "Curie", (), ["iv"], "Stage only baseline"),
        ("GSE49997", "OVCAD", ("age",), ["age10", "iv"], "No residual predictor"),
    ]
    for i, (accession, label, required, variables, scope) in enumerate(specs):
        d, steps, invalid, n_source = curated_cohort(accession, label, required)
        fit = {
            "Baseline": variables, "Progression only": ["early"],
            "Status": variables + ["early"], "Timing": variables + ["early", "timing6"],
        }
        probs = {}
        for name, cols in fit.items():
            if scope == "Full" and name != "Progression only":
                probs[name] = predict_locked(d, name)
            else:
                probs[name] = fit_predict(log_model(), dev, {label: d}, cols, 12)[label]
            rows.append({"cohort": label, "accession": accession, "scope": scope, "model": name,
                         "progressed": int(d.early.sum()), **metrics(d, probs[name], 12)})
        for n, (first, second) in enumerate(
            (("Baseline", "Status"), ("Progression only", "Status"), ("Status", "Timing"))
        ):
            c = bootstrap_metrics(d, (probs[first], probs[second]), repetitions=1000, seed=SEED + 3*i + n)
            contrasts.append({"cohort": label, "contrast": f"{second} vs {first}", **c})
        flow.append({"cohort": label, "source_n": n_source, "late_serous_inconsistent_pfs_os": invalid, **steps})
    pd.DataFrame(rows).to_csv(OUT / "additional_validation.csv", index=False)
    pd.DataFrame(contrasts).to_csv(OUT / "additional_contrasts.csv", index=False)
    pd.DataFrame(flow).to_csv(OUT / "additional_flow.csv", index=False)


def link_eligible_slide_files(dev):
    qc = ROOT / "同行预审与病理补充/data/TCGA_OV_GDC_pathology.json"
    cases = json.loads(qc.read_text(encoding="utf-8"))["data"]["hits"]
    ids = set(dev.id)
    records = []
    for case in cases:
        if case["submitter_id"] not in ids:
            continue
        for sample in case.get("samples", []):
            if sample.get("sample_type") != "Primary Tumor":
                continue
            for portion in sample.get("portions", []):
                for slide in portion.get("slides", []):
                    stroma = slide.get("percent_stromal_cells")
                    if stroma is not None:
                        records.append((slide["submitter_id"], case["submitter_id"], stroma))
    scores = pd.DataFrame(records, columns=["slide", "id", "stroma"]).drop_duplicates()
    release = json.loads((DATA / "TCGA_OV_slide_manifest.json").read_text(encoding="utf-8"))["data"]["hits"]
    files = [(entity["entity_submitter_id"], hit["file_id"], hit["file_size"]/2**20, hit["file_name"])
             for hit in release for entity in hit.get("associated_entities", [])]
    files = pd.DataFrame(files, columns=["slide", "file_id", "mb", "filename"]).drop_duplicates()
    merged = scores.merge(files, how="inner", on="slide").sort_values(["mb", "slide"])
    merged.to_csv(DATA / "eligible_qc_slide_records.csv", index=False)


def rcs(d):
    x = d.duration / 6
    curvature = (np.maximum(x, 0) ** 3 - 2 * np.maximum(x - 1, 0) ** 3 + np.maximum(x - 2, 0) ** 3) / 4
    return d.assign(curvature=curvature * d.early)


def percentile(a):
    a = np.asarray(a)
    return tuple(np.nanquantile(a, [0.025, 0.975]))


def two_stage_and_calibration(original):
    frames = {k: rcs(v) for k, v in original.items()}
    dev = frames.pop("TCGA")
    modes = {"Status": STATUS, "Timing": TIMING, "Timing spline": TIMING + ["curvature"]}
    pred = {name: fit_predict(log_model(), dev, frames, cols, 12) for name, cols in modes.items()}
    indexes = {k: (d.time.to_numpy(), d.death.to_numpy()) for k, d in frames.items()}
    limits = {k: {"calibration": [], "increment": []} for k in frames}
    rng = np.random.default_rng(SEED + 331)
    for j in range(600):
        training = dev.iloc[rng.integers(0, len(dev), len(dev))]
        # A shared bootstrap development sample preserves paired model comparisons.
        refit = {name: fit_predict(log_model(), training, frames, cols, 12) for name, cols in modes.items()}
        for k, d in frames.items():
            idx = rng.integers(0, len(d), len(d))
            tt, ee = indexes[k]
            y, w, _, controls = horizon_arrays(tt[idx], ee[idx], 12)
            if not y.sum() or not controls:
                continue
            fixed = pred["Status"][k][idx]
            ic, sl = calibration(y, w, fixed)
            ratio = np.dot(w, y) / len(d) / fixed.mean()
            limits[k]["calibration"].append([ic, sl, ratio])
            evaluations = [roc_auc_score(y, refit[name][k][idx], sample_weight=w) for name in modes]
            limits[k]["increment"].append([evaluations[1]-evaluations[0], evaluations[2]-evaluations[0]])
    out = []
    for k, d in frames.items():
        m = metrics(d, pred["Status"][k], 12)
        linear = metrics(d, pred["Timing"][k], 12)
        spline = metrics(d, pred["Timing spline"][k], 12)
        vals = {
            "calibration_intercept": (m["intercept"], np.asarray(limits[k]["calibration"])[:, 0]),
            "calibration_slope": (m["slope"], np.asarray(limits[k]["calibration"])[:, 1]),
            "oe": (m["oe"], np.asarray(limits[k]["calibration"])[:, 2]),
            "timing_delta_auc_two_stage": (linear["auc"] - m["auc"], np.asarray(limits[k]["increment"])[:, 0]),
            "spline_delta_auc_two_stage": (spline["auc"] - m["auc"], np.asarray(limits[k]["increment"])[:, 1]),
        }
        for name, (value, dist) in vals.items():
            lo, hi = percentile(dist)
            out.append({"cohort": k, "estimate": name, "value": value, "lo": lo, "hi": hi,
                        "n": m["n"], "events": m["deaths"], "resamples": len(dist)})
    pd.DataFrame(out).to_csv(OUT / "calibration_and_two_stage.csv", index=False)


def adjacent_landmarks(raw):
    rows = []
    for landmark in (11, 13):
        cohorts, _ = cohort_data(raw, landmark=landmark)
        dev = cohorts.pop("TCGA")
        results = [fit_predict(log_model(), dev, cohorts, cols, 12) for cols in (BASE, STATUS, TIMING)]
        for j, (key, data) in enumerate(cohorts.items()):
            for k, (first, second) in enumerate(((0, 1), (1, 2))):
                estimates = bootstrap_metrics(data, (results[first][key], results[second][key]),
                                              repetitions=500, seed=SEED+landmark+j+k)
                rows.append({"cohort": key, "landmark": landmark,
                             "contrast": ("Status vs Baseline", "Timing vs Status")[k], **estimates})
    pd.DataFrame(rows).to_csv(OUT / "adjacent_landmarks.csv", index=False)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    sources = load_all()
    original, _ = cohort_data(sources)
    link_eligible_slide_files(original["TCGA"])
    persist_extended_cohorts(original["TCGA"])
    print("Additional cohorts evaluated.", flush=True)
    adjacent_landmarks(sources)
    two_stage_and_calibration(original)

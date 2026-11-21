"""Shared loaders and estimators for the second-revision analyses."""
from pathlib import Path
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize
from scipy.special import expit, logit

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA, OUT, FIG = HERE / "data", HERE / "results", HERE / "figures"
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "同行预审与病理补充"))
sys.path.insert(0, str(ROOT / "修订稿"))

from analyze import (BASE, SEED, STATUS, TIMING, cohort_data, horizon_arrays,  # noqa: E402
                     horizon_data, log_model, metrics)
from prepare_data import MONTH, common, load_all  # noqa: E402
from revision_analysis import safe_phenotype  # noqa: E402

PRIMARY = ("TCGA", "AOCS", "UKE", "ICON7")
EXTRA = {"UVA": "GSE30161", "Curie": "GSE26193", "OVCAD": "GSE49997"}
SETS = {"Baseline": BASE, "Progression only": ["early"], "Status": STATUS, "Timing": TIMING}


def extra_frame(label):
    s = safe_phenotype(EXTRA[label])
    d = common(
        label, s.accession, s.age_at_initial_pathologic_diagnosis,
        s.tumorstage.map({1: "I", 2: "II", 3: "III", 4: "IV"}),
        s.debulking.map({"optimal": 0, "suboptimal": 1}),
        s.days_to_death / MONTH, s.vital_status.map({"deceased": 1, "living": 0}),
        s.days_to_tumor_recurrence / MONTH, s.recurrence_status.map({"recurrence": 1, "norecurrence": 0}),
        s.histological_type.eq("ser") & s.sample_type.eq("tumor"),
        pd.Series("Unavailable", index=s.index), s.grade)
    return d.loc[~d.pfs.gt(d.os + 2 / MONTH)].copy()


def all_sources():
    """Primary sources, Japan and the three additional curated cohorts before any landmark restriction."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        raw = load_all()
        for label in EXTRA:
            raw[label] = extra_frame(label)
    return raw


def base_population(frame):
    keep = (frame.advanced & frame.serous & ~frame.treatment.eq("Other") & frame.os.gt(0) & frame.pfs.ge(0)
            & frame.death.isin([0, 1]) & frame.progression_event.isin([0, 1]))
    d = frame.loc[keep].copy()
    d["age10"] = (d.age - 60) / 10
    return d.reset_index(drop=True)


def fit_pooled(frames, columns, horizon=12):
    """Censoring-weighted logistic model with cohort-specific censoring weights."""
    x, y, w = [], [], []
    for d in frames:
        yy, ww, _, _ = horizon_data(d, horizon)
        keep = ww > 0
        x.append(d.loc[keep, columns])
        y.append(yy[keep])
        w.append(ww[keep])
    return log_model().fit(pd.concat(x), np.concatenate(y), sample_weight=np.concatenate(w))


def weighted_logistic(design, y, w):
    def loss(beta):
        eta = design @ beta
        return np.dot(w, np.logaddexp(0, eta) - y * eta), design.T @ (w * (expit(eta) - y))
    return minimize(loss, np.zeros(design.shape[1]), jac=True, method="BFGS").x


def clip_logit(p):
    return logit(np.clip(p, 1e-6, 1 - 1e-6))


def meta(est, se):
    """DerSimonian-Laird random effects with Hartung-Knapp interval and prediction interval."""
    est, se = np.asarray(est, float), np.asarray(se, float)
    k = len(est)
    w = 1 / se ** 2
    fixed = np.sum(w * est) / w.sum()
    q = np.sum(w * (est - fixed) ** 2)
    tau2 = max(0.0, (q - (k - 1)) / (w.sum() - np.sum(w ** 2) / w.sum()))
    ws = 1 / (se ** 2 + tau2)
    mu = np.sum(ws * est) / ws.sum()
    var = max(np.sum(ws * (est - mu) ** 2) / ((k - 1) * ws.sum()), 1 / ws.sum())
    t = stats.t.ppf(.975, k - 1)
    out = {"k": k, "estimate": mu, "lo": mu - t * var ** .5, "hi": mu + t * var ** .5,
           "tau2": tau2, "i2": max(0.0, (q - (k - 1)) / q) if q > 0 else 0.0,
           "fixed": fixed, "fixed_lo": fixed - 1.96 / w.sum() ** .5, "fixed_hi": fixed + 1.96 / w.sum() ** .5}
    if k > 2:
        tp = stats.t.ppf(.975, k - 2)
        out.update(pi_lo=mu - tp * (tau2 + var) ** .5, pi_hi=mu + tp * (tau2 + var) ** .5)
    return out

"""Second-revision additions: HRD score, CCNE1 amplification and other TCGA genomic features,
stroma-adjusted transcriptome comparison, and single-cell attribution of the early-progression programme."""
import json
import warnings

import gseapy
from lifelines import CoxPHFitter
import numpy as np
import pandas as pd
import requests
from scipy import stats
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

from v2_common import (BASE, DATA, OUT, PRIMARY, SEED, STATUS, all_sources, base_population, cohort_data,
                       fit_pooled, horizon_arrays, horizon_data, metrics)
from v2_molecular import expression, with_scores

warnings.filterwarnings("ignore")
CBIO = "https://www.cbioportal.org/api"
STUDY = "ov_tcga_pan_can_atlas_2018"
GMT = DATA / "h.all.v2023.2.Hs.symbols.gmt"


# ------------------------------------------------------------------ TCGA genomic features
def genomic_sources():
    path = DATA / "tcga_ov_pancan_genomics.json"
    if not path.exists():
        cna = requests.post(f"{CBIO}/molecular-profiles/{STUDY}_gistic/discrete-copy-number/fetch"
                            "?discreteCopyNumberEventType=ALL&projection=SUMMARY",
                            json={"sampleListId": f"{STUDY}_cna", "entrezGeneIds": [898]}, timeout=180)
        cna.raise_for_status()
        profiled = requests.get(f"{CBIO}/sample-lists/{STUDY}_cna", timeout=120)
        profiled.raise_for_status()
        clin = requests.get(f"{CBIO}/studies/{STUDY}/clinical-data?clinicalDataType=SAMPLE&projection=SUMMARY"
                            "&pageSize=100000", timeout=180)
        clin.raise_for_status()
        keep = {"ANEUPLOIDY_SCORE", "FRACTION_GENOME_ALTERED", "TMB_NONSYNONYMOUS"}
        path.write_text(json.dumps({"ccne1": cna.json(), "cna_profiled": profiled.json()["sampleIds"],
                                    "clinical": [r for r in clin.json() if r["clinicalAttributeId"] in keep]}),
                        encoding="utf-8")
    return json.loads(path.read_text(encoding="utf-8"))


def genomic_table():
    src = genomic_sources()
    hrd = pd.read_csv(DATA / "TCGA.HRD_withSampleID.txt", sep="\t")
    hrd = hrd.loc[hrd.sampleID.str.endswith("-01")].assign(id=lambda d: d.sampleID.str[:12]).drop_duplicates("id")
    table = hrd.set_index("id")[["HRD"]].rename(columns={"HRD": "hrd"})
    profiled = {s[:12] for s in src["cna_profiled"] if s.endswith("-01")}
    amp = {r["patientId"] for r in src["ccne1"] if r["alteration"] == 2}
    ccne1 = pd.Series({p: float(p in amp) for p in profiled}, name="ccne1")
    clin = pd.DataFrame(src["clinical"])
    clin = clin.loc[clin.sampleId.str.endswith("-01")].pivot_table(index="patientId", columns="clinicalAttributeId",
                                                                   values="value", aggfunc="first").astype(float)
    clin.columns = ["aneuploidy", "fga", "tmb"]
    brca = json.loads((DATA / "tcga_ov_pub_brca.json").read_text(encoding="utf-8"))
    muts = pd.DataFrame(brca["mutations"])
    carriers = set(muts.loc[~muts.mutationType.eq("Silent"), "patientId"])
    b = pd.Series({s[:12]: float(s[:12] in carriers) for s in brca["sequenced"]}, name="brca")
    table = table.join(ccne1, how="outer").join(clin, how="outer").join(b, how="outer")
    table["hrd10"] = table.hrd / 10
    table["hrd_high"] = (table.hrd >= 42).where(table.hrd.notna()).astype(float)
    table["aneuploidy_sd"] = (table.aneuploidy - table.aneuploidy.mean()) / table.aneuploidy.std()
    table["fga_sd"] = (table.fga - table.fga.mean()) / table.fga.std()
    table["tmb_log2"] = np.log2(table.tmb + 1)
    return table


FEATURES = {"hrd10": "HRD score, per 10 points", "hrd_high": "HRD score ≥42", "ccne1": "CCNE1 amplification",
            "aneuploidy_sd": "Aneuploidy score, per SD", "fga_sd": "Fraction of genome altered, per SD",
            "tmb_log2": "Nonsynonymous mutation burden, per doubling"}


def genomics(raw):
    g = genomic_table()
    lm, _ = cohort_data(raw)
    d = lm["TCGA"].join(g, on="id")
    d["event24"], d["time24"] = (d.death.eq(1) & d.time.le(24)).astype(int), d.time.clip(upper=24)
    full = base_population(raw["TCGA"]).dropna(subset=["age10", "rd"]).join(g, on="id")
    prog = full.loc[full.progression_event.eq(1) & full.pfs.gt(0) & (full.os - full.pfs).gt(1e-6)].copy()
    prog["pps"], prog["ttp6"] = prog.os - prog.pfs, prog.pfs / 6
    rows = []

    def cox(label, data, cols, duration, event, terms, feature):
        data = data.dropna(subset=cols)
        try:
            m = CoxPHFitter().fit(data[cols + [duration, event]], duration, event)
        except Exception as exc:
            print("skip", label, feature, type(exc).__name__)
            return
        for term in terms:
            s = m.summary.loc[term]
            rows.append({"feature": feature, "analysis": label, "term": term, "n": len(data),
                         "events": int(data[event].sum()), "positive": int(data[term].sum()) if data[term].nunique() == 2 else "",
                         "estimate": s["exp(coef)"], "lo": s["exp(coef) lower 95%"], "hi": s["exp(coef) upper 95%"], "p": s["p"]})

    def increment(data, cols, feature):
        data = data.dropna(subset=cols).reset_index(drop=True)
        y, w, _, _ = horizon_data(data, 12)
        cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
        p0, p1 = np.empty(len(data)), np.empty(len(data))
        for learn, test in cv.split(data, y + 2 * (w == 0)):
            p0[test] = fit_pooled([data.iloc[learn]], STATUS).predict_proba(data.iloc[test][STATUS])[:, 1]
            p1[test] = fit_pooled([data.iloc[learn]], STATUS + cols).predict_proba(data.iloc[test][STATUS + cols])[:, 1]
        rng = np.random.default_rng(SEED + 31)
        t, e = data.time.to_numpy(), data.death.to_numpy()
        boot = []
        for _ in range(1000):
            i = rng.integers(0, len(t), len(t))
            yy, ww, _, controls = horizon_arrays(t[i], e[i], 12)
            if yy.sum() and controls:
                boot.append(roc_auc_score(yy, p1[i], sample_weight=ww) - roc_auc_score(yy, p0[i], sample_weight=ww))
        rows.append({"feature": feature, "analysis": "Out-of-fold change in AUC when added to Status", "term": "+".join(cols),
                     "n": len(data), "events": int(y.sum()), "positive": "",
                     "estimate": metrics(data, p1, 12)["auc"] - metrics(data, p0, 12)["auc"],
                     "lo": np.percentile(boot, 2.5), "hi": np.percentile(boot, 97.5), "p": np.nan})

    for key, label in FEATURES.items():
        cox("Progression-free survival from diagnosis", full, BASE + [key], "pfs", "progression_event", [key], label)
        cox("Overall survival from diagnosis", full, BASE + [key], "os", "death", [key], label)
        cox("Death within 24 months after the month-12 landmark, status adjusted", d, STATUS + [key], "time24", "event24", [key], label)
        cox("Death after the month-12 landmark, all follow-up, status adjusted", d, STATUS + [key], "time", "death", [key], label)
        cox("Post-progression survival, progression time adjusted", prog, BASE + ["ttp6", key], "pps", "death", [key], label)
        increment(d, [key], label)
    joint = ["brca", "hrd10", "ccne1"]
    name = "Joint model: BRCA1/2, HRD score and CCNE1"
    cox("Death after the month-12 landmark, all follow-up, status adjusted", d, STATUS + joint, "time", "death",
        ["early", *joint], name)
    cox("Post-progression survival, progression time adjusted", prog, BASE + ["ttp6", *joint], "pps", "death",
        ["ttp6", *joint], name)
    cox("Progression-free survival from diagnosis", full, BASE + joint, "pfs", "progression_event", joint, name)
    increment(d, joint, name)
    pd.DataFrame(rows).to_csv(OUT / "genomics.csv", index=False)
    lmg = d.dropna(subset=["hrd"])
    desc = {"tcga_landmark_with_hrd": len(lmg), "hrd_median": lmg.hrd.median(), "hrd_q1": lmg.hrd.quantile(.25),
            "hrd_q3": lmg.hrd.quantile(.75), "hrd_high": int(lmg.hrd_high.sum()),
            "tcga_landmark_with_cna": int(d.ccne1.notna().sum()), "ccne1_amplified": int(d.ccne1.sum()),
            "hrd_in_brca_carriers": d.loc[d.brca.eq(1), "hrd"].median(), "hrd_in_noncarriers": d.loc[d.brca.eq(0), "hrd"].median(),
            "hrd_in_ccne1_amp": d.loc[d.ccne1.eq(1), "hrd"].median(), "hrd_in_ccne1_not_amp": d.loc[d.ccne1.eq(0), "hrd"].median(),
            "early_progression_ccne1_amp": d.loc[d.ccne1.eq(1), "early"].mean(),
            "early_progression_ccne1_not_amp": d.loc[d.ccne1.eq(0), "early"].mean(),
            "early_progression_hrd_high": d.loc[d.hrd_high.eq(1), "early"].mean(),
            "early_progression_hrd_low": d.loc[d.hrd_high.eq(0), "early"].mean()}
    pd.DataFrame([desc]).to_csv(OUT / "genomics_description.csv", index=False)
    print("genomics done", flush=True)


# ------------------------------------------------------------------ stroma-adjusted transcriptome comparison
def adjusted_transcriptome(raw, expr):
    frames = with_scores(raw, expr)
    zs, weights = {}, {}
    for cohort, x in expr.items():
        d = frames[cohort]
        d = d.loc[(d.progression_event.eq(1) & d.pfs.lt(12)) | d.pfs.ge(12)].copy()
        d["fast"] = (d.progression_event.eq(1) & d.pfs.lt(12)).astype(float)
        y = x[d.id]
        y = y.loc[y.notna().mean(axis=1).gt(.9)]
        y = y.apply(lambda r: r.fillna(r.mean()), axis=1) if y.isna().any().any() else y
        design = np.column_stack((np.ones(len(d)), d.fast, d.stroma))
        inv = np.linalg.inv(design.T @ design)
        beta = y.to_numpy() @ design @ inv
        resid = y.to_numpy() - beta @ design.T
        dof = len(d) - design.shape[1]
        se = np.sqrt((resid ** 2).sum(axis=1) / dof * inv[1, 1])
        t = beta[:, 1] / se
        p = 2 * stats.t.sf(np.abs(t), dof)
        zs[cohort] = pd.Series(np.sign(t) * stats.norm.isf(np.clip(p, 1e-300, 1) / 2), index=y.index).dropna()
        n1, n0 = d.fast.sum(), (1 - d.fast).sum()
        weights[cohort] = (4 / (1 / n1 + 1 / n0)) ** .5
    sets = {"primary": list(PRIMARY), "additional": [c for c in expr if c not in PRIMARY]}
    res = {}
    for name, cohorts in sets.items():
        frame = pd.concat({c: zs[c] for c in cohorts}, axis=1).dropna()
        w = np.array([weights[c] for c in cohorts])
        rank = ((frame * w).sum(axis=1) / np.sqrt((w ** 2).sum())).sort_values(ascending=False)
        pre = gseapy.prerank(rnk=rank, gene_sets=str(GMT), min_size=15, max_size=500, permutation_num=1000,
                             seed=SEED % 2**31, outdir=None, threads=4, verbose=False)
        res[name] = pre.res2d.set_index("Term")[["NES", "NOM p-val", "FDR q-val"]].astype(float)
    adj = res["primary"].join(res["additional"], lsuffix="_primary", rsuffix="_additional", how="inner")
    raw_table = pd.read_csv(OUT / "transcriptome_hallmarks.csv", index_col=0)
    out = raw_table[["NES_primary", "FDR q-val_primary", "NES_additional", "NOM p-val_additional", "replicated"]].join(
        adj.add_prefix("adj_"), how="inner")
    out["adj_replicated"] = (out["adj_FDR q-val_primary"].lt(.05) & out["adj_NOM p-val_additional"].lt(.05)
                             & (np.sign(out.adj_NES_primary) == np.sign(out.adj_NES_additional)))
    out.sort_values("NES_primary", ascending=False).to_csv(OUT / "transcriptome_hallmarks_stroma_adjusted.csv")
    print("adjusted transcriptome done", flush=True)


# ------------------------------------------------------------------ single-cell attribution
def single_cell():
    counts = pd.read_csv(OUT / "sc_pseudobulk_counts.csv.gz", index_col=0)
    cells = pd.read_csv(OUT / "sc_cell_counts.csv")
    cells = cells.loc[cells.treatment_phase.eq("treatment-naive") & cells.cells.ge(20)]
    keep = [f"{r.patient_id}|treatment-naive|{r.cell_class}" for r in cells.itertuples()]
    cpm = np.log2(counts[keep] / counts[keep].sum() * 1e6 + 1)
    classes = cells.cell_class.to_numpy()
    usable = [c for c in pd.unique(classes) if (classes == c).sum() >= 3]
    mean = pd.DataFrame({c: cpm.loc[:, classes == c].mean(axis=1) for c in usable})
    mean = mean.loc[mean.max(axis=1) > np.log2(11)]
    spec = mean.sub(mean.mean(axis=1), axis=0).div(mean.std(axis=1), axis=0)
    mean.to_csv(OUT / "sc_celltype_logcpm.csv")
    pd.DataFrame({"cell_class": usable, "patients": [(classes == c).sum() for c in usable],
                  "cells": [int(cells.loc[cells.cell_class.eq(c), "cells"].sum()) for c in usable]}
                 ).to_csv(OUT / "sc_classes.csv", index=False)

    bulk = pd.read_csv(OUT / "transcriptome_gene_z.csv", index_col=0)["all_cohorts"].dropna()
    universe = bulk.index.intersection(spec.index)
    bulk = bulk[universe]
    hall = {}
    for line in GMT.read_text().splitlines():
        parts = line.split("\t")
        hall[parts[0]] = parts[2:]
    from v2_molecular import SIGNATURES
    order = bulk.sort_values(ascending=False).index
    gene_sets = {
        "Higher in early progression, top 200 genes": list(order[:200]),
        "Lower in early progression, top 200 genes": list(order[-200:]),
        "Hallmark epithelial–mesenchymal transition": hall["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"],
        "Hallmark angiogenesis": hall["HALLMARK_ANGIOGENESIS"],
        "Hallmark UV response down": hall["HALLMARK_UV_RESPONSE_DN"],
        "Stromal score genes": SIGNATURES["stroma"], "T-cell score genes": SIGNATURES["tcell"],
        "Proliferation score genes": SIGNATURES["prolif"],
    }
    rng = np.random.default_rng(SEED + 808)
    top = spec.loc[universe].idxmax(axis=1)
    base_share = top.value_counts(normalize=True)
    rows = []
    for name, genes in gene_sets.items():
        pool = universe if not name.endswith("score genes") else spec.index
        genes = [g for g in genes if g in pool]
        obs = spec.loc[genes].mean()
        null = np.array([spec.loc[rng.choice(pool, len(genes), replace=False)].mean().to_numpy() for _ in range(5000)])
        share = spec.loc[genes].idxmax(axis=1).value_counts(normalize=True)
        for i, c in enumerate(usable):
            rows.append({"gene_set": name, "genes": len(genes), "cell_class": c, "mean_specificity": obs[c],
                         "p_enriched": (1 + (null[:, i] >= obs[c]).sum()) / 5001,
                         "share_highest": share.get(c, 0.0), "share_highest_all_genes": base_share.get(c, 0.0)})
    pd.DataFrame(rows).to_csv(OUT / "sc_attribution.csv", index=False)
    corr = []
    for c in usable:
        r = stats.spearmanr(bulk, spec.loc[universe, c])
        corr.append({"cell_class": c, "genes": len(universe), "spearman": r.statistic, "p": r.pvalue})
    pd.DataFrame(corr).to_csv(OUT / "sc_bulk_correlation.csv", index=False)
    print("single cell done", flush=True)


if __name__ == "__main__":
    import sys
    sources = all_sources()
    if "genomics" in sys.argv or len(sys.argv) == 1:
        genomics(sources)
    if "adjusted" in sys.argv or len(sys.argv) == 1:
        adjusted_transcriptome(sources, expression())
    if "sc" in sys.argv or len(sys.argv) == 1:
        single_cell()

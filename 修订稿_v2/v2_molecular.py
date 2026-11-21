"""Second-revision molecular analyses on the tumour expression profiles of the same patients.

D1  Prespecified stromal, T-cell and proliferation scores: association and incremental value.
D2  BRCA1/2 status in the sequenced TCGA subset.
D3  Cross-cohort transcriptome of progression within 12 months with hallmark enrichment.
"""
import gzip
import io
import json
import warnings
import zipfile

import gseapy
from lifelines import CoxPHFitter
import numpy as np
import pandas as pd
import rdata
from scipy import stats
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

from v2_common import (BASE, DATA, EXTRA, OUT, PRIMARY, ROOT, SEED, STATUS, all_sources, base_population,
                       cohort_data, fit_pooled, horizon_arrays, horizon_data, meta, metrics)

warnings.filterwarnings("ignore")
RAW = ROOT / "raw"
SIGNATURES = {
    "stroma": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL11A1", "FAP", "POSTN", "FN1", "THBS2", "VCAN",
               "SPARC", "DCN", "LUM", "PDGFRB", "ACTA2", "TAGLN", "INHBA", "MMP2", "CTSK", "FBN1", "TIMP3"],
    "tcell": ["CD8A", "CD3D", "CD3E", "CD2", "CD247", "GZMA", "GZMB", "GZMK", "PRF1", "NKG7", "CXCL9", "CXCL10",
              "CXCL11", "CXCL13"],
    "prolif": ["MKI67", "TOP2A", "BUB1", "CCNB1", "CDK1", "AURKA", "BIRC5", "MCM2", "TYMS", "PCNA", "FOXM1", "CENPF"],
}
LABEL = {"stroma": "Stromal score", "tcell": "T-cell score", "prolif": "Proliferation score"}


def eset(path, key=None):
    e = next(iter(rdata.read_rda(path).values()))
    x = e.assayData["exprs"]
    frame = pd.DataFrame(np.asarray(x), index=[str(v) for v in x.coords[x.dims[0]].values],
                         columns=[str(v) for v in x.coords[x.dims[1]].values])
    if key:
        pheno = e.phenoData.data
        pheno.columns = pheno.columns.map(str)
        frame.columns = pheno[key].astype(str).to_numpy()
    return frame


def expression():
    out = {}
    out["TCGA"] = eset(RAW / "TCGA_eset.rda", "unique_patient_ID")
    out["AOCS"] = eset(DATA / "GSE9891_eset.rda")
    out["Japan"] = eset(DATA / "GSE32062.GPL6480_eset.rda")
    for label, accession in EXTRA.items():
        out[label] = eset(ROOT / "修订稿" / "data" / f"{accession}_eset.rda")
    # ICON7: Illumina probes to symbols, probe mean per gene.
    lines = (DATA / "GPL14951_full.txt").read_text(encoding="utf-8", errors="replace").splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_begin"))
    end = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_end"))
    gpl = pd.read_csv(io.StringIO("\n".join(lines[start + 1:end])), sep="\t", usecols=["ID", "Symbol"], dtype=str)
    with gzip.open(RAW / "GSE140082_series_matrix.txt.gz", "rt", errors="replace") as handle:
        text = handle.read()
    body = text.split("!series_matrix_table_begin\n", 1)[1].split("!series_matrix_table_end", 1)[0]
    icon = pd.read_csv(io.StringIO(body), sep="\t", index_col=0)
    symbol = gpl.dropna().set_index("ID").Symbol
    icon = icon.loc[icon.index.intersection(symbol.index)]
    out["ICON7"] = icon.groupby(symbol.reindex(icon.index).to_numpy()).mean()
    # UKE: raw counts to log2 counts per million, one primary sample per patient as in the clinical loader.
    with zipfile.ZipFile(RAW / "UKE_supplement.zip") as outer:
        nested = next(n for n in outer.namelist() if n.endswith("/Table_S1.zip"))
        with zipfile.ZipFile(io.BytesIO(outer.read(nested))) as inner:
            counts = pd.read_csv(inner.open("Tables_S1/Table_S1a.UKE_rawcounts.csv"), index_col=0)
            md = pd.read_csv(inner.open("Tables_S1/Table_S1b.UKE_metadata.tsv"), sep="\t")
    md = md.loc[md.patient_id.notna() & md.recurrent_tumor.eq("no")]
    md = md.sort_values(["patient_id", "tumor_sample"]).groupby("patient_id", as_index=False).first()
    md = md.loc[md.tumor_sample.isin(counts.columns)]
    cpm = np.log2(counts / counts.sum() * 1e6 + 1)
    hgnc = pd.read_csv(DATA / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna()
    ens = hgnc.drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id").symbol
    cpm = cpm.loc[cpm.index.intersection(ens.index)]
    cpm = cpm.loc[(counts.loc[cpm.index] >= 10).mean(axis=1) >= .2]
    cpm = cpm.groupby(ens.reindex(cpm.index).to_numpy()).mean()[md.tumor_sample]
    cpm.columns = md.patient_id.astype(int).astype(str).to_numpy()
    out["UKE"] = cpm
    return {k: v.loc[:, ~v.columns.duplicated()].dropna(how="all") for k, v in out.items()}


def with_scores(raw, expr):
    frames, coverage = {}, []
    for cohort, x in expr.items():
        d = base_population(raw[cohort])
        d = d.loc[d.id.isin(x.columns)].reset_index(drop=True)
        sub = x[d.id]
        z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
        row = {"cohort": cohort, "eligible_with_expression": len(d), "genes": len(x)}
        for name, genes in SIGNATURES.items():
            present = [g for g in genes if g in z.index]
            score = z.loc[present].mean()
            d[name] = ((score - score.mean()) / score.std()).to_numpy()
            row[f"{name}_genes"] = f"{len(present)}/{len(genes)}"
        coverage.append(row)
        frames[cohort] = d.assign(cohort=cohort)
    pd.DataFrame(coverage).to_csv(OUT / "molecular_coverage.csv", index=False)
    return frames


def adjusters(d, extra=()):
    return [c for c in ("age10", "iv", "rd", *extra) if d[c].notna().mean() > .8 and d[c].nunique() > 1]


def cox_term(d, duration, event, term, cols):
    d = d.dropna(subset=cols + [term])
    m = CoxPHFitter().fit(d[cols + [term, duration, event]], duration, event)
    s = m.summary.loc[term]
    return {"n": len(d), "events": int(d[event].sum()), "coef": s["coef"], "se": s["se(coef)"]}


def landmark_frame(d, landmark=12):
    known = d.os.gt(landmark) & (d.pfs.ge(landmark) | (d.progression_event.eq(1) & d.pfs.lt(landmark)))
    f = d.loc[known].copy()
    f["early"] = (f.pfs.le(landmark) & f.progression_event.eq(1)).astype(int)
    f["time"] = f.os - landmark
    f["event24"] = (f.death.eq(1) & f.time.le(24)).astype(int)
    f["time24"] = f.time.clip(upper=24)
    return f


def score_associations(frames):
    rows = []
    for cohort, d in frames.items():
        d = d.copy()
        lm = landmark_frame(d)
        prog = d.loc[d.progression_event.eq(1) & d.pfs.gt(0) & (d.os - d.pfs).gt(1e-6)].copy()
        prog["pps"], prog["ttp6"] = prog.os - prog.pfs, prog.pfs / 6
        for name in SIGNATURES:
            analyses = {
                "Progression-free survival from origin": (d, "pfs", "progression_event", adjusters(d)),
                "Overall survival from origin": (d, "os", "death", adjusters(d)),
                "Death within 24 months after the month-12 landmark, baseline adjusted": (lm, "time24", "event24", adjusters(lm)),
                "Death within 24 months after the month-12 landmark, baseline and status adjusted":
                    (lm, "time24", "event24", adjusters(lm, ("early",))),
                "Post-progression survival, baseline and progression time adjusted": (prog, "pps", "death", adjusters(prog, ("ttp6",))),
            }
            for analysis, (data, t, e, cols) in analyses.items():
                try:
                    rows.append({"cohort": cohort, "score": name, "analysis": analysis, "adjusted_for": "+".join(cols),
                                 **cox_term(data, t, e, name, cols)})
                except Exception as exc:
                    print("skip", cohort, name, analysis, type(exc).__name__)
    table = pd.DataFrame(rows)
    table["hr"], table["lo"], table["hi"] = (np.exp(table.coef), np.exp(table.coef - 1.96 * table.se),
                                             np.exp(table.coef + 1.96 * table.se))
    table.to_csv(OUT / "molecular_cohort_hr.csv", index=False)
    pooled = []
    for (score, analysis), g in table.groupby(["score", "analysis"], sort=False):
        r = meta(g.coef, g.se)
        pooled.append({"score": score, "analysis": analysis, "k": r["k"], "n": int(g.n.sum()), "events": int(g.events.sum()),
                       "hr": np.exp(r["estimate"]), "lo": np.exp(r["lo"]), "hi": np.exp(r["hi"]), "i2": r["i2"],
                       "fixed_hr": np.exp(r["fixed"]), "fixed_lo": np.exp(r["fixed_lo"]), "fixed_hi": np.exp(r["fixed_hi"])})
    pd.DataFrame(pooled).to_csv(OUT / "molecular_pooled_hr.csv", index=False)


def score_increment(raw, frames):
    """Internal-external ΔAUC from adding expression scores to the Status model in the four primary cohorts."""
    lm, _ = cohort_data(raw)
    for cohort in PRIMARY:
        scores = frames[cohort].set_index("id")[list(SIGNATURES)]
        lm[cohort] = lm[cohort].join(scores, on="id").dropna(subset=list(SIGNATURES)).reset_index(drop=True)
    rng = np.random.default_rng(SEED + 404)
    rows = []
    additions = {**{k: [k] for k in SIGNATURES}, "all three scores": list(SIGNATURES)}
    for held in PRIMARY:
        train = [d for c, d in lm.items() if c != held]
        test = lm[held]
        base = fit_pooled(train, STATUS).predict_proba(test[STATUS])[:, 1]
        base0 = fit_pooled(train, BASE).predict_proba(test[BASE])[:, 1]
        t, e = test.time.to_numpy(), test.death.to_numpy()
        for name, cols in additions.items():
            for ref_name, ref, ref_cols in (("Status", base, STATUS), ("Baseline", base0, BASE)):
                new = fit_pooled(train, ref_cols + cols).predict_proba(test[ref_cols + cols])[:, 1]
                point = metrics(test, new, 12)["auc"] - metrics(test, ref, 12)["auc"]
                boot = []
                for _ in range(500):
                    i = rng.integers(0, len(t), len(t))
                    y, w, _, controls = horizon_arrays(t[i], e[i], 12)
                    if y.sum() and controls:
                        boot.append(roc_auc_score(y, new[i], sample_weight=w) - roc_auc_score(y, ref[i], sample_weight=w))
                rows.append({"cohort": held, "added": name, "reference": ref_name, "n": len(test),
                             "deaths": metrics(test, ref, 12)["deaths"], "auc_reference": metrics(test, ref, 12)["auc"],
                             "delta": point, "se": np.std(boot, ddof=1), "lo": np.percentile(boot, 2.5),
                             "hi": np.percentile(boot, 97.5)})
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "molecular_increment.csv", index=False)
    pooled = []
    for (added, ref), g in table.groupby(["added", "reference"], sort=False):
        r = meta(g.delta, g.se)
        pooled.append({"added": added, "reference": ref, "n": int(g.n.sum()), "deaths": int(g.deaths.sum()),
                       **{k: r[k] for k in ("estimate", "lo", "hi", "i2", "fixed", "fixed_lo", "fixed_hi")}})
    pd.DataFrame(pooled).to_csv(OUT / "molecular_increment_pooled.csv", index=False)


def brca(raw):
    src = json.loads((DATA / "tcga_ov_pub_brca.json").read_text(encoding="utf-8"))
    sequenced = {s[:12] for s in src["sequenced"]}
    muts = pd.DataFrame(src["mutations"])
    muts = muts.loc[~muts.mutationType.isin(["Silent"])]
    carriers = set(muts.patientId)
    germline = set(muts.loc[muts.mutationStatus.str.lower().eq("germline"), "patientId"])
    lm, _ = cohort_data(raw)
    d = lm["TCGA"]
    d = d.loc[d.id.isin(sequenced)].reset_index(drop=True)
    d["brca"] = d.id.isin(carriers).astype(int)
    d["event24"], d["time24"] = (d.death.eq(1) & d.time.le(24)).astype(int), d.time.clip(upper=24)
    d["event"] = d.death
    d["brca_early"] = d.brca * d.early
    rows = []

    def add(label, data, cols, duration, event, terms):
        m = CoxPHFitter().fit(data[cols + [duration, event]], duration, event)
        for term in terms:
            s = m.summary.loc[term]
            rows.append({"analysis": label, "term": term, "n": len(data), "events": int(data[event].sum()),
                         "carriers": int(data.brca.sum()), "hr": s["exp(coef)"], "lo": s["exp(coef) lower 95%"],
                         "hi": s["exp(coef) upper 95%"], "p": s["p"]})

    add("Landmark, death within 24 months", d, STATUS + ["brca"], "time24", "event24", ["early", "brca"])
    add("Landmark, all follow-up", d, STATUS + ["brca"], "time", "event", ["early", "brca"])
    add("Landmark, all follow-up, with interaction", d, STATUS + ["brca", "brca_early"], "time", "event",
        ["early", "brca", "brca_early"])
    full = base_population(raw["TCGA"]).dropna(subset=["age10", "rd"])
    full = full.loc[full.id.isin(sequenced)].copy()
    full["brca"] = full.id.isin(carriers).astype(int)
    prog = full.loc[full.progression_event.eq(1) & full.pfs.gt(0) & (full.os - full.pfs).gt(1e-6)].copy()
    prog["pps"], prog["ttp6"] = prog.os - prog.pfs, prog.pfs / 6
    add("Post-progression survival", prog, BASE + ["ttp6", "brca"], "pps", "death", ["ttp6", "brca"])
    add("Progression-free survival from diagnosis", full, BASE + ["brca"], "pfs", "progression_event", ["brca"])
    # Out-of-fold increment at the 12-month horizon.
    y, w, _, _ = horizon_data(d, 12)
    cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
    p0, p1 = np.empty(len(d)), np.empty(len(d))
    for learn, test in cv.split(d, y + 2 * (w == 0)):
        p0[test] = fit_pooled([d.iloc[learn]], STATUS).predict_proba(d.iloc[test][STATUS])[:, 1]
        p1[test] = fit_pooled([d.iloc[learn]], STATUS + ["brca"]).predict_proba(d.iloc[test][STATUS + ["brca"]])[:, 1]
    rng = np.random.default_rng(SEED + 9)
    t, e = d.time.to_numpy(), d.death.to_numpy()
    boot = []
    for _ in range(1000):
        i = rng.integers(0, len(t), len(t))
        yy, ww, _, controls = horizon_arrays(t[i], e[i], 12)
        if yy.sum() and controls:
            boot.append(roc_auc_score(yy, p1[i], sample_weight=ww) - roc_auc_score(yy, p0[i], sample_weight=ww))
    rows.append({"analysis": "Out-of-fold ΔAUC from adding BRCA1/2 to Status", "term": "brca", "n": len(d),
                 "events": int(y.sum()), "carriers": int(d.brca.sum()),
                 "hr": metrics(d, p1, 12)["auc"] - metrics(d, p0, 12)["auc"],
                 "lo": np.percentile(boot, 2.5), "hi": np.percentile(boot, 97.5), "p": np.nan,
                 "note": f"AUC {metrics(d, p0, 12)['auc']:.3f} to {metrics(d, p1, 12)['auc']:.3f}; germline carriers "
                         f"{int(d.id.isin(germline).sum())}"})
    pd.DataFrame(rows).to_csv(OUT / "brca.csv", index=False)


def transcriptome(raw, expr):
    """Gene-wise comparison of progression within 12 months against progression-free at 12 months."""
    zs, sizes = {}, []
    for cohort, x in expr.items():
        d = base_population(raw[cohort])
        d = d.loc[d.id.isin(x.columns)]
        early = d.loc[d.progression_event.eq(1) & d.pfs.lt(12)]
        free = d.loc[d.pfs.ge(12)]
        a, b = x[early.id], x[free.id]
        keep = a.notna().mean(axis=1).gt(.9) & b.notna().mean(axis=1).gt(.9)
        t, p = stats.ttest_ind(a.loc[keep], b.loc[keep], axis=1, equal_var=False, nan_policy="omit")
        z = pd.Series(np.sign(t) * stats.norm.isf(np.clip(p, 1e-300, 1) / 2), index=a.index[keep]).dropna()
        zs[cohort] = z
        sizes.append({"cohort": cohort, "progressed_before_12": len(early), "progression_free_at_12": len(free),
                      "genes": len(z), "weight": (4 / (1 / len(early) + 1 / len(free))) ** .5})
    sizes = pd.DataFrame(sizes).set_index("cohort")
    sizes.to_csv(OUT / "transcriptome_sizes.csv")

    def stouffer(cohorts):
        frame = pd.concat({c: zs[c] for c in cohorts}, axis=1).dropna()
        w = sizes.loc[cohorts, "weight"].to_numpy()
        return (frame * w).sum(axis=1) / np.sqrt((w ** 2).sum())

    sets = {"Primary cohorts": list(PRIMARY), "Additional cohorts": [c for c in expr if c not in PRIMARY]}
    ranks = {name: stouffer(cohorts) for name, cohorts in sets.items()}
    genes = pd.concat(ranks, axis=1)
    genes["all_cohorts"] = stouffer(list(expr))
    both = genes.dropna()
    genes.to_csv(OUT / "transcriptome_gene_z.csv")
    res = {}
    for name, rank in ranks.items():
        pre = gseapy.prerank(rnk=rank.sort_values(ascending=False), gene_sets=str(DATA / "h.all.v2023.2.Hs.symbols.gmt"),
                             min_size=15, max_size=500, permutation_num=1000, seed=SEED % 2**31, outdir=None,
                             threads=4, verbose=False)
        res[name] = pre.res2d.set_index("Term")[["NES", "NOM p-val", "FDR q-val"]].astype(float)
    table = res["Primary cohorts"].join(res["Additional cohorts"], lsuffix="_primary", rsuffix="_additional", how="inner")
    table["replicated"] = (table["FDR q-val_primary"].lt(.05) & table["NOM p-val_additional"].lt(.05)
                           & (np.sign(table.NES_primary) == np.sign(table.NES_additional)))
    table.sort_values("NES_primary", ascending=False).to_csv(OUT / "transcriptome_hallmarks.csv")
    rho = stats.spearmanr(both["Primary cohorts"], both["Additional cohorts"])
    nes = stats.spearmanr(table.NES_primary, table.NES_additional)
    top = both["all_cohorts"].abs().sort_values(ascending=False).head(40).index
    pd.DataFrame([{"genes_in_both": len(both), "gene_z_spearman": rho.statistic, "gene_z_p": rho.pvalue,
                   "hallmark_nes_spearman": nes.statistic, "hallmark_nes_p": nes.pvalue,
                   "hallmarks_tested": len(table), "hallmarks_fdr05_primary": int(table["FDR q-val_primary"].lt(.05).sum()),
                   "hallmarks_replicated": int(table.replicated.sum()),
                   "max_abs_z_all": both["all_cohorts"].abs().max(),
                   "genes_bonferroni_all": int((both["all_cohorts"].abs() > stats.norm.isf(.025 / len(both))).sum()),
                   "top_genes": ", ".join(f"{g} ({both.loc[g, 'all_cohorts']:+.1f})" for g in top[:15])}]
                 ).to_csv(OUT / "transcriptome_summary.csv", index=False)


if __name__ == "__main__":
    sources = all_sources()
    matrices = expression()
    for name, matrix in matrices.items():
        print(name, matrix.shape, "matched", int(base_population(sources[name]).id.isin(matrix.columns).sum()),
              "of", len(base_population(sources[name])), flush=True)
    scored = with_scores(sources, matrices)
    score_associations(scored)
    print("D1 associations done", flush=True)
    score_increment(sources, scored)
    print("D1 increment done", flush=True)
    brca(sources)
    print("D2 done", flush=True)
    transcriptome(sources, matrices)
    print("molecular analyses complete")

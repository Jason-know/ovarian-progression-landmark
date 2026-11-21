"""Figures for the second revision."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RES, FIG = HERE / "results", HERE / "figures"
BLUE, TEAL, RED, INK, GRAY, GOLD = "#315B9B", "#197A6F", "#AE4346", "#1B2639", "#9AA8B4", "#B9812A"
PALETTE = [RED, GOLD, TEAL, BLUE, "#5B4B8A"]
COHORTS = ["TCGA", "AOCS", "UKE", "ICON7"]
plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "pdf.fonttype": 42, "text.color": INK})


def read(name):
    return pd.read_csv(RES / f"{name}.csv")


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=300, bbox_inches="tight", pad_inches=.12)
    plt.close(fig)


def title(ax, text):
    ax.set_title(text, loc="left", fontweight="bold", pad=7)


def forest(ax, labels, values, ref, colors=None, xlabel="", log=False, xlim=None, markers=None):
    ax.axvline(ref, color=GRAY, linestyle="--", linewidth=1)
    for i, (v, lo, hi) in enumerate(values):
        ax.errorbar(v, i, xerr=[[max(0, v - lo)], [max(0, hi - v)]], fmt=(markers or ["o"] * len(values))[i],
                    capsize=2.6, lw=1.4, ms=5.5, color=(colors or [BLUE] * len(values))[i])
    if log:
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set(yticks=np.arange(len(labels)), yticklabels=labels, ylim=(len(labels) - .45, -.6), xlabel=xlabel)
    if xlim:
        ax.set_xlim(*xlim)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#E1E6EA", linewidth=.6)


def performance():
    c = read("iecv_cohort")
    c = c.loc[c.design.eq("IECV")].set_index("cohort")
    p = read("iecv_pooled")
    p = p.loc[p.design.eq("IECV")].set_index("measure")
    fig, ax = plt.subplots(2, 2, figsize=(10.6, 6.6), layout="constrained")
    labels = COHORTS + ["Random-effects summary"]
    colors, markers = [BLUE] * 4 + [INK], ["o"] * 4 + ["D"]

    def block(key, measure):
        rows = [c.loc[k, [key, f"{key}_lo", f"{key}_hi"]] for k in COHORTS]
        rows.append(p.loc[measure, ["estimate", "lo", "hi"]])
        return rows

    forest(ax[0, 0], labels, block("d_status_baseline", "ΔAUC Status − Baseline"), 0, colors,
           "Change in time-dependent AUC", markers=markers)
    title(ax[0, 0], "A   Progression status added to the clinical baseline")
    forest(ax[0, 1], labels, block("d_status_progonly", "ΔAUC Status − Progression only"), 0, colors,
           "Change in time-dependent AUC", markers=markers)
    title(ax[0, 1], "B   Clinical baseline added to progression status")
    locked = read("iecv_cohort")
    locked = locked.loc[locked.design.eq("TCGA only (locked)")].set_index("cohort")
    rec = read("recalibration").query("method == 'Intercept update'").set_index("cohort")
    a = ax[1, 0]
    a.axvline(1, color=GRAY, linestyle="--", linewidth=1)
    for i, k in enumerate(COHORTS):
        r = c.loc[k]
        a.errorbar(r.oe, i - .2, xerr=[[r.oe - r.oe_lo], [r.oe_hi - r.oe]], fmt="o", color=BLUE, capsize=2.6, ms=5.5,
                   label="Internal–external (held-out cohort)" if i == 0 else None)
        if k in locked.index:
            r = locked.loc[k]
            a.errorbar(r.oe, i, xerr=[[r.oe - r.oe_lo], [r.oe_hi - r.oe]], fmt="s", color=RED, capsize=2.6, ms=5,
                       label="TCGA-only equation" if i == 1 else None)
            a.plot(rec.loc[k, "oe"], i + .2, "^", color=TEAL, ms=6,
                   label="TCGA-only, cross-validated intercept update" if i == 1 else None)
    a.set(yticks=range(4), yticklabels=COHORTS, ylim=(3.5, -.6), xlabel="Observed-to-expected 12-month mortality",
          xlim=(.3, 2.95))
    a.tick_params(axis="y", length=0)
    a.grid(axis="x", color="#E1E6EA", linewidth=.6)
    a.legend(frameon=False, fontsize=7.4, loc="upper right")
    title(a, "C   Calibration in the large")
    curves = read("calibration_curves").query("model == 'IECV'")
    a = ax[1, 1]
    a.plot([0, .8], [0, .8], color=GRAY, linestyle="--", linewidth=1)
    for k, color in zip(COHORTS, PALETTE):
        d = curves.loc[curves.cohort.eq(k)]
        a.plot(d.predicted, d.observed, color=color, lw=1.8, label=k)
    a.set(xlabel="Predicted 12-month mortality", ylabel="Observed (smoothed, censoring-weighted)", xlim=(0, .8),
          ylim=(0, 1))
    a.legend(frameon=False, fontsize=8)
    a.grid(color="#E1E6EA", linewidth=.6)
    title(a, "D   Calibration curves, held-out cohorts")
    save(fig, "Figure_2_performance")


def landmark():
    dur = read("landmark_durations")
    hr = read("supermodel_hr")
    pooled = read("landmark_iecv_pooled")
    dist = read("landmark_timing_distribution").groupby("landmark")[["progressed", "progressed_ge6"]].sum()
    lms = sorted(dur.lm.unique())
    fig, ax = plt.subplots(2, 2, figsize=(10.6, 6.8), layout="constrained")
    a = ax[0, 0]
    parts = a.boxplot([dur.loc[dur.lm.eq(l), "duration"] for l in lms], positions=lms, widths=1.7, patch_artist=True,
                      showfliers=False, medianprops={"color": INK})
    for b, l in zip(parts["boxes"], lms):
        b.set(facecolor=RED if l == 12 else "#C9D6E6", edgecolor=INK, linewidth=.8)
    for l in lms:
        a.text(l, a.get_ylim()[1] if False else 25.3, f"{dist.loc[l, 'progressed']}\n({dist.loc[l, 'progressed_ge6']})",
               ha="center", va="top", fontsize=7.2)
    a.axhline(6, color=GRAY, linestyle=":", linewidth=1)
    a.set(xlabel="Landmark, months from origin", ylabel="Months since first progression", xticks=lms, ylim=(-.5, 26))
    title(a, "A   Spread of the timing variable among progressed patients")
    single = hr.loc[hr.model.eq("Single landmark")].assign(landmark=lambda d: d.landmark.astype(int))
    a = ax[0, 1]
    t = single.loc[single.term.eq("timing6") & single.landmark.ge(9)]
    a.axhline(1, color=GRAY, linestyle="--", linewidth=1)
    a.errorbar(t.landmark, t.hr, yerr=[t.hr - t.lo, t.hi - t.hr], fmt="o", color=BLUE, capsize=2.6, ms=5.5)
    sm = hr.loc[hr.model.str.startswith("Supermodel, landmarks 6") & hr.term.eq("timing6")].iloc[0]
    a.axhspan(sm.lo, sm.hi, color=TEAL, alpha=.14)
    a.axhline(sm.hr, color=TEAL, lw=1.3)
    a.text(27.3, sm.hi * 1.04, "supermodel, 95% CI", color=TEAL, fontsize=7.5, va="bottom", ha="right")
    a.set(yscale="log", xlabel="Landmark, months from origin", xticks=lms[1:], xlim=(7.5, 27.5),
          ylabel="Hazard ratio per 6 months\nlonger since progression")
    a.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
    a.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    a.set_yticks([.5, 1, 2, 4, 6])
    title(a, "B   Timing association at each landmark")
    for a, name, text in ((ax[1, 0], "status", "C   Status added to baseline, by landmark"),
                          (ax[1, 1], "timing", "D   Timing added to status, by landmark")):
        d = pooled.loc[pooled.contrast.eq(name)]
        a.axhline(0, color=GRAY, linestyle="--", linewidth=1)
        a.errorbar(d.landmark, d.estimate, yerr=[d.estimate - d.lo, d.hi - d.estimate], fmt="D", color=INK,
                   capsize=2.6, ms=5)
        a.set(xlabel="Landmark, months from origin", xticks=lms,
              ylabel="Pooled held-out change in AUC")
        a.grid(axis="y", color="#E1E6EA", linewidth=.6)
        title(a, text)
    save(fig, "Figure_3_landmarks")


def post_progression():
    km, bands = read("pps_km"), read("pps_bands")
    spline, cohort, pooled, conc = read("pps_spline"), read("pps_cohort_hr"), read("pps_pooled"), read("pps_concordance")
    fig, ax = plt.subplots(2, 2, figsize=(10.6, 7.2), layout="constrained")
    a = ax[0, 0]
    for (label, d), color in zip(km.groupby("band", sort=False), PALETTE):
        n = int(bands.loc[bands.band.eq(label), "n_all"].iloc[0])
        a.step(d.months, d.survival, where="post", color=color, lw=1.6, label=f"{label} months (n = {n})")
    a.set(xlim=(0, 60), ylim=(0, 1.02), xlabel="Months since first progression", ylabel="Survival after progression")
    a.legend(frameon=False, fontsize=7.6, title="Origin to progression", title_fontsize=7.6)
    title(a, "A   Post-progression survival, eight cohorts")
    a = ax[0, 1]
    a.axhline(1, color=GRAY, linestyle="--", linewidth=1)
    a.fill_between(spline.months, spline.lo, spline.hi, color=BLUE, alpha=.16)
    a.plot(spline.months, spline.hr, color=BLUE, lw=1.8)
    a.axvspan(0, 12, color=RED, alpha=.07)
    a.text(7.2, .22, "visible at\nmonth 12", ha="center", fontsize=7.4, color=RED)
    a.set(yscale="log", xlabel="Months from origin to first progression", xlim=(2, 48),
          ylabel="Adjusted hazard ratio for death\n(reference 12 months)")
    a.set_yticks([.2, .5, 1, 2, 4])
    a.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
    a.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    title(a, "B   Shape of the association, primary cohorts")
    full = cohort.loc[~cohort.adjusted_for.str.contains("within")]
    summary = pooled.loc[pooled.analysis.eq("Two-stage random effects, All progressions") & pooled.scope.eq("All cohorts")].iloc[0]
    labels = [f"{r.cohort} ({r.deaths}/{r.n})" for r in full.itertuples()] + ["Random-effects summary"]
    values = [(r.hr, r.lo, r.hi) for r in full.itertuples()] + [(summary.hr, summary.lo, summary.hi)]
    forest(ax[1, 0], labels, values, 1, [BLUE] * len(full) + [INK],
           "Hazard ratio per 6 months later progression", log=True, markers=["o"] * len(full) + ["D"])
    ax[1, 0].set_xticks([.25, .5, .75, 1])
    title(ax[1, 0], "C   Cohort-specific estimates (deaths/patients)")
    labels = list(conc.cohort)
    values = [(r.delta, r.lo, r.hi) for r in conc.itertuples()]
    forest(ax[1, 1], labels, values, 0, [BLUE] * 4 + [INK], "Change in concordance index (held-out cohort)",
           markers=["o"] * 4 + ["D"])
    title(ax[1, 1], "D   Discrimination gained from progression time")
    save(fig, "Figure_4_post_progression")


def molecular():
    hr, inc = read("molecular_pooled_hr"), read("molecular_increment_pooled")
    hall, brca = read("transcriptome_hallmarks"), read("brca")
    names = {"stroma": "Stromal", "tcell": "T-cell", "prolif": "Proliferation"}
    short = {"Progression-free survival from origin": "PFS",
             "Overall survival from origin": "OS",
             "Death within 24 months after the month-12 landmark, baseline and status adjusted": "Post-landmark death | status",
             "Post-progression survival, baseline and progression time adjusted": "Post-progression survival"}
    fig, ax = plt.subplots(2, 2, figsize=(10.8, 7.4), layout="constrained")
    labels, values, colors = [], [], []
    for (score, color) in zip(names, (GOLD, TEAL, BLUE)):
        for analysis, tag in short.items():
            r = hr.loc[hr.score.eq(score) & hr.analysis.eq(analysis)].iloc[0]
            labels.append(f"{names[score]}: {tag}")
            values.append((r.hr, r.lo, r.hi))
            colors.append(color)
    forest(ax[0, 0], labels, values, 1, colors, "Pooled hazard ratio per SD (eight cohorts)", log=True)
    ax[0, 0].set_xticks([.75, 1, 1.25])
    ax[0, 0].tick_params(axis="y", labelsize=7.8)
    title(ax[0, 0], "A   Expression scores and outcome")
    d = inc.loc[inc.reference.eq("Status")]
    forest(ax[0, 1], [("All three scores" if k.startswith("all") else names[k] + " score") for k in d.added],
           [(r.estimate, r.lo, r.hi) for r in d.itertuples()], 0, [GOLD, TEAL, BLUE, INK],
           "Pooled held-out change in AUC when added to Status", xlim=(-.05, .05))
    title(ax[0, 1], "B   Incremental discrimination at month 12")
    a = ax[1, 0]
    a.axhline(0, color=GRAY, lw=.8)
    a.axvline(0, color=GRAY, lw=.8)
    rep = hall.replicated.astype(bool)
    a.scatter(hall.NES_primary[~rep], hall.NES_additional[~rep], s=18, color="#B8C3CE")
    a.scatter(hall.NES_primary[rep], hall.NES_additional[rep], s=34, color=RED)
    for r in hall.loc[rep].itertuples():
        a.annotate(r.Term.replace("HALLMARK_", "").replace("_", " ").title(), (r.NES_primary, r.NES_additional),
                   xytext=(-7, -3), textcoords="offset points", fontsize=7.2, ha="right", va="center",
                   bbox={"facecolor": "white", "edgecolor": "none", "alpha": .85, "pad": .8})
    a.set(xlabel="Enrichment score, primary cohorts", ylabel="Enrichment score, additional cohorts", ylim=(-2.1, 3.4))
    title(a, "C   Hallmarks of progression within 12 months")
    keep = [("Progression-free survival from diagnosis", "brca", "PFS from diagnosis"),
            ("Landmark, death within 24 months", "brca", "Death ≤24 months after landmark"),
            ("Landmark, all follow-up", "brca", "Death after landmark, all follow-up"),
            ("Post-progression survival", "brca", "Post-progression survival"),
            ("Landmark, all follow-up", "early", "Progression status, BRCA adjusted")]
    rows = [brca.loc[brca.analysis.eq(a_) & brca.term.eq(t)].iloc[0] for a_, t, _ in keep]
    forest(ax[1, 1], [k[2] for k in keep], [(r.hr, r.lo, r.hi) for r in rows], 1, [TEAL] * 4 + [RED],
           "Hazard ratio (TCGA sequenced subset)", log=True)
    ax[1, 1].set_xticks([.33, .5, 1, 2, 4])
    ax[1, 1].tick_params(axis="y", labelsize=7.8)
    title(ax[1, 1], "D   BRCA1/2 mutation, TCGA")
    save(fig, "Figure_5_molecular")


if __name__ == "__main__":
    FIG.mkdir(exist_ok=True)
    performance()
    landmark()
    post_progression()
    molecular()


def genomic_singlecell():
    g, b = read("genomics"), read("brca")
    att, hall = read("sc_attribution"), read("transcriptome_hallmarks_stroma_adjusted")
    fig, ax = plt.subplots(2, 2, figsize=(11.2, 8.0), layout="constrained")
    short, long_ = "Death within 24 months after the month-12 landmark, status adjusted", \
        "Death after the month-12 landmark, all follow-up, status adjusted"
    labels, values, colors = [], [], []
    rows = [("BRCA1/2 mutation", b.loc[b.term.eq("brca") & b.analysis.eq("Landmark, death within 24 months")].iloc[0],
             b.loc[b.term.eq("brca") & b.analysis.eq("Landmark, all follow-up")].iloc[0], "hr")]
    for feature in ("HRD score ≥42", "CCNE1 amplification"):
        d = g.loc[g.feature.eq(feature)]
        rows.append((feature, d.loc[d.analysis.eq(short)].iloc[0], d.loc[d.analysis.eq(long_)].iloc[0], "estimate"))
    for name, first, second, key in rows:
        labels += [f"{name}: ≤24 months", f"{name}: all follow-up"]
        values += [(first[key], first.lo, first.hi), (second[key], second.lo, second.hi)]
        colors += [GOLD, TEAL]
    forest(ax[0, 0], labels, values, 1, colors, "Hazard ratio for death after the landmark, status adjusted", log=True)
    ax[0, 0].set_xticks([.5, .75, 1, 1.5, 2])
    ax[0, 0].tick_params(axis="y", labelsize=7.8)
    title(ax[0, 0], "A   Genomic features and death after month 12, TCGA")
    inc = [("BRCA1/2 mutation", b.loc[b.analysis.str.startswith("Out-of-fold")].iloc[0], "hr")]
    for feature in ("HRD score ≥42", "CCNE1 amplification", "Joint model: BRCA1/2, HRD score and CCNE1"):
        inc.append((feature.replace("Joint model: ", "").replace("BRCA1/2, HRD score and CCNE1", "All three features"),
                    g.loc[g.feature.eq(feature) & g.analysis.str.startswith("Out-of-fold")].iloc[0], "estimate"))
    forest(ax[0, 1], [k[0] for k in inc], [(r[key], r.lo, r.hi) for _, r, key in inc], 0, [TEAL] * 3 + [INK],
           "Out-of-fold change in AUC when added to Status", xlim=(-.07, .07))
    title(ax[0, 1], "B   Incremental discrimination at month 12")
    order = ["Fibroblast", "Mesothelial", "Malignant", "Myeloid", "T/NK", "B/plasma"]
    sets = ["Higher in early progression, top 200 genes", "Lower in early progression, top 200 genes",
            "Hallmark epithelial–mesenchymal transition", "Hallmark angiogenesis", "Hallmark UV response down",
            "Stromal score genes", "T-cell score genes"]
    spec = att.pivot(index="gene_set", columns="cell_class", values="mean_specificity").loc[sets, order]
    share = att.pivot(index="gene_set", columns="cell_class", values="share_highest").loc[sets, order]
    a = ax[1, 0]
    im = a.imshow(spec.to_numpy(), cmap="RdBu_r", vmin=-1.5, vmax=1.5, aspect="auto")
    for i in range(len(sets)):
        for j in range(len(order)):
            a.text(j, i, f"{100 * share.iloc[i, j]:.0f}%", ha="center", va="center", fontsize=7.4,
                   color="white" if abs(spec.iloc[i, j]) > .9 else INK)
    a.set_xticks(range(len(order)), order, rotation=30, ha="right")
    a.set_yticks(range(len(sets)), [s.replace(", top 200 genes", " (200)").replace("Hallmark ", "").capitalize()
                                    for s in sets], fontsize=7.8)
    for side in ("left", "bottom"):
        a.spines[side].set_visible(False)
    a.tick_params(length=0)
    fig.colorbar(im, ax=a, shrink=.75, label="Mean cell-class specificity (z)")
    title(a, "C   Cell classes expressing each gene set")
    keep = hall.loc[hall.replicated.astype(bool) | hall.adj_replicated.astype(bool)].sort_values("NES_primary", ascending=False)
    a = ax[1, 1]
    a.axvline(0, color=GRAY, linestyle="--", linewidth=1)
    y = np.arange(len(keep))
    a.hlines(y, keep.NES_primary, keep.adj_NES_primary, color="#C5CED8", lw=1.5)
    a.scatter(keep.NES_primary, y, color=RED, s=30, label="Unadjusted", zorder=3)
    a.scatter(keep.adj_NES_primary, y, color=BLUE, s=30, marker="D", label="Adjusted for stromal score", zorder=3)
    a.set(yticks=y, ylim=(len(keep) - .5, -.6), xlabel="Normalised enrichment score, primary cohorts")
    a.set_yticklabels([t.replace("HALLMARK_", "").replace("_", " ").capitalize() for t in keep.Term], fontsize=7.6)
    a.tick_params(axis="y", length=0)
    a.grid(axis="x", color="#E1E6EA", linewidth=.6)
    a.legend(frameon=False, fontsize=7.6, loc="upper center", bbox_to_anchor=(.45, -.13), ncol=2)
    title(a, "D   Hallmarks before and after stromal adjustment")
    save(fig, "Figure_6_genomic_singlecell")


if __name__ == "__main__":
    genomic_singlecell()

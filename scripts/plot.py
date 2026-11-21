"""Render original manuscript figures from the landmark analysis."""
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Circle
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

from analyze import (COHORTS, SEED, TIMING, FEATURE_LABELS, cohort_data,
                     fit_predict, horizon_data, log_model, algorithms)
from prepare_data import load_all

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)
INK, GRID, RED, BLUE, TEAL, GOLD = "#1B2638", "#DDE5ED", "#A63E42", "#365E9F", "#14786B", "#AD6718"
GROUPS = ["No progression by month 12", "Progression <6 months ago", "Progression ≥6 months ago"]
GROUP_COLORS = [TEAL, GOLD, RED]
MODEL_COLORS = {"Baseline": "#626D80", "Status": BLUE, "Timing": RED}
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10,
    "axes.titleweight": "bold", "axes.labelsize": 9, "axes.edgecolor": "#9AAABE",
    "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
    "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "savefig.facecolor": "white",
})


def read(name):
    return pd.read_csv(OUT / (name + ".csv"))


def panel(ax, title):
    ax.set_title(title, loc="left", pad=9)


def save(fig, filename):
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{filename}.{ext}", dpi=300, bbox_inches="tight", pad_inches=.12)
    plt.close(fig)


def spline_design(frame):
    t = frame.duration / 6
    spline = (np.maximum(t, 0) ** 3 - 2 * np.maximum(t - 1, 0) ** 3 + np.maximum(t - 2, 0) ** 3) / 4
    return frame.assign(curvature=spline * frame.early)


def figure_1(frames):
    fig = plt.figure(figsize=(10.5, 5.6), layout="constrained")
    grid = fig.add_gridspec(2, 3, width_ratios=[1.4, 1, 1])
    ax = fig.add_subplot(grid[:, 0])
    train = spline_design(frames["TCGA"])
    q = train.loc[train.early.eq(1), "duration"].quantile([.05, .95])
    time = np.linspace(*q, 95)
    scenario = spline_design(pd.DataFrame({
        "duration": np.r_[time, 0], "age10": 0, "iv": 0, "rd": 0,
        "early": np.r_[np.ones(len(time)), 0], "timing6": np.r_[(time - 6) / 6, 0],
    }))
    cols = TIMING + ["curvature"]
    risk = fit_predict(log_model(), train, {"profile": scenario}, cols, 12)["profile"]
    rng = np.random.default_rng(SEED)
    band = np.array([fit_predict(log_model(), train.iloc[rng.integers(0, len(train), len(train))],
                                 {"profile": scenario}, cols, 12)["profile"] for _ in range(350)])
    lo, hi = np.percentile(band, [2.5, 97.5], axis=0)
    ax.fill_between(time, lo[:-1], hi[:-1], color=RED, alpha=.15)
    ax.plot(time, risk[:-1], lw=2, color=RED, label="Earlier progression")
    ax.axhline(risk[-1], color=TEAL, lw=1.4, ls="--", label="No progression by landmark")
    ax.scatter(train.loc[train.early.eq(1), "duration"], np.repeat(.014, train.early.sum()),
               marker="|", s=22, color=RED, alpha=.23, clip_on=True)
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    panel(ax, "A   Timing and adjusted mortality")
    ax.set(xlim=(0, 12.1), ylim=(0, max(.55, hi[:-1].max() * 1.12)),
           xlabel="Months since first progression at landmark", ylabel="Predicted 12-month death risk")
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(1))
    for pos, (cohort, d) in enumerate(frames.items()):
        ax = fig.add_subplot(grid[pos // 2, 1 + pos % 2])
        for group, part in d.groupby("group"):
            km = KaplanMeierFitter().fit(part.time, part.death)
            times = km.survival_function_.index[km.survival_function_.index <= 24].to_numpy()
            times = np.r_[times, 24]
            ax.step(times, km.predict(times), where="post", color=GROUP_COLORS[group], lw=1.7,
                    label=f"{GROUPS[group]} (n={len(part)})")
            cens = part.loc[part.death.eq(0) & part.time.le(24), "time"]
            ax.plot(cens, km.predict(cens), ls="", marker="+", ms=3.4, color=GROUP_COLORS[group], alpha=.48)
        panel(ax, f"{chr(66 + pos)}   {cohort}")
        ax.set(xlim=(0, 24), ylim=(0, 1.03), xticks=[0, 6, 12, 18, 24],
               ylabel="Overall survival" if pos % 2 == 0 else None,
               xlabel="Months after landmark" if pos // 2 == 1 else None)
        ax.grid(axis="y", color=GRID, lw=.65)
        ax.legend(loc="lower left", fontsize=6.5, frameon=False, handlelength=1.25)
    save(fig, "Figure_1_survival")


def figure_2(frames):
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), layout="constrained",
                             gridspec_kw={"width_ratios": [1.12, 1]})
    ax = axes[0]
    panel(ax, "A   Information available at the landmark")
    for j, (name, first, color) in enumerate([("No progression recorded", None, TEAL),
                                             ("Progression at month 9", 9, GOLD),
                                             ("Progression at month 3", 3, RED)]):
        level = 3 - j
        ax.hlines(level, 0, 12, color="#9AAABE", lw=2)
        ax.hlines(level, 12, 24, color=BLUE, lw=1.7, ls=(0, (3, 2)))
        if first:
            ax.hlines(level, first, 12, color=color, lw=2.8)
            ax.scatter(first, level, s=38, color=color, zorder=4)
        ax.text(.05, level + .16, name, size=8, color=color)
    ax.axvline(12, color=INK, lw=.8)
    ax.text(12, 3.52, "Predict at month 12", ha="center", fontsize=8)
    ax.annotate("Ascertain death", (17.8, .64), color=BLUE, ha="center", fontsize=8)
    ax.set(xlabel="Months from the cohort-specific origin", xlim=(0, 24.2), ylim=(.45, 3.9),
           xticks=[0, 3, 6, 9, 12, 18, 24], yticks=[])
    ax.spines["left"].set_visible(False)
    ax = axes[1]
    panel(ax, "B   Observed 12-month mortality")
    for idx, (cohort, d) in enumerate(frames.items()):
        for early, color, off in [(0, TEAL, .14), (1, RED, -.14)]:
            group = d.loc[d.early.eq(early)]
            km = KaplanMeierFitter().fit(group.time, group.death)
            point = 1 - km.predict(12)
            row = km.confidence_interval_.loc[:12].iloc[-1]
            lo, hi = 1 - row.iloc[1], 1 - row.iloc[0]
            ax.errorbar(point, idx + off, xerr=[[point - lo], [hi - point]], fmt="o",
                        color=color, ms=5, capsize=2.5,
                        label=("No progression" if early == 0 else "Progression") if not idx else None)
    ax.set(xlim=(0, .9), ylim=(3.65, -.63), yticks=range(len(frames)), yticklabels=list(frames),
           xlabel="Risk of death within 12 months after landmark")
    ax.xaxis.set_major_formatter(ticker.PercentFormatter(1))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(.2))
    ax.grid(axis="x", color=GRID, lw=.6)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    save(fig, "Figure_2_landmark_evidence")


def figure_3(data):
    fig, ax = plt.subplots(2, 2, figsize=(10.5, 6.5), layout="constrained")
    cs = np.geomspace(.005, 20, 26)
    y, w, _, _ = horizon_data(data, 12)
    X = data[TIMING].to_numpy()
    transform = StandardScaler().fit(X)
    X_std = transform.transform(X)
    obs = w > 0
    strata = y + 2 * (~obs)
    splits = list(StratifiedKFold(5, shuffle=True, random_state=SEED).split(X, strata))

    def lasso(c):
        return LogisticRegression(C=c, l1_ratio=1, solver="saga", max_iter=4000, random_state=SEED)

    coefficients, briers = [], []
    for c in cs:
        coefficients.append(lasso(c).fit(X_std[obs], y[obs], sample_weight=w[obs]).coef_[0])
        losses = []
        for tr, ts in splits:
            yt, wt, _, _ = horizon_data(data.iloc[tr], 12)
            observed = wt > 0
            scale = StandardScaler().fit(X[tr])
            fit = lasso(c).fit(scale.transform(X[tr][observed]), yt[observed], sample_weight=wt[observed])
            p = fit.predict_proba(scale.transform(X[ts]))[:, 1]
            yn, wn, _, _ = horizon_data(data.iloc[ts], 12)
            losses.append(np.dot(wn, (yn - p) ** 2) / len(ts))
        briers.append(losses)
    coefficients, briers = np.asarray(coefficients), np.asarray(briers)
    colors = [BLUE, TEAL, GOLD, RED, "#79578B"]
    best = briers.mean(axis=1).argmin()
    lasso_set = set(np.flatnonzero(np.abs(coefficients[best]) > 1e-6))
    for column, color, coef in zip(TIMING, colors, coefficients.T, strict=True):
        ax[0, 0].plot(cs, coef, color=color, lw=1.65, label=FEATURE_LABELS[column])
    ax[0, 0].axvline(cs[best], color=INK, lw=.8, ls=":")
    ax[0, 0].legend(frameon=False, fontsize=6.7, loc="upper left")
    panel(ax[0, 0], "A   LASSO coefficient paths")
    ax[0, 0].set(xscale="log", xlabel="Inverse regularization strength C", ylabel="Standardized coefficient")
    mean = briers.mean(axis=1)
    sem = briers.std(axis=1, ddof=1) / np.sqrt(5)
    ax[0, 1].plot(cs, mean, color=BLUE, marker=".", ms=4, lw=1.5)
    ax[0, 1].fill_between(cs, mean - sem, mean + sem, color=BLUE, alpha=.15)
    ax[0, 1].axvline(cs[best], color=INK, lw=.8, ls=":")
    ax[0, 1].set(xscale="log", xlabel="Inverse regularization strength C", ylabel="Validation Brier score")
    panel(ax[0, 1], "B   Five-fold penalty comparison")
    rng = np.random.default_rng(SEED)
    imp = []
    for trial in range(90):
        orig = X_std[obs]
        shadow = np.column_stack([rng.permutation(c) for c in orig.T])
        model = RandomForestClassifier(n_estimators=90, min_samples_leaf=12, max_features="sqrt",
                                       random_state=SEED + trial, n_jobs=2)
        model.fit(np.column_stack([orig, shadow]), y[obs], sample_weight=w[obs])
        imp.append(model.feature_importances_)
    imp = np.asarray(imp)
    shadow_max = imp[:, len(TIMING):].max(axis=1)
    hits = (imp[:, :len(TIMING)] > shadow_max[:, None]).mean(axis=0)
    screen = set(np.flatnonzero(hits >= .9))
    labels = [FEATURE_LABELS[k] for k in TIMING] + ["Maximum shadow"]
    box = ax[1, 0].boxplot([*imp[:, :len(TIMING)].T, shadow_max], tick_labels=labels,
                            vert=False, patch_artist=True, showfliers=False, widths=.53)
    for item, color in zip(box["boxes"], [*colors, "#9BAABD"], strict=True):
        item.set(facecolor=color, alpha=.72, edgecolor=color)
    for median in box["medians"]:
        median.set(color=INK)
    ax[1, 0].tick_params(axis="y", labelsize=7.5)
    ax[1, 0].set_xlabel("Weighted forest impurity importance")
    panel(ax[1, 0], "C   Shadow-feature comparison")
    ax[1, 1].set(xlim=(0, 10), ylim=(0, 5.4), aspect="equal")
    ax[1, 1].axis("off")
    for x, col in [(3.75, BLUE), (6.25, TEAL)]:
        ax[1, 1].add_patch(Circle((x, 2.55), 2.05, facecolor=col, edgecolor=col, lw=1.1, alpha=.19))
    for x, label, set_ in [(2.55, "LASSO", lasso_set - screen),
                           (5, "Shared", lasso_set & screen),
                           (7.45, "Shadow", screen - lasso_set)]:
        ax[1, 1].text(x, 3.23, label, weight="bold", ha="center", fontsize=8)
        ax[1, 1].text(x, 2.34, "\n".join(FEATURE_LABELS[TIMING[i]].replace("Time since progression", "Progression timing")
                                             for i in sorted(set_)) or "—",
                      ha="center", va="center", fontsize=7)
    panel(ax[1, 1], "D   Diagnostic screening overlap")
    pd.DataFrame({"feature": TIMING, "lasso_coef": coefficients[best], "shadow_hit_rate": hits,
                  "c_min_brier": cs[best]}).to_csv(OUT / "screen_diagnostics.csv", index=False)
    save(fig, "Figure_3_predictor_diagnostics")


def figure_4(frames, predictions, performance):
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.7), layout="constrained")
    names = list(algorithms())
    palette = [BLUE, RED, TEAL, "#AA6918", "#714A84", "#4286A0", "#748A3D", "#9C5842"]
    for pos, (cohort, data) in enumerate(frames.items()):
        ax = axes.flat[pos]
        y, w, _, _ = horizon_data(data, 12)
        ax.plot([0, 1], [0, 1], color="#AAB5C0", ls="--", lw=1)
        for name, col in zip(names, palette, strict=True):
            p = predictions.loc[predictions.cohort.eq(cohort) & predictions.model.eq(name), "risk"]
            x, z, _ = roc_curve(y, p, sample_weight=w)
            auc = performance.loc[performance.cohort.eq(cohort) & performance.model.eq(name), "auc"].iloc[0]
            ax.plot(x, z, color=col, lw=1.9 if name == "Logistic" else 1.05,
                    label=f"{name}  {auc:.3f}", alpha=.92)
        panel(ax, f"{chr(65 + pos)}   {cohort}" + (" out of fold" if cohort == "TCGA" else " external"))
        ax.set(xlim=(0, 1.005), ylim=(0, 1.005), xlabel="False positive rate", ylabel="True positive rate")
        ax.grid(color=GRID, lw=.55)
        ax.legend(loc="lower right", frameon=False, fontsize=7.1, labelspacing=.27)
    save(fig, "Figure_4_roc")


def figure_5(frames, predictions):
    fig, axes = plt.subplots(2, 4, figsize=(11.2, 5.9), layout="constrained")
    for i, (cohort, data) in enumerate(frames.items()):
        ax, dc = axes[:, i]
        ax.plot([0, .75], [0, .75], c="#9BAABD", lw=.85, ls="--")
        y, w, _, _ = horizon_data(data, 12)
        th = np.linspace(.05, .45, 101)
        dc.axhline(0, c="#9BAABD", lw=.8, ls="--")
        treat_all = np.sum(w * y) / len(y) - np.sum(w * (1 - y)) / len(y) * th / (1 - th)
        dc.plot(th, treat_all, c="#778899", lw=.8, ls=":", label="Flag all")
        for model, color in MODEL_COLORS.items():
            p = predictions.loc[predictions.cohort.eq(cohort) & predictions.model.eq(model), "risk"].to_numpy()
            for_bo = data.assign(risk=p, quantile=pd.qcut(p, 4, duplicates="drop"))
            cal = []
            for _, part in for_bo.groupby("quantile", observed=True):
                km = KaplanMeierFitter().fit(part.time, part.death)
                cal.append((part.risk.mean(), 1 - km.predict(12)))
            x, z = np.array(cal).T
            ax.plot(x, z, "o-", c=color, lw=1, ms=3, label=model)
            action = p[:, None] >= th
            benefit = (w * y) @ action / len(y) - (w * (1 - y)) @ action / len(y) * th / (1 - th)
            dc.plot(th, benefit, c=color, lw=1.2, label=model)
        panel(ax, f"{chr(65 + i)}   {cohort}")
        panel(dc, f"{chr(69 + i)}   {cohort}")
        ax.set(xlim=(0, .75), ylim=(0, .75), xticks=[0, .25, .5, .75], yticks=[0, .25, .5, .75],
               xlabel="Mean predicted risk", ylabel="Observed risk" if i == 0 else None)
        dc.set(xlim=(.05, .45), ylim=(-.09, .25), xticks=[.1, .2, .3, .4],
               xlabel="Risk threshold", ylabel="Net benefit" if i == 0 else None)
        for panel_ax in (ax, dc):
            panel_ax.xaxis.set_major_formatter(ticker.PercentFormatter(1, decimals=0))
            panel_ax.grid(axis="y", color=GRID, lw=.55)
        ax.yaxis.set_major_formatter(ticker.PercentFormatter(1, decimals=0))
        if i == 0:
            ax.legend(fontsize=6.7, loc="upper left", frameon=False)
            dc.legend(fontsize=6.4, loc="upper right", frameon=False, ncol=1)
    save(fig, "Figure_5_calibration_and_decision_curves")


def figure_6(frames):
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 6.2), layout="constrained")
    cox = read("cox_associations")
    sub = read("subgroups")
    for idx, (cohort, d) in enumerate(frames.items()):
        ax = axes.flat[idx]
        summary = cox.loc[cox.cohort.eq(cohort) & cox.model.eq("Status") & cox.term.eq("early")].copy()
        summary = summary.assign(group="Overall", n=len(d))
        f = pd.concat([summary, sub.loc[sub.cohort.eq(cohort)]], ignore_index=True)
        for i, row in f.iterrows():
            if np.isfinite(row.hr):
                color = RED if i == 0 else BLUE
                lo, hi = max(.27, row.lo), min(39, row.hi)
                x = np.clip(row.hr, .27, 39)
                ax.errorbar(x, i, xerr=[[max(0, x - lo)], [max(0, hi - x)]],
                            c=color, fmt="D" if i == 0 else "o", ms=4.7, capsize=2)
                if row.hi > 39:
                    ax.annotate("", xy=(39.5, i), xytext=(31, i),
                                arrowprops={"arrowstyle": "->", "color": color, "lw": .8})
            else:
                ax.text(1.25, i, "Insufficient overlap", va="center", fontsize=7, color="#697386")
        ax.set(xscale="log", xlim=(.25, 40), ylim=(len(f) - .45, -.62), yticks=range(len(f)),
               yticklabels=[f"{r.group}   {r.events}/{r.n}" for r in f.itertuples()],
               xlabel="Adjusted hazard ratio (95% CI)")
        ax.set_xticks([.25, .5, 1, 2, 4, 8, 16, 32], labels=["0.25", "0.5", "1", "2", "4", "8", "16", "32"])
        ax.tick_params(axis="y", labelsize=7.5, length=0)
        ax.axvline(1, c="#9BAABD", lw=.9, ls="--")
        ax.grid(axis="x", color=GRID, lw=.6)
        panel(ax, f"{chr(65 + idx)}   {cohort}")
    save(fig, "Figure_6_subgroups")


def supplements(frames):
    coefficients = json.loads((OUT / "equations.json").read_text(encoding="utf-8"))["Status"]
    fig, ax = plt.subplots(figsize=(9.4, 4.5), layout="constrained")
    ax.set(xlim=(-67, 262), ylim=(-.54, 7.6))
    ax.axis("off")
    points = 50
    reference = coefficients["intercept"] - 2.5 * coefficients["age10"] + coefficients["iv"]

    def scale(y, name, x, labels, color=INK):
        ax.text(-66, y, name, size=9, va="center", weight="bold")
        ax.hlines(y, min(x), max(x), lw=1.05, color=color)
        ax.vlines(x, y - .12, y + .12, color=color, lw=.85)
        for value, label in zip(x, labels, strict=True):
            ax.text(value, y - .33, label, ha="center", va="top", fontsize=7.5, color=color)

    scale(7, "Points per predictor", list(range(0, 101, 20)), list(range(0, 101, 20)))
    age = np.arange(35, 86, 10)
    scale(6, "Age at diagnosis", points * coefficients["age10"] * (age - 35) / 10, age)
    scale(5, "FIGO stage", [0, -points * coefficients["iv"]], ["IV", "III"])
    scale(4, "Residual disease", [0, points * coefficients["rd"]], ["≤1 cm", ">1 cm"])
    scale(3, "Progression by month 12", [0, points * coefficients["early"]], ["No", "Yes"], RED)
    scale(1.65, "Total points", list(range(0, 261, 26)), list(range(0, 261, 26)), BLUE)
    probs = np.array([.02, .05, .1, .2, .3, .4, .5, .6])
    scale(.6, "Next 12-month death risk", points * (np.log(probs / (1 - probs)) - reference),
          [f"{100 * v:.0f}%" for v in probs], BLUE)
    ax.text(-66, -.46, "TCGA equation without external recalibration. Values are for research illustration.",
            fontsize=7.6, color="#59677B")
    panel(ax, "Nomogram for the progression status model")
    save(fig, "Figure_S1_nomogram")

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.7), layout="constrained", sharey=True)
    rng = np.random.default_rng(SEED + 1)
    for i, cohort in enumerate(COHORTS[1:]):
        ax = axes[i]
        source = spline_design(frames[cohort])
        duration = source.loc[source.early.eq(1), "duration"]
        grid = np.linspace(*duration.quantile([.05, .95]), 65)
        scenario = spline_design(pd.DataFrame({
            "age10": 0, "iv": 0, "rd": 0, "early": np.r_[np.ones(len(grid)), 0],
            "duration": np.r_[grid, 0], "timing6": np.r_[(grid - 6) / 6, 0],
        }))
        cols = TIMING + ["curvature"]
        fit = fit_predict(log_model(), source, {"scenario": scenario}, cols, 12)["scenario"]
        boot = np.array([
            fit_predict(log_model(), source.iloc[rng.integers(0, len(source), len(source))],
                        {"scenario": scenario}, cols, 12)["scenario"]
            for _ in range(250)
        ])
        lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
        ax.fill_between(grid, lo[:-1], hi[:-1], color=RED, alpha=.16)
        ax.plot(grid, fit[:-1], color=RED, lw=1.8)
        ax.axhline(fit[-1], color=TEAL, lw=1, ls="--")
        ax.scatter(duration, np.repeat(.014, len(duration)), marker="|", s=20, color=RED, alpha=.3)
        ax.set(xlim=(0, 12), ylim=(0, 1), xlabel="Months since progression")
        ax.yaxis.set_major_formatter(ticker.PercentFormatter(1, decimals=0))
        if i == 0:
            ax.set_ylabel("Predicted 12-month death risk")
        panel(ax, f"{chr(65 + i)}   {cohort}")
    save(fig, "Figure_S2_external_timing")

    c = read("model_contrasts")
    c = c.loc[c.cohort.ne("TCGA") & c.contrast.isin(["Status vs baseline", "Timing vs status"])]
    variants = list(c.variant.unique())
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 5.6), layout="constrained", sharey=True)
    for i, contrast in enumerate(c.contrast.unique()):
        ax = axes[i]
        part = c.loc[c.contrast.eq(contrast)]
        for k, (key, color) in enumerate(zip(["AOCS", "UKE", "ICON7", "Japan"], [TEAL, GOLD, BLUE, RED], strict=True)):
            cols = part.loc[part.cohort.eq(key)]
            y = [variants.index(v) + (k - 1.5) * .16 for v in cols.variant]
            ax.errorbar(cols.delta, y, xerr=[cols.delta - np.minimum(cols.delta, cols.lo),
                                            np.maximum(cols.delta, cols.hi) - cols.delta],
                        c=color, fmt="o", ms=4.2, lw=1.05, capsize=1.8, label=key)
        ax.axvline(0, c="#9BAABD", ls="--", lw=.85)
        ax.grid(axis="x", color=GRID, lw=.55)
        ax.set(yticks=np.arange(len(variants)), yticklabels=variants,
               ylim=(len(variants) - .5, -.51), xlabel="Change in time-dependent AUC")
        panel(ax, f"{chr(65 + i)}   {contrast}")
        ax.legend(loc="lower right", frameon=False, fontsize=7, ncol=2)
    save(fig, "Figure_S3_sensitivity")

    flow = read("flow")
    fig, ax = plt.subplots(figsize=(10.2, 3.5), layout="constrained")
    ax.axis("off")
    fields = ["Source patients", "Stage III–IV serous", "Treatment eligible", "Complete required fields",
              "Alive and observed at landmark", "Progression status known"]
    tbl = flow.pivot(index="step", columns="cohort", values="n").loc[list(flow.step.unique()), list(COHORTS)]
    for col, key in enumerate(COHORTS):
        ax.text(.42 + col * .15, .91, key, ha="center", weight="bold", fontsize=10)
    for i, label in enumerate(fields):
        y = .81 - .137 * i
        ax.axhspan(y - .064, y + .054, color="#EFF4F8" if i % 2 == 0 else "white")
        ax.text(.015, y, label, va="center", fontsize=9)
        for col, count in enumerate(tbl.iloc[i]):
            ax.text(.42 + col * .15, y, f"{count:,}", ha="center", va="center",
                    weight="bold" if i == len(fields) - 1 else "normal")
    panel(ax, "Eligibility at the month 12 landmark")
    save(fig, "Figure_S4_cohort_flow")


if __name__ == "__main__":
    frames, _ = cohort_data(load_all())
    preds = read("predictions").loc[lambda d: d.variant.eq("Primary")]
    performance = read("model_metrics").loc[lambda d: d.variant.eq("Primary")]
    figure_1(frames)
    figure_2(frames)
    figure_3(frames["TCGA"])
    figure_4(frames, preds, performance)
    figure_5(frames, preds)
    figure_6(frames)
    supplements(frames)

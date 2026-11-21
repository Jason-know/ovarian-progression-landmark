"""Plot the added benchmark and fixed-model calibration estimates."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
benchmark = pd.read_csv(OUT / "clinical_benchmark.csv")
apollo = pd.read_csv(OUT / "apollo_clinical_contrasts.csv").query(
    "subset == 'All sites' and contrast == 'Status vs Baseline'").iloc[0]
pathology = pd.read_csv(OUT / "tcga_stroma_oof.csv").iloc[0]
metrics = pd.read_csv(HERE.parent / "results" / "model_metrics.csv").query(
    "variant == 'Primary' and model == 'Status'").set_index("cohort")
external = pd.read_csv(OUT / "apollo_clinical_transport.csv").query(
    "subset == 'All sites' and model == 'Status'").iloc[0]

plt.rcParams.update({"font.family": "DejaVu Sans", "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "font.size": 9})
ink = "#203049"
blue, teal, red = "#325B9C", "#12756D", "#A6433D"
fig, axes = plt.subplots(1, 2, figsize=(10.1, 3.55), layout="constrained")


def point(ax, value, lo, hi, row, color):
    ax.errorbar(value, row, xerr=[[max(0, value - lo)], [max(0, hi - value)]], fmt="o",
                color=color, capsize=2.7, markersize=5, lw=1.45)


ax = axes[0]
for row, r in enumerate(benchmark.itertuples()):
    point(ax, r.delta_auc, r.delta_auc_lo, r.delta_auc_hi, row, blue)
point(ax, pathology.delta_auc, pathology.delta_auc_lo, pathology.delta_auc_hi, 3, teal)
ax.set(yticks=range(4), yticklabels=["AOCS (30)", "UKE (29)", "ICON7 (33)", "TCGA stroma (49)"],
       ylim=(3.65, -.63), xlabel="Change in time-dependent AUC", xlim=(-.05, .155))
ax.set_title("A   Performance of additions to simpler models", loc="left", fontweight="bold")
ax.axvline(0, ls="--", color="#9AA8B5", lw=.8)

ax = axes[1]
for row, r in enumerate(benchmark.itertuples()):
    point(ax, metrics.loc[r.cohort, "oe"], r.oe_added_lo, r.oe_added_hi, row, blue)
point(ax, external.oe, apollo.oe_added_lo, apollo.oe_added_hi, 3, red)
ax.set(yticks=range(4), yticklabels=["AOCS (30)", "UKE (29)", "ICON7 (33)", "APOLLO-2 (5)"],
       ylim=(3.65, -.63), xlabel="Observed-to-expected death risk ratio", xlim=(0, 3.05))
ax.set_title("B   Calibration of the locked Status model", loc="left", fontweight="bold")
ax.axvline(1, ls="--", color="#9AA8B5", lw=.8)
for ax in axes:
    ax.tick_params(axis="y", length=0, colors=ink)
    ax.grid(axis="x", alpha=.25)
for suffix in ("pdf", "png"):
    fig.savefig(HERE / f"Review_Evidence.{suffix}", dpi=300, bbox_inches="tight", pad_inches=.12)

"""Build revised figures and traceable H&E panels from source Aperio slides."""
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle
import numpy as np
import openslide
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA, RESULTS, FIG = HERE / "data", HERE / "results", HERE / "figures"
LEGACY = ROOT / "同行预审与病理补充" / "results"
BLUE, TEAL, RED, INK = "#315B9B", "#197A6F", "#AE4346", "#1B2639"
GRAY = "#9AA8B4"
plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "pdf.fonttype": 42, "text.color": INK})
COHORTS = ["AOCS", "UKE", "ICON7"]
SLIDES = [
    ("TCGA-09-0366-01A-01-BS1", "e9206264-5233-4f7d-8535-c1bdca572edc", 5, (5050, 4800)),
    ("TCGA-24-2033-01A-01-BS1", "9b83823f-bd88-4b26-92c8-9af5165ecc69", 35, (8380, 6500)),
]


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=300, bbox_inches="tight", pad_inches=.12)
    plt.close(fig)


def forest(ax, title, xlabel, labels, values, ref, colors=None, xlim=None):
    ax.axvline(ref, color=GRAY, linestyle="--", linewidth=1)
    for i, (value, lo, hi) in enumerate(values):
        ax.errorbar(value, i, xerr=[[max(0, value-lo)], [max(0, hi-value)]],
                    fmt="o", capsize=2.6, lw=1.4, ms=5, color=(colors or [BLUE]*len(labels))[i])
    ax.set(yticks=np.arange(len(labels)), yticklabels=labels,
           ylim=(len(labels)-.47, -.57), xlabel=xlabel, xlim=xlim)
    ax.set_title(title, loc="left", fontweight="bold", pad=7)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#E1E6EA", linewidth=.6)


def efficacy():
    contrasts = pd.read_csv(ROOT / "results/model_contrasts.csv")
    a = contrasts.loc[contrasts.variant.eq("Primary") & contrasts.contrast.eq("Status vs baseline")].set_index("cohort")
    b = pd.read_csv(LEGACY / "clinical_benchmark.csv").set_index("cohort")
    fixed = pd.read_csv(RESULTS / "calibration_and_two_stage.csv").set_index(["estimate", "cohort"])
    apollo = pd.read_csv(LEGACY / "apollo_clinical_contrasts.csv").query(
        "subset == 'All sites' and contrast == 'Status vs Baseline'").iloc[0]
    apm = pd.read_csv(LEGACY / "apollo_clinical_transport.csv").query(
        "subset == 'All sites' and model == 'Status'").iloc[0]

    fig, ax = plt.subplots(2, 2, figsize=(10.6, 5.9), layout="constrained")
    forest(ax[0, 0], "A   Status model versus clinical baseline", "Change in time-dependent AUC", COHORTS,
           [a.loc[c, ["delta", "lo", "hi"]] for c in COHORTS], 0, xlim=(-.05, .44))
    forest(ax[0, 1], "B   Status model versus progression alone", "Change in time-dependent AUC", COHORTS,
           [b.loc[c, ["delta_auc", "delta_auc_lo", "delta_auc_hi"]] for c in COHORTS], 0, xlim=(-.05, .16))
    forest(ax[1, 0], "C   Adding the interval since first progression", "Change in AUC with two-stage resampling", COHORTS,
           [fixed.loc[("timing_delta_auc_two_stage", c), ["value", "lo", "hi"]] for c in COHORTS], 0,
           xlim=(-.072, .037))
    values = [fixed.loc[("oe", c), ["value", "lo", "hi"]] for c in COHORTS]
    values.append((apm.oe, apollo.oe_added_lo, apollo.oe_added_hi))
    forest(ax[1, 1], "D   Locked status model calibration", "Observed-to-expected 12-month mortality", 
           COHORTS + ["APOLLO-2 (exploratory)"], values, 1, colors=[BLUE, BLUE, BLUE, TEAL], xlim=(0, 2.9))
    save(fig, "Figure_2_model_performance")


def slide_panels():
    manifest = json.loads((DATA / "TCGA_OV_slide_manifest.json").read_text(encoding="utf-8"))["data"]["hits"]
    data = pd.read_csv(DATA / "eligible_qc_slide_records.csv")
    source_list = []
    detail = []
    for name, file_id, stroma, center in SLIDES:
        entry = next(e for e in manifest if e["file_id"] == file_id)
        assert ((data.slide == name) & (data.file_id == file_id) & (data.stroma == stroma)).any()
        path = DATA / (name + ".svs")
        with path.open("rb") as handle:
            assert hashlib.file_digest(handle, "md5").hexdigest() == entry["md5sum"]
        with openslide.OpenSlide(str(path)) as slide:
            mpp_x = float(slide.properties["openslide.mpp-x"])
            mpp_y = float(slide.properties["openslide.mpp-y"])
            width, height = round(760/mpp_x), round(540/mpp_y)
            x, y = center[0] - width//2, center[1] - height//2
            assert x >= 0 and y >= 0 and x+width <= slide.dimensions[0] and y+height <= slide.dimensions[1]
            crop = slide.read_region((x, y), 0, (width, height)).convert("RGB")
            base = slide.get_thumbnail((840, 840)).convert("RGB")
            image_w, image_h = slide.dimensions
            crop.save(FIG / (name + "_field.png"))
            detail.append((name, stroma, mpp_x, crop, base,
                           (x/image_w, y/image_h, width/image_w, height/image_h)))
            source_list.append({"slide": name, "case": name[:12], "file_id": file_id,
                                "file_md5": entry["md5sum"], "file_page": f"https://portal.gdc.cancer.gov/files/{file_id}",
                                "stroma_pct_slide": stroma, "mpp_x": mpp_x, "mpp_y": mpp_y, "x": x,
                                "y": y, "width_px": width, "height_px": height})
    pd.DataFrame(source_list).to_csv(RESULTS / "slide_provenance.csv", index=False)
    return detail


def annotate_scale(ax, shape, length):
    h, w = shape[:2]
    x, y = w*.055, h*.91
    ax.plot([x, x+length], [y, y], color="white", linewidth=5, solid_capstyle="butt")
    ax.plot([x, x+length], [y, y], color="#19202B", linewidth=2.2, solid_capstyle="butt")
    ax.text(x+length/2, y-h*.038, "100 µm", ha="center", va="bottom", fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": .88, "pad": 1.6})


def pathology():
    detail = slide_panels()
    fig, axes = plt.subplots(2, 2, figsize=(10.1, 7.4), layout="constrained",
                             gridspec_kw={"height_ratios": [1.2, .74]})
    for i, (name, stroma, mx, photo, _, _) in enumerate(detail):
        ax = axes[0, i]
        ax.imshow(photo)
        ax.axis("off")
        ax.set_title(f"{'AB'[i]}   TCGA tissue slide  |  QC stroma {stroma}%", loc="left", fontweight="bold")
        annotate_scale(ax, np.asarray(photo).shape, 100/mx)
    r = pd.read_csv(LEGACY / "tcga_stroma_oof.csv").iloc[0]
    forest(axes[1, 0], "C   Incremental discrimination from TCGA stroma",
           "Change in out-of-fold time-dependent AUC", ["TCGA"],
           [(r.delta_auc, r.delta_auc_lo, r.delta_auc_hi)], 0, xlim=(-.0065, .0065), colors=[TEAL])
    axes[1, 0].xaxis.set_major_formatter(ticker.FormatStrFormatter("%.3f"))
    axes[1, 0].text(.02, .07, "342 patients, 49 deaths  |  AUC 0.7721 to 0.7713",
                    transform=axes[1, 0].transAxes, size=7.6)
    ax = axes[1, 1]
    risk = pd.read_csv(LEGACY / "pathology_subset_landmark_risk.csv")
    labels = ["TCGA with stroma", "APOLLO-2 primary site"]
    for i, key in enumerate(risk.cohort.unique()):
        for r in risk.loc[risk.cohort.eq(key)].itertuples():
            color = RED if r.early else TEAL
            ax.errorbar(r.risk12, i+(.11 if r.early == 0 else -.11),
                        xerr=[[r.risk12-r.lower], [r.upper-r.risk12]], color=color, fmt="o", capsize=2,
                        label=("Progression" if r.early else "No progression") if i == 0 else None)
    ax.set(yticks=[0, 1], yticklabels=labels, ylim=(1.45, -.53), xlim=(0, .87),
           xlabel="Observed 12-month mortality")
    ax.xaxis.set_major_formatter(ticker.PercentFormatter(1))
    ax.set_title("D   Observed mortality in pathology subsets", loc="left", fontweight="bold")
    ax.grid(axis="x", color="#E1E6EA", linewidth=.6)
    ax.legend(loc="upper right", frameon=False, fontsize=7.5)
    save(fig, "Figure_3_pathology")

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.5), layout="constrained")
    for i, (name, stroma, _, _, base, coordinates) in enumerate(detail):
        ax = axes[i]
        ax.imshow(base)
        bw, bh = base.size
        x, y, w, h = coordinates
        ax.add_patch(Rectangle((x*bw, y*bh), w*bw, h*bh, fill=False, edgecolor="#087888", lw=2))
        ax.axis("off")
        ax.set_title(f"{'AB'[i]}   {name}  |  QC stroma {stroma}%", loc="left", fontweight="bold")
    save(fig, "Figure_S10_whole_slide_overviews")


if __name__ == "__main__":
    FIG.mkdir(exist_ok=True)
    efficacy()
    pathology()

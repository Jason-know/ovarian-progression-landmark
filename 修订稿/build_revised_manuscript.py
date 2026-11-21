"""Render the revised article while leaving the first draft intact."""
from pathlib import Path
import csv
import sys

from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OLD = ROOT / "results"
NEW = HERE / "results"
PEER = ROOT / "同行预审与病理补充" / "results"
sys.path.insert(0, str(ROOT / "scripts"))
import make_manuscript as layout

DOC = layout.DOC
COHORTS = ["TCGA", "AOCS", "UKE", "ICON7"]
EXTERNAL = COHORTS[1:]


def rows(folder, file):
    with (folder / f"{file}.csv").open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def fmt(s, decimals=3):
    return f"{float(s):.{decimals}f}".replace("-", "−") if s else "—"


def interval(value, lo, hi, decimals=3):
    return f"{fmt(value, decimals)} ({fmt(lo, decimals)} to {fmt(hi, decimals)})"


def content():
    layout.setup()
    DOC.core_properties.title = "Incremental Value of Progression Timing for Mortality Prediction in Advanced Serous Tubo-Ovarian Carcinoma"
    DOC.core_properties.subject = "Landmark prediction, external calibration and tissue pathology"
    DOC.styles["Source Reference"].paragraph_format.space_after = Pt(3)
    in_refs = False
    for line in (HERE / "revised_manuscript.txt").read_text(encoding="utf-8").splitlines():
        kind, text = line.split("|", 1)
        if kind == "REF":
            if not in_refs:
                DOC.add_heading("References", 1)
                in_refs = True
            body, url = text.rsplit(" https://", 1)
            p = layout.paragraph(body + " ", "Source Reference")
            label = "doi:" + url.removeprefix("doi.org/") if url.startswith("doi.org/") else "Source"
            layout.hyperlink(p, "https://" + url, label)
        elif kind.startswith("H"):
            DOC.add_heading(text, int(kind[1:]))
        elif kind in ("TITLE", "SUBTITLE"):
            DOC.add_paragraph(text, style=kind.title())
        else:
            layout.paragraph(text)


def main_tables():
    DOC.add_page_break()
    DOC.add_heading("Main tables", 1)
    characteristics = {r["cohort"]: r for r in rows(OLD, "cohort_characteristics") if r["group"] == "All"}
    clinical = [
        ("TCGA", "Diagnosis", "CDR PFI. Recurrence, progression and other new tumor events or death with tumor.",
         "Historical grade. Residual diameter from Firehose. No harmonized treatment filter."),
        ("AOCS", "Diagnosis", "Relapse time and vital status from the original annotations.",
         "Historical grade. Primary surgery without recorded NACT and with platinum."),
        ("UKE", "Diagnosis", "PFS and OS from the published patient-level archive.",
         "Serous label, no available grade. Platinum with or without bevacizumab."),
        ("ICON7", "Randomization", "Trial PFS and OS. Clinical or radiographic progression. Scheduled assessments.",
         "AGO-OVAR11 tissue subset. Source grade. Standard chemotherapy or bevacizumab arm."),
    ]
    table = []
    for name, origin, event, histology in clinical:
        r = characteristics[name]
        table.append([name, f'{r["n"]} / {r["early"]} / {r["dead12"]}', origin, event, histology])
    layout.table("Table 1 Primary cohort clocks and endpoint ascertainment",
                 ["Cohort", "N / progressed / deaths", "Origin", "Progression record", "Case mix and treatment"],
                 table, [.65, .94, .94, 2.06, 2.12],
                 "Progressed means recorded by the month-12 landmark. Deaths refer to the subsequent 12-month horizon.",
                 left_columns=5)
    m = {(r["cohort"], r["model"]): r for r in rows(OLD, "model_metrics") if r["variant"] == "Primary"}
    benchmark = {r["cohort"]: r for r in rows(PEER, "clinical_benchmark")}
    tcga = next(r for r in rows(OLD, "progression_only_benchmark") if r["cohort"] == "TCGA")
    contrasts = {r["cohort"]: r for r in rows(NEW, "calibration_and_two_stage") if r["estimate"] == "timing_delta_auc_two_stage"}
    table = []
    for name in COHORTS:
        alone = tcga["auc"] if name == "TCGA" else benchmark[name]["auc_reference"]
        status_only_delta = "—" if name == "TCGA" else interval(
            *(benchmark[name][k] for k in ("delta_auc", "delta_auc_lo", "delta_auc_hi")))
        timing_delta = fmt(float(m[name, "Timing"]["auc"]) - float(m[name, "Status"]["auc"])) if name == "TCGA" else interval(
            *(contrasts[name][k] for k in ("value", "lo", "hi")))
        table.append([name, fmt(m[name, "Baseline"]["auc"]), fmt(alone), fmt(m[name, "Status"]["auc"]),
                      fmt(m[name, "Timing"]["auc"]), status_only_delta, timing_delta])
    layout.table("Table 2 Discrimination with four predictor sets",
                 ["Cohort", "Baseline", "Progression only", "Status", "Timing",
                  "Status minus progression only", "Timing minus Status"],
                 table, [.64, .72, .78, .63, .63, 1.66, 1.64],
                 "Values are fixed-horizon AUCs. TCGA results are out of fold. External contrasts contain 95% "
                 "confidence intervals. The Timing contrast uses two-stage resampling.")


MAIN = [
    ("1", ROOT / "figures/Figure_2_landmark_evidence.png", "Landmark design and observed mortality",
     "Panel A shows the observation and outcome windows for three example progression histories. Panel B shows "
     "Kaplan–Meier mortality during the next 12 months for patients with and without progression at the landmark."),
    ("2", HERE / "figures/Figure_2_model_performance.png", "Additional predictors and external calibration",
     "A shows the gain from adding progression status to the clinical baseline. B uses progression alone as the "
     "benchmark. C adds timing and resamples both development and evaluation patients. D shows observed-to-expected "
     "risk for the locked Status model. APOLLO-2 is exploratory. All error bars indicate 95% bootstrap intervals."),
    ("3", HERE / "figures/Figure_3_pathology.png", "TCGA H&E fields and tissue pathology analyses",
     "A and B are unnormalized fields from TCGA-09-0366-01A-01-BS1 and TCGA-24-2033-01A-01-BS1, respectively. "
     "Percentages are source slide QC estimates, not estimates from the selected fields. Scanner metadata determined "
     "100 µm scale bars. C shows the AUC change after adding patient-median stromal percentage in matched "
     "out-of-fold comparisons (95% evaluation-bootstrap CI). D shows observed mortality by progression status "
     "in the TCGA tissue subset and the APOLLO-2 primary-site subset (pointwise 95% Kaplan–Meier CIs). "
     "Full-slide overviews, exact field coordinates and GDC file IDs accompany the analysis."),
]


def add_figure(number, path, title, caption, width=6.67):
    DOC.add_heading(f"Figure {number} {title}", 2)
    p = DOC.add_paragraph()
    p.paragraph_format.keep_with_next = True
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    pic = p.add_run().add_picture(str(path), width=Inches(width))
    pic._inline.docPr.set("descr", title + ". " + caption)
    layout.paragraph(caption, "Figure Caption")


def figures():
    DOC.add_page_break()
    DOC.add_heading("Main figures", 1)
    for i, item in enumerate(MAIN):
        if i == 2:
            DOC.add_page_break()
        add_figure(*item)
    DOC.add_page_break()
    DOC.add_heading("Supplementary figures", 1)
    prior = {num: (name, title, desc) for num, name, title, desc in layout.FIGURES}
    keep = ["1", "3", "4", "5", "6", "S1", "S2", "S3", "S4"]
    for i, key in enumerate(keep, 1):
        if i in (2, 3, 4, 6, 8):
            DOC.add_page_break()
        name, title, desc = prior[key]
        add_figure(f"S{i}", ROOT / f"figures/{name}.png", title, desc)
    DOC.add_page_break()
    add_figure("S10", HERE / "figures/Figure_S10_whole_slide_overviews.png", "Source whole-slide overviews",
               "Outline boxes identify the regions displayed in Figure 3. No slide-level color transformation "
               "other than downsampling for display was applied.")
    s = rows(NEW, "slide_provenance")
    layout.table("Table S1 H&E source files and crop coordinates",
                 ["Slide barcode", "Stroma QC", "Source file UUID", "Level 0 origin and size"],
                 [[r["slide"], r["stroma_pct_slide"] + "%",
                   r["file_id"], f'{r["x"]}, {r["y"]} | {r["width_px"]} × {r["height_px"]} px'] for r in s],
                 [1.79, .74, 2.56, 1.61],
                 "Image-file checksum and both microns-per-pixel values are in slide_provenance.csv.",
                 left_columns=4)


def supplementary_tables():
    DOC.add_page_break()
    DOC.add_heading("Supplementary validation tables", 1)
    ext = rows(NEW, "additional_validation")
    apollo = rows(PEER, "apollo_clinical_transport")
    subset = [r for r in apollo if r["subset"] == "All sites"]
    japan = [r for r in rows(OLD, "model_metrics") if r["variant"] == "No age" and r["cohort"] == "Japan"]
    all_out = {"APOLLO-2": (subset, "Full clinical model"),
               "UVA": ([r for r in ext if r["cohort"] == "UVA"], "Full clinical model"),
               "Curie": ([r for r in ext if r["cohort"] == "Curie"], "Stage-only baseline"),
               "OVCAD": ([r for r in ext if r["cohort"] == "OVCAD"], "Age and stage baseline"),
               "Japan": (japan, "Stage and residual baseline")}
    results = []
    for name, (dataset, features) in all_out.items():
        m = {r["model"]: r for r in dataset}
        status = m["Status"]
        results.append([name, features, f'{status["n"]} / {status["deaths"]}',
                        *(fmt(m[k]["auc"]) if k in m else "—" for k in ("Baseline", "Progression only", "Status", "Timing")),
                        fmt(status["oe"], 2)])
    layout.table("Table S2 Separate additional cohort evaluations",
                 ["Cohort", "Predictor availability", "N / deaths", "Baseline", "Progression only", "Status", "Timing", "Status O/E"],
                 results, [.8, 1.47, .78, .72, .8, .67, .67, .79],
                 "Each row uses the source-specific predictor set. Estimates are exploratory. "
                 "OVCAD includes recruitment at Hamburg and possible overlap with UKE. N is the month-12 risk set.")
    c = {(r["cohort"], r["estimate"]): r for r in rows(NEW, "calibration_and_two_stage")}
    layout.table("Table S3 External calibration of the locked Status model",
                 ["Cohort", "Intercept", "Slope", "Observed to expected ratio"],
                 [[name, *(interval(*(c[name, kind][k] for k in ("value", "lo", "hi")), 2)
                   for kind in ("calibration_intercept", "calibration_slope", "oe"))] for name in EXTERNAL],
                 [1, 1.9, 1.9, 1.9],
                 "Intercept is estimated with slope fixed at one. Parentheses contain 95% fixed-prediction bootstrap intervals.")
    layout.table("Table S4 Timing sensitivity with nonlinear specification",
                 ["Cohort", "Linear timing ΔAUC", "Spline timing ΔAUC"],
                 [[name, *(interval(*(c[name, kind][k] for k in ("value", "lo", "hi")))
                   for kind in ("timing_delta_auc_two_stage", "spline_delta_auc_two_stage"))] for name in EXTERNAL],
                 [1.16, 2.77, 2.77],
                 "Intervals use 600 joint development and evaluation bootstrap replicates. Both comparisons use Status as the reference.")
    DOC.add_page_break()
    b = {r["cohort"]: r for r in rows(PEER, "clinical_benchmark")}
    prior = {(r["cohort"], r["model"]): r for r in rows(OLD, "model_metrics") if r["variant"] == "Primary"}
    alone = {r["cohort"]: r for r in rows(OLD, "progression_only_benchmark")}
    layout.table("Table S5 Brier scores and paired differences",
                 ["Cohort", "Baseline", "Progression only", "Status", "Timing", "Status minus progression only"],
                 [[name, fmt(prior[name, "Baseline"]["brier"]), fmt(alone[name]["brier"]),
                   fmt(prior[name, "Status"]["brier"]), fmt(prior[name, "Timing"]["brier"]),
                   interval(*(b[name][k] for k in ("delta_brier", "delta_brier_lo", "delta_brier_hi")), 4)
                   if name in b else "—"] for name in COHORTS],
                 [.74, .93, 1.01, .83, .82, 2.37],
                 "Lower scores indicate smaller squared probability errors. Parentheses contain paired 95% confidence intervals.")
    s = {r["cohort"]: r for r in rows(OLD, "cohort_characteristics") if r["group"] == "All"}
    layout.table("Table S6 Baseline characteristics by primary cohort",
                 ["Cohort", "Age median (IQR)", "Stage IV", "Suboptimal cytoreduction",
                  "Median post-landmark follow-up"],
                 [[name, f'{fmt(s[name]["age_median"], 1)} ({fmt(s[name]["age_q1"], 1)}–{fmt(s[name]["age_q3"], 1)})',
                   s[name]["stage_iv"], s[name]["suboptimal"], fmt(s[name]["followup_median"], 1)]
                  for name in COHORTS],
                 [.8, 1.7, .9, 1.6, 1.7],
                 "Age is in years and follow-up is in months. Median follow-up uses reverse Kaplan–Meier.")


if __name__ == "__main__":
    content()
    main_tables()
    figures()
    supplementary_tables()
    DOC.save(HERE / "Ovarian_Cancer_Revised_Manuscript.docx")

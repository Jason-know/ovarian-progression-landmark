"""Render the second-revision manuscript. Earlier drafts are left untouched."""
from pathlib import Path
import csv
import sys

from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OLD, REV, PEER, NEW = ROOT / "results", ROOT / "修订稿" / "results", ROOT / "同行预审与病理补充" / "results", HERE / "results"
sys.path.insert(0, str(ROOT / "scripts"))
import make_manuscript as layout  # noqa: E402

DOC = layout.DOC
COHORTS = ["TCGA", "AOCS", "UKE", "ICON7"]
MINUS = "−"


def rows(folder, name):
    with (folder / f"{name}.csv").open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def f(value, d=3):
    if value in ("", None):
        return "—"
    text = f"{float(value):.{d}f}"
    if float(text) == 0:
        text = text.lstrip("-")
    return text.replace("-", MINUS)


def ci(v, lo, hi, d=3, sep=" "):
    return f"{f(v, d)}{sep}({f(lo, d)} to {f(hi, d)})"


def pval(p):
    p = float(p)
    return "<0.001" if p < .001 else f"{p:.3f}"


def content():
    layout.setup()
    DOC.core_properties.title = ("Progression Status and Progression Timing as Updated Prognostic Information "
                                 "in Advanced Serous Tubo-Ovarian Carcinoma")
    DOC.core_properties.subject = "Multi-cohort landmark, multistate and transcriptomic analysis"
    DOC.core_properties.keywords = "serous tubo-ovarian carcinoma, progression, landmark, multistate, external validation"
    DOC.styles["Source Reference"].paragraph_format.space_after = Pt(3)
    in_refs = False
    for line in (HERE / "v2_manuscript.txt").read_text(encoding="utf-8").splitlines():
        kind, text = line.split("|", 1)
        if kind == "REF":
            if not in_refs:
                DOC.add_heading("References", 1)
                in_refs = True
            body, url = text.rsplit(" https://", 1)
            p = layout.paragraph(body + " ", "Source Reference")
            label = "doi:" + url.removeprefix("doi.org/") if url.startswith("doi.org/") else "https://" + url
            layout.hyperlink(p, "https://" + url, label)
        elif kind.startswith("H"):
            DOC.add_heading(text, int(kind[1:]))
        elif kind in ("TITLE", "SUBTITLE"):
            DOC.add_paragraph(text, style=kind.title())
        else:
            layout.paragraph(text)


def main_tables():
    DOC.add_page_break()
    DOC.add_heading("Tables", 1)
    ch = {r["cohort"]: r for r in rows(OLD, "cohort_characteristics") if r["group"] == "All"}
    clinical = [
        ("TCGA", "Diagnosis", "CDR progression-free interval: recurrence, progression, other new tumour event or death with tumour.",
         "Historical three-tier grade. Residual diameter from Firehose. No treatment filter available."),
        ("AOCS", "Diagnosis", "Relapse time and vital status from the original clinical annotation.",
         "Historical three-tier grade. Primary surgery, no recorded neoadjuvant therapy, platinum given."),
        ("UKE", "Diagnosis", "PFS and OS from the published patient-level table.",
         "Serous label, no grade. Platinum with or without bevacizumab."),
        ("ICON7", "Randomisation", "Trial PFS and OS. Clinical or radiological progression at scheduled assessments.",
         "AGO-OVAR11 tissue subset. Two-level source grade. Chemotherapy with or without bevacizumab."),
    ]
    layout.table("Table 1 Primary cohorts at the month-12 landmark",
                 ["Cohort", "N / progressed / deaths", "Time origin", "Progression record", "Case mix and treatment"],
                 [[c, f'{ch[c]["n"]} / {ch[c]["early"]} / {ch[c]["dead12"]}', o, e, h] for c, o, e, h in clinical],
                 [.65, .94, .94, 2.06, 2.12],
                 "Progressed means progression recorded by the landmark. Deaths are those within the following 12 months.",
                 left_columns=5)
    c = {(r["design"], r["cohort"]): r for r in rows(NEW, "iecv_cohort")}
    p = {(r["design"], r["measure"]): r for r in rows(NEW, "iecv_pooled")}
    body = []
    for k in COHORTS:
        r = c["IECV", k]
        lock = c.get(("TCGA only (locked)", k))
        body.append([k, f'{r["n"]} / {r["deaths"]}', f(r["auc_Baseline"]), f(r["auc_Progression only"]), f(r["auc_Status"]),
                     ci(r["d_status_baseline"], r["d_status_baseline_lo"], r["d_status_baseline_hi"], sep="\n"),
                     ci(r["d_timing_status"], r["d_timing_status_lo"], r["d_timing_status_hi"], sep="\n"),
                     ci(r["oe"], r["oe_lo"], r["oe_hi"], 2, "\n"),
                     ci(lock["oe"], lock["oe_lo"], lock["oe_hi"], 2, "\n") if lock else "Development"])
    g = lambda m: p["IECV", m]  # noqa: E731
    body.append(["Summary", "850 / 142", f(g("AUC Baseline")["estimate"]), "—", f(g("AUC Status")["estimate"]),
                 *(ci(g(m)["estimate"], g(m)["lo"], g(m)["hi"], sep="\n") for m in
                   ("ΔAUC Status − Baseline", "ΔAUC Timing − Status")),
                 ci(g("O/E Status")["estimate"], g("O/E Status")["lo"], g("O/E Status")["hi"], 2, "\n"), "—"])
    layout.table("Table 2 Prediction of death within 12 months of the month-12 landmark",
                 ["Held-out cohort", "N / deaths", "Baseline AUC", "Progression-only AUC", "Status AUC",
                  "Status minus Baseline", "Timing minus Status", "Status O/E, held out", "Status O/E, TCGA-only equation"],
                 body, [.64, .58, .68, .85, .5, .86, .86, .82, .85],
                 "Cohort rows are internal–external estimates: the models were fitted in the other three cohorts. "
                 "Parentheses contain 95% bootstrap intervals. The summary row gives random-effects estimates with "
                 "Hartung–Knapp 95% confidence intervals. The last column applies the prespecified equation fitted in TCGA alone.")
    DOC.add_page_break()
    two = {r["cohort"]: r for r in rows(REV, "calibration_and_two_stage") if r["estimate"] == "timing_delta_auc_two_stage"}
    lm = {(int(r["landmark"]), r["contrast"]): r for r in rows(NEW, "landmark_iecv_pooled")}
    sm = next(r for r in rows(NEW, "supermodel_hr") if r["model"].startswith("Supermodel, landmarks 6") and r["term"] == "timing6")
    ill = next(r for r in rows(NEW, "illness_death") if r["term"] == "prog_ttp")
    pp = {(r["analysis"], r["scope"]): r for r in rows(NEW, "pps_pooled")}
    conc = next(r for r in rows(NEW, "pps_concordance") if r["cohort"].startswith("Random"))
    one = pp["One-stage stratified Cox, age, stage and residual adjusted", "Primary cohorts"]
    allc = pp["Two-stage random effects, All progressions", "All cohorts"]
    short = pp["Two-stage random effects, Progression within 12 months", "All cohorts"]
    t3 = [
        ["Month-12 landmark, TCGA-only equation (prespecified)", "AOCS, UKE, ICON7", "503 / 92", "Change in AUC",
         "; ".join(f'{k} {ci(two[k]["value"], two[k]["lo"], two[k]["hi"])}' for k in COHORTS[1:])],
        *[[f"Month-{m} landmark, internal–external", "Four primary cohorts", f'{lm[m, "timing"]["n"]} / {lm[m, "timing"]["deaths"]}',
           "Change in AUC", ci(lm[m, "timing"]["estimate"], lm[m, "timing"]["lo"], lm[m, "timing"]["hi"])] for m in (12, 18, 24)],
        ["Landmark supermodel, 6 to 24 months", "Four primary cohorts", f'{sm["n"]} / {sm["deaths"]} landmark-specific deaths',
         "HR per 6 months longer since progression", ci(sm["hr"], sm["lo"], sm["hi"], 2)],
        ["Illness–death model", "Four primary cohorts", f'{ill["patients"]} / {ill["deaths_after_progression"]} deaths after progression',
         "HR per 6 months later progression", ci(ill["hr"], ill["lo"], ill["hi"], 2)],
        ["Post-progression survival, one-stage", "Four primary cohorts", f'{int(float(one["n"]))} / {int(float(one["deaths"]))}',
         "HR per 6 months later progression", ci(one["hr"], one["lo"], one["hi"], 2)],
        ["Post-progression survival, two-stage", "Eight cohorts", f'{int(float(allc["n"]))} / {int(float(allc["deaths"]))}',
         "HR per 6 months later progression", ci(allc["hr"], allc["lo"], allc["hi"], 2)],
        ["Post-progression survival, progression within 12 months", "Eight cohorts",
         f'{int(float(short["n"]))} / {int(float(short["deaths"]))}', "HR per 6 months later progression",
         ci(short["hr"], short["lo"], short["hi"], 2)],
        ["Post-progression survival, held-out cohorts", "Four primary cohorts", f'{conc["n"]} / {conc["deaths"]}',
         "Change in concordance index", ci(conc["delta"], conc["lo"], conc["hi"])],
    ]
    layout.table("Table 3 Progression timing across analysis designs",
                 ["Design", "Cohorts", "Patients / deaths", "Measure", "Estimate (95% CI)"], t3,
                 [1.95, 1.05, 1.15, 1.2, 1.42],
                 "HR, hazard ratio. At a landmark, a longer time since progression means earlier progression, so hazard ratios "
                 "above one in the supermodel and below one in the post-progression analyses point in the same direction. "
                 "The first row uses two-stage resampling of development and evaluation patients. The two-stage "
                 "post-progression models use cohort-specific complete cases.", left_columns=4)


FIGURES = [
    ("1", ROOT / "figures/Figure_2_landmark_evidence.png", "Landmark design and observed mortality",
     "A shows three illustrative histories before the month-12 landmark and the 12-month window in which deaths were "
     "counted. B shows Kaplan–Meier 12-month mortality with 95% intervals by recorded progression at the landmark."),
    ("2", HERE / "figures/Figure_2_performance.png", "Discrimination and calibration at the month-12 landmark",
     "A and B show changes in AUC for each held-out cohort in internal–external cross-validation, with the random-effects "
     "summary. C shows observed-to-expected 12-month mortality for the Status model fitted without the held-out cohort, "
     "for the prespecified TCGA-only equation, and for the TCGA-only equation after a cross-validated intercept update. "
     "D shows smoothed calibration curves for the held-out predictions. Bars are 95% intervals."),
    ("3", HERE / "figures/Figure_3_landmarks.png", "Progression timing across landmarks",
     "A shows the distribution of months since first progression among progressed patients at each landmark in the four "
     "primary cohorts. Numbers give progressed patients and, in parentheses, those who progressed more than six months "
     "earlier. The dotted line marks six months. B shows the hazard ratio for death within 12 months per six months longer "
     "since progression at each landmark and from the supermodel. The month-6 estimate is omitted because no patient had "
     "progressed more than six months earlier. C and D show pooled held-out changes in AUC. Bars are 95% intervals."),
    ("4", HERE / "figures/Figure_4_post_progression.png", "Survival after first progression by time to progression",
     "A shows Kaplan–Meier survival from first progression in eight cohorts by the interval from time origin to progression. "
     "B shows the adjusted hazard ratio from a natural cubic spline in the four primary cohorts, stratified by cohort, with "
     "the range observable at a month-12 landmark shaded. C shows cohort-specific hazard ratios per six months of later "
     "progression with the random-effects summary. D shows the gain in concordance in held-out cohorts when time to "
     "progression is added to age, stage and residual disease. Bars and bands are 95% intervals."),
    ("5", HERE / "figures/Figure_5_molecular.png", "Expression scores, hallmark programmes and BRCA1/2",
     "A shows pooled hazard ratios per standard deviation of three prespecified expression scores in eight cohorts. "
     "Post-landmark death is adjusted for baseline factors and progression status, and post-progression survival for "
     "baseline factors and time to progression. B shows pooled held-out changes in AUC when scores are added to the Status "
     "model. C compares hallmark enrichment for progression within 12 months between the primary and additional cohorts, "
     "with replicated hallmarks in red. D shows hazard ratios for BRCA1/2 mutation in sequenced TCGA cases and for "
     "progression status after adjustment for BRCA1/2. Bars are 95% intervals."),
    ("6", HERE / "figures/Figure_6_genomic_singlecell.png", "Genomic features and cellular origin of the early-progression programme",
     "A shows hazard ratios for death after the month-12 landmark in TCGA, adjusted for age, stage, residual disease and "
     "progression status, within 24 months and over all follow-up. B shows out-of-fold changes in AUC when each feature is "
     "added to the Status model. C shows, for each gene set, the mean specificity of its genes for six cell classes in "
     "single-cell data from treatment-naive tumours; numbers give the percentage of genes with their highest expression in "
     "that class. D shows hallmark enrichment for progression within 12 months in the primary cohorts before and after "
     "adjustment for the stromal score, for hallmarks replicated in either analysis. Bars are 95% intervals."),
]
SUPP_FIGURES = [
    ("S1", ROOT / "figures/Figure_S4_cohort_flow.png", "Cohort eligibility by step",
     "Counts show successive restrictions for the month-12 risk sets. TCGA starts with the 542 unredacted records and UKE "
     "with the 206 patients who had a primary sample and patient identifier."),
    ("S2", ROOT / "figures/Figure_1_survival.png", "Mortality and survival by progression history",
     "A shows an exploratory TCGA spline of predicted 12-month mortality by time since progression for age 60, stage III "
     "and optimal cytoreduction, with a 95% bootstrap band. B to E show landmark Kaplan–Meier survival by time since "
     "progression. Legends give the starting risk sets."),
    ("S3", ROOT / "figures/Figure_S3_sensitivity.png", "Prespecified sensitivity analyses of the TCGA-only equations",
     "Paired patient-bootstrap 95% intervals for AUC contrasts under alternative landmarks, horizon, predictor sets and "
     "treatment subsets. Japan contributes to the model without age."),
    ("S4", ROOT / "figures/Figure_6_subgroups.png", "Progression status and death in subgroups",
     "Cox hazard ratios for progression by month 12 through 24 months after the landmark, adjusted for the remaining "
     "baseline covariates. Counts give deaths and patients. Intervals are descriptive."),
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
    DOC.add_heading("Figures", 1)
    for i, item in enumerate(FIGURES):
        if i > 0:
            DOC.add_page_break()
        add_figure(*item)
    DOC.add_page_break()
    DOC.add_heading("Supplementary figures", 1)
    for i, item in enumerate(SUPP_FIGURES):
        if i in (1, 2, 3):
            DOC.add_page_break()
        add_figure(*item)


TRIPOD = [
    ("Title and abstract", "Title page; structured abstract"),
    ("Background and objectives", "Introduction"),
    ("Data sources, dates and settings", "Methods, Cohorts and participants; Table 1; source manifests"),
    ("Participants and eligibility", "Methods, Cohorts and participants; Supplementary Figure S1"),
    ("Outcome definition and assessment", "Methods, Predictors and outcome at the landmark; Table 1"),
    ("Predictors and their measurement", "Methods, Predictors and outcome at the landmark"),
    ("Sample size", "Methods, Development and validation strategies; Supplementary Table S9"),
    ("Missing data", "Methods, Development and validation strategies; Supplementary Table S9"),
    ("Analytical methods, model type and evaluation", "Methods, all subsections"),
    ("Class imbalance", "Not applicable; no resampling of outcome classes"),
    ("Fairness", "Not assessed; race and ethnicity unavailable (Discussion, limitations)"),
    ("Model output", "Probability of death within 12 months; Supplementary Table S6"),
    ("Development versus evaluation differences", "Table 1; Methods; Discussion, limitations"),
    ("Ethics", "Declarations"),
    ("Funding, competing interests", "Declarations (to be completed by authors)"),
    ("Protocol and registration", "Methods, Study design; not registered; Supplementary Table S2"),
    ("Data and code sharing", "Declarations"),
    ("Patient and public involvement", "Declarations; none"),
    ("Participant flow and characteristics", "Results; Table 1; Supplementary Figure S1; Supplementary Table S5"),
    ("Model development and specification", "Results; Supplementary Table S6"),
    ("Model performance with intervals", "Table 2; Figure 2; Supplementary Tables S7 and S8"),
    ("Model updating", "Results; Supplementary Table S8"),
    ("Interpretation, limitations, usability and future research", "Discussion"),
]
STATUS = [
    ("Month-12 landmark, 12-month horizon; Baseline, Status and Timing models fitted in TCGA and applied unchanged", "Prespecified"),
    ("Landmarks at 6 and 18 months, 24-month horizon, treatment and predictor-set sensitivity analyses", "Prespecified"),
    ("Cox associations and subgroup hazard ratios after the landmark", "Prespecified"),
    ("Progression-only benchmark; two-stage bootstrap; spline timing term; grade subsets", "Post hoc"),
    ("Internal–external cross-validation, random-effects summaries and prediction intervals", "Post hoc"),
    ("Recalibration, smoothed calibration, multiple imputation, precision targets", "Post hoc"),
    ("Landmarks at 9, 15, 21 and 24 months; landmark supermodel; illness–death model", "Post hoc"),
    ("Post-progression survival in eight cohorts", "Post hoc"),
    ("Expression scores, BRCA1/2 status, transcriptome-wide comparison and hallmark enrichment", "Post hoc, exploratory"),
    ("HRD score, CCNE1 amplification and other TCGA genomic features; single-cell attribution; stroma-adjusted enrichment", "Post hoc, exploratory"),
    ("TCGA slide stromal percentage; APOLLO-2 and additional-cohort application of the TCGA equation", "Post hoc, exploratory"),
    ("Eight machine-learning algorithms, penalisation diagnostics, nomogram, decision curves", "In the plan as secondary; removed"),
]
GENES = __import__("ast").literal_eval(
    (HERE / "v2_molecular.py").read_text(encoding="utf-8").split("SIGNATURES = ", 1)[1].split("\nLABEL", 1)[0])


def supplementary_tables():
    DOC.add_page_break()
    DOC.add_heading("Supplementary tables", 1)
    layout.table("Table S1 TRIPOD+AI reporting checklist, condensed", ["Item", "Where reported"],
                 [list(r) for r in TRIPOD], [2.9, 3.85],
                 "Condensed from the TRIPOD+AI checklist. The full journal form should be completed at submission.", 2)
    DOC.add_page_break()
    layout.table("Table S2 Status of each analysis", ["Analysis", "Status"], [list(r) for r in STATUS], [5.1, 1.65],
                 "The plan was written after inspection of source fields and before model fitting. It was not registered.", 2)
    names = {"stroma": "Stromal", "tcell": "T-cell", "prolif": "Proliferation"}
    layout.table("Table S3A Prespecified expression gene sets", ["Score", "Genes"],
                 [[names[k], ", ".join(v)] for k, v in GENES.items()], [1.0, 5.75],
                 "Scores are means of within-cohort gene z-scores. Gene coverage by cohort is in molecular_coverage.csv.", 2)
    hall = rows(NEW, "transcriptome_hallmarks")
    shown = [r for r in hall if float(r["FDR q-val_primary"]) < .05]
    layout.table("Table S3B Hallmark gene sets with false discovery rate below 0.05 in the primary cohorts",
                 ["Hallmark", "NES, primary", "FDR, primary", "NES, additional", "P, additional", "Replicated"],
                 [[r["Term"].replace("HALLMARK_", "").replace("_", " ").capitalize(), f(r["NES_primary"], 2),
                   pval(r["FDR q-val_primary"]), f(r["NES_additional"], 2), pval(r["NOM p-val_additional"]),
                   "Yes" if r["replicated"] == "True" else "No"] for r in shown],
                 [2.35, .9, .9, .95, .9, .75],
                 "NES, normalised enrichment score. Positive values indicate higher expression in tumours that progressed "
                 "within 12 months. All 50 hallmarks are in transcriptome_hallmarks.csv.")
    DOC.add_page_break()
    st = rows(PEER, "tcga_stroma_oof")[0]
    ap = {r["model"]: r for r in rows(PEER, "apollo_clinical_transport") if r["subset"] == "All sites"}
    ext = rows(REV, "additional_validation")
    japan = [r for r in rows(OLD, "model_metrics") if r["variant"] == "No age" and r["cohort"] == "Japan"]
    sets = {"APOLLO-2": (ap.values(), "Age, stage, residual"), "UVA": ([r for r in ext if r["cohort"] == "UVA"], "Age, stage, residual"),
            "Curie": ([r for r in ext if r["cohort"] == "Curie"], "Stage"),
            "OVCAD": ([r for r in ext if r["cohort"] == "OVCAD"], "Age, stage"), "Japan": (japan, "Stage, residual")}
    body = []
    for name, (data, feats) in sets.items():
        m = {r["model"]: r for r in data}
        body.append([name, feats, f'{m["Status"]["n"]} / {m["Status"]["deaths"]}',
                     *(f(m[k]["auc"]) if k in m else "—" for k in ("Baseline", "Progression only", "Status", "Timing")),
                     f(m["Status"]["oe"], 2)])
    layout.table("Table S4 Earlier supplementary analyses with the TCGA-only equations",
                 ["Cohort", "Baseline predictors", "N / deaths", "Baseline", "Progression only", "Status", "Timing", "Status O/E"],
                 body, [.8, 1.45, .8, .72, .8, .67, .67, .8],
                 "AUCs in additional cohorts using TCGA equations with the predictors each source provides. In 342 TCGA "
                 f"patients with a slide-level stromal-cell percentage from biospecimen quality control, adding it to the Status "
                 f"model changed the out-of-fold AUC by {ci(st['delta_auc'], st['delta_auc_lo'], st['delta_auc_hi'], 4)}. "
                 "These analyses are imprecise and are not used for inference.")
    s = {r["cohort"]: r for r in rows(OLD, "cohort_characteristics") if r["group"] == "All"}
    pc = lambda a, b: f"{a} ({100 * int(a) / int(b):.1f}%)"  # noqa: E731
    layout.table("Table S5 Characteristics of the month-12 risk sets",
                 ["Cohort", "N", "Age, median (IQR)", "Stage IV", "Residual >1 cm", "Bevacizumab", "Follow-up after landmark, months"],
                 [[k, s[k]["n"], f'{f(s[k]["age_median"], 1)} ({f(s[k]["age_q1"], 1)}–{f(s[k]["age_q3"], 1)})',
                   pc(s[k]["stage_iv"], s[k]["n"]), pc(s[k]["suboptimal"], s[k]["n"]),
                   pc(s[k]["bevacizumab"], s[k]["n"]) if k in ("UKE", "ICON7") else "—",
                   f(s[k]["followup_median"], 1)] for k in COHORTS],
                 [.7, .5, 1.3, 1.0, 1.1, 1.0, 1.15],
                 "Residual >1 cm includes the source's suboptimal category. Follow-up is the reverse Kaplan–Meier median.")
    DOC.add_page_break()
    co = {(r["fit"], r["model"]): r for r in rows(NEW, "iecv_coefficients")}
    terms = [("Intercept", "intercept"), ("Age per 10 years, centred at 60", "age10"), ("Stage IV", "iv"),
             ("Residual >1 cm", "rd"), ("Progression by month 12", "early"), ("Months since progression, (m − 6)/6", "timing6")]
    body = []
    for label, key in terms:
        row = [label]
        for model in ("Status", "Timing"):
            a, b = co["TCGA only", model], co["All four cohorts", model]
            lo, hi = co["All four cohorts, 2.5%", model], co["All four cohorts, 97.5%", model]
            row += [f(a[key]) if a[key] else "—", ci(b[key], lo[key], hi[key], 2) if b[key] else "—"]
        body.append(row)
    layout.table("Table S6 Model coefficients",
                 ["Term", "Status, TCGA only", "Status, four cohorts", "Timing, TCGA only", "Timing, four cohorts"],
                 body, [2.15, 1.0, 1.3, 1.0, 1.3],
                 "Log odds of death within 12 months of the landmark. Four-cohort intervals are 95% percentile intervals "
                 "from 500 bootstrap resamples within cohort. The four-cohort equations are given for transparency and "
                 "need local recalibration.")
    m = {(r["cohort"], r["model"]): r for r in rows(OLD, "model_metrics") if r["variant"] == "Primary"}
    po = {r["cohort"]: r for r in rows(OLD, "progression_only_benchmark")}
    cal = {(r["cohort"], r["estimate"]): r for r in rows(REV, "calibration_and_two_stage")}
    layout.table("Table S7 External performance of the prespecified TCGA-only equations",
                 ["Cohort", "Baseline AUC", "Progression-only AUC", "Status AUC", "Timing AUC", "Status Brier",
                  "Calibration intercept", "Calibration slope"],
                 [[k, f(m[k, "Baseline"]["auc"]), f(po[k]["auc"]), f(m[k, "Status"]["auc"]), f(m[k, "Timing"]["auc"]),
                   f(m[k, "Status"]["brier"]),
                   *(ci(cal[k, e]["value"], cal[k, e]["lo"], cal[k, e]["hi"], 2, "\n")
                     for e in ("calibration_intercept", "calibration_slope"))] for k in COHORTS[1:]],
                 [.7, .75, .95, .7, .7, .75, 1.1, 1.1],
                 "Intercept is estimated with the slope fixed at one. Parentheses contain 95% bootstrap intervals.")
    rec = rows(NEW, "recalibration")
    ici = {(r["model"], r["cohort"]): r for r in rows(NEW, "calibration_ici")}
    layout.table("Table S8 Updating the TCGA-only Status equation within each external cohort",
                 ["Cohort", "Method", "O/E", "Slope", "Brier", "AUC"],
                 [[r["cohort"], r["method"], f(r["oe"], 2), f(r["slope"], 2), f(r["brier"]), f(r["auc"])] for r in rec],
                 [.9, 2.2, .9, .9, .9, .9],
                 "Updated estimates are means over 20 repeats of five-fold cross-validation. Integrated calibration index "
                 "without updating: " + ", ".join(f'{k} {f(ici["TCGA only (locked)", k]["ici"])}' for k in COHORTS[1:])
                 + ". For held-out internal–external predictions: "
                 + ", ".join(f'{k} {f(ici["IECV", k]["ici"])}' for k in COHORTS) + ".", 2)
    DOC.add_page_break()
    imp = {r["cohort"]: r for r in rows(NEW, "imputation")}
    c = {(r["design"], r["cohort"]): r for r in rows(NEW, "iecv_cohort")}
    body = [[f"Multiple imputation, {k}", f'{int(float(imp[k]["n_mean"]))} / {int(float(imp[k]["deaths_mean"]))}',
             f(imp[k]["auc_status_mean"]), f(imp[k]["d_status_baseline_mean"]), f(imp[k]["d_timing_status_mean"]),
             f(imp[k]["oe_mean"], 2)] for k in COHORTS[1:]]
    for k in ("UKE", "ICON7"):
        r = c["IECV without German training cohorts", k]
        body.append([f"Training on TCGA and AOCS only, {k}", f'{r["n"]} / {r["deaths"]}', f(r["auc_Status"]),
                     f(r["d_status_baseline"]), f(r["d_timing_status"]), f(r["oe"], 2)])
    layout.table("Table S9A Sensitivity to missing residual disease and to possible cohort overlap",
                 ["Analysis and evaluation cohort", "N / deaths", "Status AUC", "Status minus Baseline", "Timing minus Status", "Status O/E"],
                 body, [2.5, .85, .8, .95, .9, .75],
                 "Imputation rows use TCGA-only equations refitted in each of 30 imputed data sets (386 TCGA patients, 39 "
                 "imputed) and show means across imputations.", 1)
    pr = rows(NEW, "precision")
    layout.table("Table S9B Events needed for precise external validation", ["Target or cohort", "Patients", "Deaths"],
                 [[r["target"], r["required_n"], r["required_events"]] for r in pr], [3.6, 1.5, 1.5],
                 f'Targets assume an outcome proportion of {f(pr[0]["outcome_proportion"], 3)} and the observed distribution of '
                 "the linear predictor. O/E and slope targets follow Riley and colleagues [28] and the C statistic target "
                 "uses the Hanley–McNeil variance [29].")
    DOC.add_page_break()
    dist = rows(NEW, "landmark_timing_distribution")
    hr = {(r["landmark"], r["term"]): r for r in rows(NEW, "supermodel_hr") if r["model"] == "Single landmark"}
    lm = {(int(r["landmark"]), r["contrast"]): r for r in rows(NEW, "landmark_iecv_pooled")}
    body = []
    for month in (6, 9, 12, 15, 18, 21, 24):
        d = [r for r in dist if int(r["landmark"]) == month]
        t = hr[str(month), "timing6"]
        body.append([month, sum(int(r["n"]) for r in d), sum(int(r["deaths12"]) for r in d),
                     sum(int(r["progressed"]) for r in d), sum(int(r["progressed_ge6"]) for r in d),
                     ci(t["hr"], t["lo"], t["hi"], 2, "\n") if month > 6 else "Not estimable",
                     *(ci(lm[month, k]["estimate"], lm[month, k]["lo"], lm[month, k]["hi"], sep="\n") for k in ("status", "timing"))])
    layout.table("Table S10 Landmark-specific risk sets and estimates in the four primary cohorts",
                 ["Landmark, months", "Patients", "Deaths in 12 months", "Progressed", "Progressed >6 months earlier",
                  "HR per 6 months since progression", "Status minus Baseline AUC", "Timing minus Status AUC"],
                 body, [.72, .68, .72, .78, .95, 1.0, .95, .95],
                 "Hazard ratios are from cohort-stratified Cox models at each landmark. AUC contrasts are random-effects "
                 "summaries of held-out estimates with 95% confidence intervals.")
    DOC.add_page_break()
    flow = {r["cohort"]: r for r in rows(NEW, "pps_flow")}
    ch = {r["cohort"]: r for r in rows(NEW, "pps_cohort_hr") if "within" not in r["adjusted_for"]}
    order = ["TCGA", "AOCS", "UKE", "ICON7", "Japan", "UVA", "Curie", "OVCAD"]
    layout.table("Table S11A Post-progression survival by cohort",
                 ["Cohort", "Progressed", "Excluded", "Analysed / deaths", "Months to progression, median (IQR)",
                  "Adjusted for", "HR per 6 months later progression"],
                 [[k, flow[k]["progressed"], flow[k]["death_or_last_contact_at_progression"],
                   f'{flow[k]["analysed"]} / {flow[k]["deaths"]}',
                   f'{f(flow[k]["ttp_median"], 1)} ({f(flow[k]["ttp_q1"], 1)}–{f(flow[k]["ttp_q3"], 1)})',
                   ch[k]["adjusted_for"].replace("age10", "age").replace("iv", "stage").replace("rd", "residual").replace("+", ", "),
                   ci(ch[k]["hr"], ch[k]["lo"], ch[k]["hi"], 2, "\n")] for k in order],
                 [.7, .8, .7, .95, 1.35, 1.2, 1.05],
                 "Excluded patients had a progression date equal to the date of death or last contact. Cohort models use "
                 "complete cases for the listed adjusters.")
    bands = rows(NEW, "pps_bands")
    layout.table("Table S11B Post-progression outcome by interval from time origin to progression",
                 ["Months to progression", "Patients / deaths, eight cohorts", "Median survival after progression, months",
                  "12-month mortality after progression", "Adjusted HR, primary cohorts"],
                 [[r["band"], f'{r["n_all"]} / {r["deaths_all"]}', f(r["median_pps"], 1),
                   f'{100 * float(r["mortality12"]):.1f}% ({100 * float(r["mortality12_lo"]):.1f} to {100 * float(r["mortality12_hi"]):.1f})',
                   ci(r["hr"], r["lo"], r["hi"], 2) if r["lo"] else "Reference"] for r in bands],
                 [1.1, 1.3, 1.45, 1.5, 1.4],
                 "Hazard ratios are adjusted for age, stage and residual disease and stratified by cohort.")
    DOC.add_page_break()
    mh = rows(NEW, "molecular_pooled_hr")
    short = {"Progression-free survival from origin": "Progression-free survival",
             "Overall survival from origin": "Overall survival",
             "Death within 24 months after the month-12 landmark, baseline adjusted": "Death ≤24 months after landmark, baseline adjusted",
             "Death within 24 months after the month-12 landmark, baseline and status adjusted": "Death ≤24 months after landmark, also status adjusted",
             "Post-progression survival, baseline and progression time adjusted": "Post-progression survival, also adjusted for time to progression"}
    layout.table("Table S12 Expression scores and outcome, random-effects summaries of eight cohorts",
                 ["Score", "Outcome", "Patients / events", "HR per SD (95% CI)", "I²"],
                 [[names[r["score"]], short[r["analysis"]], f'{r["n"]} / {r["events"]}', ci(r["hr"], r["lo"], r["hi"], 2),
                   f'{100 * float(r["i2"]):.0f}%'] for r in mh],
                 [.95, 3.0, 1.05, 1.25, .5],
                 "Cohort models were adjusted for the available baseline factors among age, stage and residual disease.", 2)
    br = rows(NEW, "brca")
    layout.table("Table S13 BRCA1/2 mutation in sequenced TCGA cases",
                 ["Analysis", "Term", "Patients / events", "Carriers", "Estimate (95% CI)"],
                 [[r["analysis"], {"brca": "BRCA1/2 mutation", "early": "Progression by month 12", "ttp6": "Later progression, per 6 months",
                                   "brca_early": "Mutation × progression"}[r["term"]], f'{r["n"]} / {r["events"]}',
                   r["carriers"], ci(r["hr"], r["lo"], r["hi"], 3 if "AUC" in r["analysis"] else 2)] for r in br],
                 [2.55, 1.55, .9, .6, 1.15],
                 "Hazard ratios are adjusted for age, stage and residual disease. The last row is a change in AUC. "
                 "Mutations are the germline and somatic calls of the TCGA marker study [32].", 2)

    DOC.add_page_break()
    gen = rows(NEW, "genomics")
    short = {"Progression-free survival from diagnosis": "Progression-free survival",
             "Overall survival from diagnosis": "Overall survival",
             "Death within 24 months after the month-12 landmark, status adjusted": "Death ≤24 months after landmark, status adjusted",
             "Death after the month-12 landmark, all follow-up, status adjusted": "Death after landmark, all follow-up, status adjusted",
             "Post-progression survival, progression time adjusted": "Post-progression survival, adjusted for time to progression",
             "Out-of-fold change in AUC when added to Status": "Change in AUC when added to Status"}
    term = {"early": "Progression by month 12", "brca": "BRCA1/2 mutation", "hrd10": "HRD score, per 10 points",
            "ccne1": "CCNE1 amplification", "ttp6": "Later progression, per 6 months"}
    body = []
    for r in gen:
        joint = r["feature"].startswith("Joint")
        label = "Joint model" if joint else r["feature"]
        what = short[r["analysis"]] + (f': {term.get(r["term"], "all three features")}' if joint else "")
        body.append([label, what, f'{r["n"]} / {r["events"]}',
                     ci(r["estimate"], r["lo"], r["hi"], 3 if "AUC" in r["analysis"] else 2)])
    layout.table("Table S14 HRD score, CCNE1 amplification and other genomic features in TCGA",
                 ["Feature", "Analysis", "Patients / events", "Estimate (95% CI)"], body, [1.75, 2.9, .95, 1.15],
                 "Estimates are hazard ratios adjusted for age, stage and residual disease, except rows labelled change in "
                 "AUC. The joint model contains BRCA1/2 mutation, HRD score and CCNE1 amplification. HRD scores are from "
                 "the TCGA pan-cancer DNA damage repair analysis [34] and the other features from the PanCancer Atlas "
                 "release in cBioPortal [33].", 2)
    DOC.add_page_break()
    att = rows(NEW, "sc_attribution")
    order = ["Fibroblast", "Mesothelial", "Malignant", "Myeloid", "T/NK", "B/plasma"]
    names_sc = list(dict.fromkeys(r["gene_set"] for r in att))
    cell = {(r["gene_set"], r["cell_class"]): r for r in att}
    body = [[n, cell[n, order[0]]["genes"],
             *(f'{100 * float(cell[n, c]["share_highest"]):.0f}% ({f(cell[n, c]["mean_specificity"], 2)})' for c in order)]
            for n in names_sc]
    body.append(["All genes in the comparison", "5678",
                 *(f'{100 * float(cell[names_sc[0], c]["share_highest_all_genes"]):.0f}%' for c in order)])
    layout.table("Table S15 Cell classes expressing each gene set in single-cell data",
                 ["Gene set", "Genes", *order], body, [1.75, .5, .78, .78, .75, .7, .7, .75],
                 "Cells are the percentage of genes with their highest expression in the class and, in parentheses, the mean "
                 "specificity z-score. Treatment-naive cells from GSE165897 [38]. Permutation P values are in sc_attribution.csv.")
    adj = rows(NEW, "transcriptome_hallmarks_stroma_adjusted")
    keep = [r for r in adj if r["replicated"] == "True" or r["adj_replicated"] == "True"]
    layout.table("Table S16 Hallmark enrichment before and after adjustment for the stromal score",
                 ["Hallmark", "NES, primary", "NES, additional", "Adjusted NES, primary", "Adjusted FDR, primary",
                  "Adjusted NES, additional", "Replicated after adjustment"],
                 [[r["Term"].replace("HALLMARK_", "").replace("_", " ").capitalize(), f(r["NES_primary"], 2),
                   f(r["NES_additional"], 2), f(r["adj_NES_primary"], 2), pval(r["adj_FDR q-val_primary"]),
                   f(r["adj_NES_additional"], 2), "Yes" if r["adj_replicated"] == "True" else "No"] for r in keep],
                 [1.9, .75, .8, .85, .85, .85, .75],
                 "Hallmarks replicated in the unadjusted or the adjusted comparison. Negative values indicate lower "
                 "expression in tumours that progressed within 12 months. Adjustment used the stromal score as a covariate "
                 "in gene-wise linear models within cohort.")


if __name__ == "__main__":
    content()
    main_tables()
    figures()
    supplementary_tables()
    DOC.save(HERE / "Ovarian_Cancer_Manuscript_v2.docx")
    print("saved")

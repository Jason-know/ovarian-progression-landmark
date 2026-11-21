"""Build the English manuscript from the text, figure and table inputs."""
from pathlib import Path
import csv
import json
import re

from docx import Document
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = ROOT / "figures"
DOC = Document()


def style(name, size, before=0, after=0, bold=False, color="171A21"):
    styles = DOC.styles
    s = styles[name] if name in styles else styles.add_style(name, WD_STYLE_TYPE.PARAGRAPH)
    s.font.name = "Calibri"
    s.font.size = Pt(size)
    s.font.bold = bold
    s.font.color.rgb = RGBColor.from_string(color)
    pf = s.paragraph_format
    pf.space_before, pf.space_after = Pt(before), Pt(after)
    pf.line_spacing = 1.13
    node = s.element.get_or_add_pPr()
    for b in node.findall(qn("w:pBdr")):
        node.remove(b)
    return s


def setup():
    section = DOC.sections[0]
    section.page_width, section.page_height = Inches(8.27), Inches(11.69)
    section.top_margin, section.bottom_margin = Inches(.73), Inches(.73)
    section.left_margin = section.right_margin = Inches(.75)
    section.footer_distance = Inches(.35)
    style("Normal", 11, after=7)
    style("Title", 25, after=9, bold=True, color="000000")
    style("Subtitle", 13, after=16, color="424A55")
    style("Heading 1", 15, before=12, after=6, bold=True)
    style("Heading 2", 11.5, before=8, after=4, bold=True)
    style("Figure Caption", 9.5, before=2, after=9)
    style("Source Reference", 9, after=3)
    style("Table Note", 9.3, before=1, after=7)
    style("Table Text", 9.3, after=1)
    style("Footer", 8.5, color="657186")
    DOC.styles["Subtitle"].paragraph_format.keep_with_next = True
    DOC.styles["Title"].paragraph_format.keep_with_next = True
    DOC.styles["Source Reference"].paragraph_format.left_indent = Inches(.22)
    DOC.styles["Source Reference"].paragraph_format.first_line_indent = Inches(-.22)
    DOC.styles["Table Text"].paragraph_format.line_spacing = 1.06
    DOC.core_properties.title = "Progression Status and Timing for Updating Mortality Risk in Advanced Serous Ovarian Cancer"
    DOC.core_properties.subject = "Landmark prediction with independent external validation"
    DOC.core_properties.author = ""
    DOC.core_properties.keywords = "serous ovarian cancer, progression, landmark validation"
    footer = section.footer.paragraphs[0]
    footer.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    field = OxmlElement("w:fldSimple")
    field.set(qn("w:instr"), "PAGE")
    footer._p.append(field)


def paragraph(text, style_name="Normal"):
    p = DOC.add_paragraph(style=style_name)
    for n, val in enumerate(re.split(r"\*\*(.*?)\*\*", text)):
        if val:
            p.add_run(val).bold = bool(n % 2)
    return p


def hyperlink(p, address, label):
    url_kind = "http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink"
    link = OxmlElement("w:hyperlink")
    link.set(qn("r:id"), p.part.relate_to(address, url_kind, is_external=True))
    run = OxmlElement("w:r")
    props = OxmlElement("w:rPr")
    color = OxmlElement("w:color")
    color.set(qn("w:val"), "2A5885")
    props.append(color)
    run.append(props)
    t = OxmlElement("w:t")
    t.text = label
    run.append(t)
    link.append(run)
    p._p.append(link)


def main_text():
    in_refs = False
    for line in (ROOT / "scripts" / "manuscript.txt").read_text(encoding="utf-8").splitlines():
        kind, text = line.split("|", 1)
        if kind == "REF":
            if not in_refs:
                DOC.add_heading("References", 1)
                in_refs = True
            body, url = text.rsplit(" https://", 1)
            p = paragraph(body + " ", "Source Reference")
            label = "doi:" + url.removeprefix("doi.org/") if url.startswith("doi.org/") else "Source"
            hyperlink(p, "https://" + url, label)
        elif kind in ("H1", "H2"):
            DOC.add_heading(text, int(kind[-1]))
        elif kind == "TITLE":
            DOC.add_paragraph(text, style="Title")
        elif kind == "SUBTITLE":
            DOC.add_paragraph(text, style="Subtitle")
        else:
            paragraph(text)


def read(name):
    with (RES / f"{name}.csv").open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def number(s, d=3):
    return f"{float(s):.{d}f}" if s else "—"


def pct(n, d):
    return f"{100 * int(n) / int(d):.1f}%"


def cell_fill(cell, color):
    element = OxmlElement("w:shd")
    element.set(qn("w:fill"), color)
    cell._tc.get_or_add_tcPr().append(element)


def table(title, headers, rows, widths, note, left_columns=1):
    DOC.add_heading(title, 2)
    t = DOC.add_table(rows=1, cols=len(headers))
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = False
    for column, width in zip(t.columns, widths, strict=True):
        column.width = Inches(width)
    for i, values in enumerate([headers, *rows]):
        row = t.rows[0] if i == 0 else t.add_row()
        trpr = row._tr.get_or_add_trPr()
        trpr.append(OxmlElement("w:cantSplit"))
        if i == 0:
            trpr.append(OxmlElement("w:tblHeader"))
        for j, (value, width) in enumerate(zip(values, widths, strict=True)):
            cell = row.cells[j]
            cell.width = Inches(width)
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            p = cell.paragraphs[0]
            p.style = "Table Text"
            p.paragraph_format.space_before = Pt(4)
            p.paragraph_format.space_after = Pt(4)
            if j >= left_columns:
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            run = p.add_run(str(value))
            if i == 0:
                run.bold = True
                cell_fill(cell, "E7EDF3")
            elif i % 2 == 0:
                cell_fill(cell, "F4F6F8")
    paragraph(note, "Table Note")


def main_tables():
    DOC.add_page_break()
    DOC.add_heading("Main tables", 1)
    summary = [r for r in read("cohort_characteristics") if r["group"] == "All"]
    fields = [
        ("Participants", lambda r: r["n"]),
        ("Age in years, median (IQR)",
         lambda r: f'{number(r["age_median"], 1)} ({number(r["age_q1"], 1)}–{number(r["age_q3"], 1)})'),
        ("Stage IV", lambda r: f'{r["stage_iv"]} ({pct(r["stage_iv"], r["n"])})'),
        ("Suboptimal cytoreduction", lambda r: f'{r["suboptimal"]} ({pct(r["suboptimal"], r["n"])})'),
        ("Progression by month 12", lambda r: f'{r["early"]} ({pct(r["early"], r["n"])})'),
        ("Deaths in next 12 months", lambda r: r["dead12"]),
        ("Observed beyond next 12 months", lambda r: r["known_alive12"]),
        ("Estimated mortality at horizon",
         lambda r: f'{100*float(r["risk12"]):.1f}%'),
        ("Median follow-up after landmark, months", lambda r: number(r["followup_median"], 1)),
    ]
    table("Table 1 Landmark cohort characteristics", ["Characteristic", *(r["cohort"] for r in summary)],
          [[label, *(fn(r) for r in summary)] for label, fn in fields], [2.47, 1.07, 1.07, 1.07, 1.07],
          "Mortality uses Kaplan–Meier estimates. Median follow-up uses reverse Kaplan–Meier.")
    metric = { (r["cohort"], r["model"]): r for r in read("model_metrics") if r["variant"] == "Primary" }
    delta = { (r["cohort"], r["contrast"]): r for r in read("model_contrasts") if r["variant"] == "Primary" }
    rows = []
    for r in summary:
        c = r["cohort"]
        vals = [metric[c, name] for name in ("Baseline", "Status", "Timing")]
        s, t = (delta[c, name] for name in ("Status vs baseline", "Timing vs status"))
        change = [f'{number(d["delta"])}\n({number(d["lo"])} to {number(d["hi"])})' if d["lo"]
                  else number(d["delta"]) for d in (s, t)]
        rows.append([c, *(number(m["auc"]) for m in vals), *change, number(vals[1]["oe"], 2)])
    table("Table 2 Discrimination and calibration at the primary horizon",
          ["Cohort", "Baseline\nAUC", "Status\nAUC", "Timing\nAUC", "Status minus baseline\nAUC change",
           "Timing minus status\nAUC change", "Status\nO/E"],
          rows, [.75, .75, .75, .75, 1.47, 1.48, .8],
          "The two AUC contrasts use paired estimates on the same patients. Parentheses contain 95% external "
          "bootstrap CIs. TCGA estimates are out of fold. O/E is observed-to-expected mortality.")


FIGURES = [
    ("1", "Figure_1_survival", "Mortality and survival by progression history",
     "Panel A shows an exploratory TCGA spline on the predicted 12-month mortality scale for a fixed age 60, stage III and optimal cytoreduction profile. The band shows 95% bootstrap limits. The rug shows time since progression. Panels B–E show unadjusted landmark Kaplan–Meier survival by time since progression. Marks denote censoring. Numbers in the legends are the starting risk sets."),
    ("2", "Figure_2_landmark_evidence", "Landmark design and observed mortality",
     "Panel A depicts three illustrative pre-landmark histories and the later death ascertainment window. Panel B shows cohort-specific Kaplan–Meier 12-month mortality with pointwise 95% intervals, classified by documented progression at month 12."),
    ("3", "Figure_3_predictor_diagnostics", "Predictor diagnostics in the development cohort",
     "Panel A shows LASSO coefficient paths on standardized predictors. Panel B shows five-fold mean Brier scores and standard errors. The vertical line marks the minimum. Panel C compares observed predictor importance with the maximum shuffled feature importance. Panel D compares nonzero LASSO terms with predictors exceeding that maximum in at least 90% of permutations."),
    ("4", "Figure_4_roc", "Time-dependent ROC curves across eight algorithms",
     "The models used the same five Timing predictors. TCGA predictions were out of fold. External curves used the full TCGA fits without any update. Legends report cumulative/dynamic AUC with censoring weights."),
    ("5", "Figure_5_calibration_and_decision_curves", "Calibration and decision curves",
     "Panels A–D show predicted mean risk and Kaplan–Meier mortality by prediction quartile. Panels E–H show IPCW net benefit for hypothetical risk flagging. The horizontal dashed line represents flagging nobody. The dotted curve represents flagging everybody."),
    ("6", "Figure_6_subgroups", "Adjusted associations with early progression",
     "Cox hazard ratios compare progression by month 12 with no documented progression, using follow-up through 24 months after the landmark. Subgroup counts give deaths and participants. Models included the remaining baseline covariates with variation. Intervals are descriptive. Arrows indicate intervals extending beyond the plot."),
    ("S1", "Figure_S1_nomogram", "Nomogram for the TCGA Status model",
     "Map each predictor to the points axis and sum the four values. Total points map to 12-month mortality. The stage coefficient is negative in this selected landmark cohort. The nomogram displays the uncalibrated locked model for reproducibility."),
    ("S2", "Figure_S2_external_timing", "Exploratory timing curves in external cohorts",
     "Separate weighted spline fits in each external cohort for a fixed age 60, stage III and optimal cytoreduction profile. Shading shows 95% bootstrap limits. Ranges cover the central 90% of progression times in each cohort. Dashed lines show fitted risks for a profile without progression. These fits illustrate associations within cohorts."),
    ("S3", "Figure_S3_sensitivity", "Sensitivity analyses for AUC contrasts",
     "The horizontal intervals are paired patient-bootstrap 95% intervals. Japan contributes to the model without age. The primary analysis used 1,000 resamples. Other analyses used 350 resamples."),
    ("S4", "Figure_S4_cohort_flow", "Cohort eligibility by step",
     "Counts show successive restrictions. TCGA starts with the 542 unredacted OV rows. UKE starts with the 206 unique patients having an available primary sample and patient identifier."),
]


def figure(number_, name, title, caption):
    DOC.add_heading(f"Figure {number_} {title}", 2)
    p = DOC.add_paragraph()
    p.paragraph_format.keep_with_next = True
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    picture = p.add_run().add_picture(str(FIG / f"{name}.png"), width=Inches(6.7))
    picture._inline.docPr.set("descr", title + ". " + caption)
    paragraph(caption, "Figure Caption")


def figures():
    # Give the denser four-panel figures a full page and pair shorter figures.
    for i, row in enumerate(FIGURES):
        if i in (0, 2, 3, 4, 5, 6, 8):
            DOC.add_page_break()
        if i == 0:
            DOC.add_heading("Main figures", 1)
        if i == 6:
            DOC.add_heading("Supplementary figures", 1)
        figure(*row)


def supplement_tables():
    DOC.add_page_break()
    DOC.add_heading("Supplementary tables", 1)
    sources = [
        ["TCGA", "CDR OS and PFI events and days, clinical stage and age.",
         "Broad Firehose patient.tumor_residual_disease. 0–10 mm or >10 mm. Diagnosis."],
        ["AOCS", "Tothill clinical RData. MosToDeath, MosToRelapse, Status, Age and Stage.",
         "Nil or <1 cm versus >1 cm. Diagnosis. Recorded platinum and no neoadjuvant chemotherapy."],
        ["UKE", "Table S1b. OS_time, OS_event, PFS_time, PFS_event, age and FIGO_stage.",
         "No tumor or <1 cm versus >1 cm. Diagnosis. Platinum with or without bevacizumab."],
        ["ICON7", "GSE140082 final_ostm, final_osid, final_pfstm, final_pfsid, age and figo_stage.",
         "Source optimal versus suboptimal. Trial randomization. Both assigned treatment arms."],
        ["Japan", "GSE32062 GPL6480. OS and PFS in months, death and recurrence flags, stage.",
         "Source optimal versus suboptimal. Published PFS and OS clocks. No individual age field."],
    ]
    table("Table S1 Data source and field mapping", ["Cohort", "Predictor and outcome fields", "Residual or origin details"],
          sources, [.8, 2.98, 2.97], "Cohort-specific protocols determined the definitions and ascertainment of progression.",
          left_columns=3)
    metrics = [r for r in read("model_metrics") if r["variant"] == "Primary"]
    models = list(dict.fromkeys(r["model"] for r in metrics if r["model"] not in ("Baseline", "Status", "Timing")))
    status = {r["cohort"]: r for r in read("progression_only_benchmark")}
    ids = ["TCGA", "AOCS", "UKE", "ICON7"]
    rows = [[mod, *(number(next(r["auc"] for r in metrics if r["model"] == mod and r["cohort"] == code))
                    for code in ids)] for mod in models]
    rows += [["Progression only (post hoc)", *(number(status[c]["auc"]) for c in ids)]]
    table("Table S2 AUC for eight algorithms and the progression-only benchmark", ["Model", *ids], rows,
          [2.45, 1.075, 1.075, 1.075, 1.075], "All eight algorithms outside the benchmark used the Timing predictors.")
    DOC.add_page_break()
    eqs = json.loads((RES / "equations.json").read_text(encoding="utf-8"))
    row_names = [("Intercept", "intercept"), ("Age per 10 years centered at 60", "age10"), ("Stage IV", "iv"),
                 ("Suboptimal cytoreduction", "rd"), ("Progression by month 12", "early"),
                 ("Timing term among progressed", "timing6")]
    table("Table S3 Locked equation coefficients", ["Term", *eqs],
          [[label, *(number(str(terms[key]), 5) if key in terms else "—" for terms in eqs.values())]
           for label, key in row_names],
          [3.24, 1.17, 1.17, 1.17],
          "Death risk equals 1 divided by 1 plus exp of the negative linear predictor. The Timing term is the "
          "elapsed interval since progression, centered at six months and divided by six, among patients "
          "with documented progression. Its value is zero in others.")
    p = [r for r in metrics if r["model"] in ("Baseline", "Status", "Timing")]
    table("Table S4 Absolute prediction performance", ["Cohort", "Model", "Brier", "Observed", "Expected", "O/E", "Intercept", "Slope"],
          [[r["cohort"], r["model"], *(number(r[k]) for k in
            ["brier", "observed", "expected", "oe", "intercept", "slope"])] for r in p],
          [.76, .88, .67, .85, .87, .76, .98, .98],
          "Observed risk uses censoring weights. Expected risk is mean predicted probability. Calibration "
          "intercept holds slope at one. The accompanying results provide AUC CIs and sensitivity results.")


if __name__ == "__main__":
    setup()
    main_text()
    main_tables()
    figures()
    supplement_tables()
    DOC.save(ROOT / "Ovarian_Cancer_Landmark_Manuscript.docx")

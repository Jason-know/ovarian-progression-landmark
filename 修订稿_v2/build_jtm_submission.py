"""Build the Journal of Translational Medicine submission files from the second-revision sources.

python build_jtm_submission.py main        -> JTM投稿包/Manuscript.docx (double spaced, line numbers, tables at end)
python build_jtm_submission.py additional  -> JTM投稿包/Additional_file_1.docx (supplementary figures and tables)
python build_jtm_submission.py letter      -> JTM投稿包/Cover_letter.docx and figure files
"""
from pathlib import Path
import re
import shutil
import sys

from docx import Document
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PACK = ROOT / "JTM投稿包"
MODE = sys.argv[1] if len(sys.argv) > 1 else "main"
REPO = "[repository URL to be inserted after the code is deposited]"
ARCHIVE = "[Zenodo DOI to be inserted after the release is archived]"

ABSTRACT = [
    "**Background.** After first-line treatment of advanced serous tubo-ovarian carcinoma, progression status at a follow-up "
    "visit is known to predict survival. Whether the timing of progression adds prognostic information, and why its apparent "
    "value differs between analyses, is unclear.",
    "**Methods.** We analysed public individual-patient data from eight cohorts with stage III–IV serous carcinoma. In four "
    "cohorts with age, stage and residual disease, censoring-weighted logistic models predicted death within 12 months of a "
    "month-12 landmark, using prespecified equations fitted in TCGA and post hoc internal–external cross-validation. Further "
    "post hoc analyses used landmarks from 6 to 24 months, an illness–death model, post-progression survival, tumour expression "
    "scores, hallmark enrichment with single-cell attribution, and BRCA1/2, homologous recombination deficiency (HRD) and CCNE1 "
    "status in TCGA.",
    "**Results.** At month 12 (850 patients, 142 deaths), progression status raised the held-out AUC over baseline factors by "
    "0.180 (95% CI 0.095 to 0.265), whereas time since progression changed it by −0.000 (−0.004 to 0.004). Only 19 of 221 "
    "progressed patients had progressed more than six months earlier. The increment from timing rose to 0.071 (0.006 to 0.136) "
    "at month 24. Among 1126 patients followed from first progression, each six months of later progression was associated "
    "with a pooled hazard ratio for death of 0.73 (0.66 to 0.82). Observed-to-expected mortality in held-out cohorts ranged "
    "from 0.67 to 1.85. A T-cell expression score remained associated with death after adjustment for progression status "
    "(hazard ratio per SD 0.86, 0.75 to 0.98) but added 0.007 to the AUC. BRCA1/2 mutation, HRD score and CCNE1 amplification "
    "were associated with long-term but not 24-month mortality. Epithelial–mesenchymal transition genes enriched in "
    "early-progressing tumours were expressed mainly by fibroblasts, and the enrichment disappeared after adjustment for a "
    "stromal score.",
    "**Conclusions.** Progression timing is strongly prognostic, but a month-12 landmark observes too narrow a range of it to "
    "improve prediction beyond progression status. Predicted probabilities required local recalibration. Baseline expression "
    "and genomic features added little to short-term prediction once progression status was known.",
]
ABBREVIATIONS = ("AOCS: Australian Ovarian Cancer Study; AUC: area under the receiver operating characteristic curve; CI: "
                 "confidence interval; FDR: false discovery rate; HR: hazard ratio; HRD: homologous recombination deficiency; "
                 "NES: normalised enrichment score; O/E: observed-to-expected ratio; OS: overall survival; PFS: progression-free "
                 "survival; SD: standard deviation; TCGA: The Cancer Genome Atlas; UKE: University Medical Center "
                 "Hamburg-Eppendorf.")
DECLARATIONS = [
    ("H2", "Ethics approval and consent to participate"),
    ("P", "Not applicable. All analyses used existing de-identified public data. The original studies report their ethics "
          "approvals and consent procedures. No new human participants were involved."),
    ("H2", "Consent for publication"), ("P", "Not applicable."),
    ("H2", "Availability of data and materials"),
    ("P", "The datasets analysed during the current study are publicly available: the TCGA Pan-Cancer Clinical Data Resource "
          "and Broad Firehose clinical data [9, 10], TCGA HRD scores [34] and TCGA PanCancer Atlas data in cBioPortal [33]; "
          "Gene Expression Omnibus accessions GSE9891, GSE140082, GSE32062, GSE30161, GSE26193, GSE49997 and GSE165897; the "
          "curatedOvarianData collection [15]; and the UKE patient-level tables in the supplementary material of reference "
          "[12]. File-level URLs and checksums are listed in the source manifests in the code repository. Project name: "
          "Progression timing in serous ovarian carcinoma. Project home page: " + REPO + ". Archived version: " + ARCHIVE +
          ". Operating system: platform independent. Programming language: Python 3.12. Other requirements: numpy, pandas, "
          "scipy, scikit-learn, lifelines, patsy, gseapy, matplotlib, python-docx. License: MIT. Any restrictions to use by "
          "non-academics: none."),
    ("H2", "Competing interests"), ("P", "The authors declare that they have no competing interests."),
    ("H2", "Funding"),
    ("P", "This work was supported by the Natural Science Foundation of Guangdong Province (2026A1515012040) and the Project "
          "of Administration of Traditional Chinese Medicine of Guangdong Province of China (20261063) to Ye Chen; the Youth "
          "Program of Medical Science Research at Wuhan Eighth Hospital (YX24Q02) to Huang Chen; and the General Project of "
          "Medical Science Research at Wuhan Eighth Hospital (YX24Y02) to Lijuan Zhao. The funders had no role in study "
          "design, analysis, interpretation or writing."),
    ("H2", "Authors' contributions"),
    ("P", "YC conceived and supervised the study. HC performed data curation, formal analysis and figure preparation and "
          "wrote the original manuscript draft. LZ and XL contributed to data interpretation and literature review. All "
          "authors critically revised the manuscript, read and approved the final manuscript and agree to be accountable "
          "for the work."),
    ("H2", "Acknowledgements"),
    ("P", "The results are in part based on data generated by the TCGA Research Network and by the investigators of the "
          "AOCS, UKE, ICON7/AGO-OVAR11, Japanese, UVA, Curie and OVCAD studies and of the GSE165897 single-cell study."),
]
AI_USE = ("An AI coding assistant (Devin, Cognition) was used to write analysis code and to draft and edit the manuscript "
          "text. The authors reviewed and verified the code, results, text and references and take full responsibility for "
          "the content.")


def journal_text(text):
    text = re.sub(r"Supplementary Tables S(\d+) and S(\d+)", r"Additional file 1: Tables S\1 and S\2", text)
    text = text.replace("Supplementary Table S", "Additional file 1: Table S")
    text = text.replace("Supplementary Figures S", "Additional file 1: Figures S")
    return text.replace("Supplementary Figure S", "Additional file 1: Figure S")


def source_lines():
    lines = (HERE / "v2_manuscript.txt").read_text(encoding="utf-8").splitlines()
    out, skip = [], False
    for line in lines:
        kind, text = line.split("|", 1)
        if text.startswith(("**Word count.**", "**Use of artificial intelligence.**")):
            continue
        if line == "H1|Abstract":
            out += [line, *("P|" + p for p in ABSTRACT)]
            skip = True
            continue
        if skip and text.startswith("**Keywords.**"):
            skip = False
        if skip:
            continue
        if line == "H1|Introduction":
            line = "H1|Background"
        if text.startswith("Two earlier supplementary analyses are retained"):
            line = line + " " + AI_USE
        if line == "H1|Declarations":
            out += ["H1|List of abbreviations", "P|" + ABBREVIATIONS, line, *(f"{k}|{t}" for k, t in DECLARATIONS)]
            skip = True
            continue
        if kind == "REF":
            skip = False
        if skip:
            continue
        out.append(line)
    return out


def main_document():
    import build_v2_manuscript as b
    layout, doc = b.layout, b.DOC
    layout.cell_fill = lambda *a, **k: None
    type(doc).add_page_break = lambda self: None
    paragraph, table = layout.paragraph, layout.table
    layout.paragraph = lambda text, style="Normal": paragraph(journal_text(text), style)

    def plain_table(title, headers, rows, widths, note, left_columns=1):
        table(title, headers, [[journal_text(str(c)) for c in r] for r in rows], widths, journal_text(note), left_columns)
        doc.tables[-1].style = "Table Grid"
    layout.table = plain_table
    layout.setup()
    doc.core_properties.title = "Progression status and progression timing in advanced serous tubo-ovarian carcinoma"
    in_refs = False
    for line in source_lines():
        kind, text = line.split("|", 1)
        if kind == "REF":
            if not in_refs:
                doc.add_heading("References", 1)
                in_refs = True
            body, url = text.rsplit(" https://", 1)
            if url.startswith("doi.org/"):
                layout.paragraph(f"{body} https://{url}", "Source Reference")
            else:
                layout.paragraph(f"{body} https://{url}. Accessed 20 Nov 2026.", "Source Reference")
        elif kind.startswith("H"):
            doc.add_heading(text, int(kind[1:]))
        elif kind in ("TITLE", "SUBTITLE"):
            doc.add_paragraph(text, style=kind.title())
        else:
            layout.paragraph(text)
    doc.add_heading("Figure legends", 1)
    for number, _, title, caption in b.FIGURES:
        layout.paragraph(f"**Fig. {number} {title}.** {caption}")
    doc.add_heading("Additional files", 1)
    layout.paragraph("**Additional file 1 (.docx).** Supplementary figures and tables. Figure S1 cohort eligibility. Figure S2 "
                     "survival by progression history. Figure S3 prespecified sensitivity analyses. Figure S4 subgroup "
                     "hazard ratios. Table S1 TRIPOD+AI checklist. Table S2 status of each analysis. Table S3 gene sets and "
                     "hallmark enrichment. Table S4 earlier supplementary analyses. Table S5 characteristics of the risk "
                     "sets. Table S6 model coefficients. Table S7 performance of the TCGA-only equations. Table S8 "
                     "recalibration. Table S9 imputation, cohort-overlap sensitivity and precision targets. Table S10 "
                     "landmark-specific estimates. Table S11 post-progression survival. Table S12 expression scores. Table "
                     "S13 BRCA1/2. Table S14 HRD score, CCNE1 and other genomic features. Table S15 single-cell "
                     "attribution. Table S16 stroma-adjusted hallmark enrichment.")
    b.main_tables()
    for style in ("Normal", "Source Reference"):
        doc.styles[style].paragraph_format.line_spacing = 2
    numbering = OxmlElement("w:lnNumType")
    numbering.set(qn("w:countBy"), "1")
    numbering.set(qn("w:restart"), "continuous")
    doc.sections[0]._sectPr.append(numbering)
    doc.save(PACK / "Manuscript.docx")


def additional_document():
    import build_v2_manuscript as b
    layout, doc = b.layout, b.DOC
    layout.cell_fill = lambda *a, **k: None
    table = layout.table

    def plain_table(*args, **kwargs):
        table(*args, **kwargs)
        doc.tables[-1].style = "Table Grid"
    layout.table = plain_table
    layout.setup()
    doc.add_paragraph("Additional file 1", style="Title")
    doc.add_paragraph("Progression status and progression timing as updated prognostic information in advanced serous "
                      "tubo-ovarian carcinoma. Supplementary figures S1 to S4 and tables S1 to S16.", style="Subtitle")
    doc.add_heading("Supplementary figures", 1)
    for i, item in enumerate(b.SUPP_FIGURES):
        if i:
            doc.add_page_break()
        b.add_figure(*item)
    b.supplementary_tables()
    doc.save(PACK / "Additional_file_1.docx")


def letter_and_figures():
    sources = [ROOT / "figures/Figure_2_landmark_evidence", HERE / "figures/Figure_2_performance",
               HERE / "figures/Figure_3_landmarks", HERE / "figures/Figure_4_post_progression",
               HERE / "figures/Figure_5_molecular", HERE / "figures/Figure_6_genomic_singlecell"]
    for n, src in enumerate(sources, 1):
        for ext in (".pdf", ".png"):
            shutil.copyfile(str(src) + ext, PACK / f"Figure_{n}{ext}")
    doc = Document()
    style = doc.styles["Normal"]
    style.font.name, style.font.size = "Calibri", Pt(11)
    for text in [
        "20 November 2026",
        "Francesco Marincola, MD\nEditor-in-Chief, Journal of Translational Medicine",
        "Dear Dr Marincola,",
        "We submit the manuscript “Progression Status and Progression Timing as Updated Prognostic Information in Advanced "
        "Serous Tubo-Ovarian Carcinoma: a multi-cohort landmark, multistate and transcriptomic analysis of public "
        "individual-patient data” for consideration as a Research article in Journal of Translational Medicine. We suggest "
        "the Medical Bioinformatics section.",
        "Clinicians already know that a patient who has relapsed one year after diagnosis has a worse outlook than one who "
        "has not. What is less clear is whether the date of relapse adds anything, and published analyses disagree. Using "
        "individual-patient data from eight public cohorts, we show that both statements can be true. The timing of first "
        "progression is strongly associated with survival after progression (pooled hazard ratio 0.73 per six months of "
        "later progression in 1126 patients), yet at a month-12 visit it adds nothing to prediction because almost all "
        "progressed patients have progressed within the preceding few months. Its contribution becomes measurable as the "
        "landmark moves later.",
        "We believe the work fits the journal for three reasons. It joins clinical follow-up information with tumour "
        "expression, BRCA1/2, HRD and CCNE1 data from the same patients, and shows that these baseline molecular features "
        "relate to long-term outcome but add little to short-term prediction once progression status is known. It uses "
        "single-cell data to show that the mesenchymal expression signal of early progression comes mainly from fibroblasts "
        "and disappears after adjustment for stromal content. And it reports transport of predicted risks between cohorts "
        "with internal–external cross-validation, prediction intervals and recalibration, following TRIPOD+AI.",
        "Issues relating to journal policies. The study uses only public, de-identified data, so ethics approval and consent "
        "were not required. The prespecified analysis and all post hoc analyses are labelled as such in the manuscript and "
        "in Additional file 1. The cohorts predate maintenance PARP inhibitors, which we state as a limitation. An AI coding "
        "assistant was used for analysis code and for drafting and editing the text, as documented in the Methods. All code "
        "and intermediate results are deposited in a public repository.",
        "The authors declare that they have no competing interests. All authors have approved the manuscript for "
        "submission. The content has not been published and is not under consideration elsewhere.",
        "Yours sincerely,",
        "Ye Chen\nDepartment of Gynecology and Biomedical Innovation Center, The Sixth Affiliated Hospital, Sun Yat-sen "
        "University, 26 Yuancun Erheng Road, Tianhe District, Guangzhou 510655, Guangdong, China\n"
        "E-mail: cheny2599@mail.sysu.edu.cn",
    ]:
        doc.add_paragraph(text).paragraph_format.space_after = Pt(10)
    doc.save(PACK / "Cover_letter.docx")


if __name__ == "__main__":
    PACK.mkdir(exist_ok=True)
    {"main": main_document, "additional": additional_document, "letter": letter_and_figures}[MODE]()
    print("built", MODE)

# Progression status and progression timing in advanced serous tubo-ovarian carcinoma

Analysis code and intermediate results for the manuscript "Progression Status and Progression Timing as Updated
Prognostic Information in Advanced Serous Tubo-Ovarian Carcinoma".

All input data are public. They are not stored here; the download scripts fetch them and record URLs and checksums.

## Layout

The folder names are kept as used by the scripts, which locate each other by relative path.

| Folder | Content |
|---|---|
| `scripts/`, `results/`, `figures/` | Prespecified month-12 landmark analysis with equations fitted in TCGA (`Analysis_plan.md`). |
| `同行预审与病理补充/` | First set of post hoc checks: progression-only benchmark, grade subsets, TCGA slide stroma, APOLLO-2. |
| `修订稿/` | Second set of post hoc checks: two-stage bootstrap, adjacent landmarks, additional cohorts. |
| `修订稿_v2/` | Analyses reported in the submitted manuscript: internal-external cross-validation, multi-landmark supermodel, illness-death model, post-progression survival, expression scores, BRCA1/2, HRD, CCNE1, single-cell attribution, and the manuscript build scripts. |

## Reproduce

Python 3.12. Install `requirements.txt`, then run in this order from the repository root:

```text
python scripts/download_inputs.py
python scripts/prepare_data.py
python scripts/analyze.py
python scripts/status_benchmark.py
python scripts/plot.py
python 同行预审与病理补充/peer_review_checks.py
python 修订稿/fetch_revision_sources.py
python 修订稿/revision_analysis.py
python 修订稿_v2/fetch_v2_sources.py
python 修订稿_v2/v2_clinical.py
python 修订稿_v2/v2_molecular.py
python 修订稿_v2/v2_singlecell_pseudobulk.py
python 修订稿_v2/v2_genomics.py
python 修订稿_v2/v2_figures.py
python 修订稿_v2/build_v2_manuscript.py
python 修订稿_v2/build_jtm_submission.py main
```

Some inputs of the first two post hoc folders (GDC pathology JSON, APOLLO-2 workbook, TCGA slide files) were retrieved
manually; their URLs and checksums are in the `source_manifest.csv` files of those folders.

## License

MIT, see `LICENSE`. The source data keep the licences of their providers.

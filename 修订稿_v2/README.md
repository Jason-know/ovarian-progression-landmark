# Second revision: progression status, progression timing and molecular analyses

This folder holds the second revision. The first draft (`../Ovarian_Cancer_Landmark_Manuscript.docx`) and the first
revision (`../修订稿/`) are unchanged and are still needed as inputs.

## Contents

- `Ovarian_Cancer_Manuscript_v2.docx`: revised manuscript with tables, figures and supplementary material.
- `审稿意见逐条回应.md`: point-by-point response to the review comments.
- `v2_manuscript.txt`: manuscript text and references. Tables are generated from `results/`.
- `results/`: CSV output of every new analysis. `figures/`: main Figures 2 to 5 (PDF and PNG).
- `qa/`: PDF rendering and page images used to check layout.

## Reproduce

Run from any directory. The original `scripts/`, `raw/`, `results/`, `修订稿/` and `同行预审与病理补充/` folders must be
present next to this folder.

```text
python "修订稿_v2/fetch_v2_sources.py"     # public downloads, writes results/v2_source_manifest.csv
python "修订稿_v2/v2_clinical.py"          # post-progression survival, landmarks, illness-death, IECV, recalibration
python "修订稿_v2/v2_molecular.py"         # expression scores, BRCA1/2, transcriptome and hallmark enrichment
python "修订稿_v2/v2_singlecell_pseudobulk.py"  # GSE165897 pseudobulk by patient and cell class
python "修订稿_v2/v2_genomics.py"         # HRD, CCNE1, stroma-adjusted enrichment, single-cell attribution
python "修订稿_v2/v2_figures.py"
python "修订稿_v2/build_v2_manuscript.py"
```

Additional packages beyond `../requirements.txt`: `patsy` (spline bases), `gseapy` (pre-ranked
enrichment), `requests`.

## Before submission

TCGA copy-number and tumour-level features are fetched from cBioPortal by `v2_genomics.py` on first run and cached in `data/tcga_ov_pancan_genomics.json`.

The repository URL and the AI-use statement contain
bracketed placeholders that the authors must complete. The TRIPOD+AI table is condensed and the journal's full form
is still required.

# Peer review checks and pathology supplement

The review report reflects three independent AI reviewer perspectives. Named human experts in the report are prospective candidates only. No manuscript was sent to them.

## Files

- `peer_review_checks.py` runs the post hoc clinical benchmark, progression-only subgroup analysis, conditional censoring-weight check, grade-restricted check, and pathology analyses.
- `plot_review.py` draws the report figure. `build_review_report.py` creates the review document.
- `results/*.csv` records the paired estimates, failure to improve the model with TCGA stroma, and exploratory APOLLO-2 transport. `results/source_manifest.csv` records data provenance.
- `data/TCGA_OV_GDC_pathology.json` contains public GDC slide QC attributes. `data/APOLLO2_MOESM5.xlsx` contains the published APOLLO-2 pathology and outcomes. `data/SIGNPOST_MOESM3.xlsx` and `data/PTRC_HGSOC.xlsx` were inspected for data feasibility.

## Reproduction

Extract this supplement into the original study directory so that its subdirectory `同行预审与病理补充` sits alongside the original `scripts`, `raw` and `results` directories. Use the dependencies listed in the original `requirements.txt`.

```text
python "同行预审与病理补充/peer_review_checks.py"
python "同行预审与病理补充/plot_review.py"
python "同行预审与病理补充/build_review_report.py"
```

Analyses use the existing 12-month landmark and predict death within the next 12 months. The TCGA pathology comparison is internal and post hoc. The paired bootstrap conditions on predictions, including out-of-fold predictions where specified. APOLLO-2 applies the original locked clinical equations without refitting. The large uncertainty in this small external pathology cohort precludes a reliable calibration claim.

# Revised landmark manuscript and pathology analyses

The revised Word article incorporates peer-review changes, post hoc clinical comparisons, additional cohort analyses and two public TCGA H&E slide fields. The initial manuscript remains available in the study directory.

## Reproduce

Use the original study directory, with its `scripts/`, `figures/`, `raw/` and `results/`, and the earlier review supplement under `同行预审与病理补充/`. Place this revision folder next to those directories. The previous research and review archives provide these inputs.

The original requirements still apply. Rendering the H&E fields additionally uses openslide-python 1.4.7 and openslide-bin 4.0.1.2. The new source SVS and ExpressionSet files are referenced by URL and checksum in `results/source_manifest.csv`. Both raw SVS files are also stored in the shared workspace's revision data folder. The revision archive contains only the analysis-specific extracts of the ExpressionSet clinical data and omits administrative dates and unrelated source metadata.

```text
python "修订稿/fetch_revision_sources.py"
python "修订稿/revision_analysis.py"
python "修订稿/revision_figures.py"
python "修订稿/build_revised_manuscript.py"
```

The selected GDC image records are rechecked against clinical analysis eligibility and slide-level quality-control records. Figure 3 fields come from the source image at level 0. The PNG fields contain the cropped colors as read. The figure overlays scanner-derived scale bars. `results/slide_provenance.csv` gives GDC file identifiers, recorded stroma percentages and field coordinates.

## Analysis scope

The original 850-patient development and primary external evaluation is unchanged. Paired status-only benchmarks, separate progression-group censoring checks, subgroup estimates and exploratory pathology analyses use the earlier peer-review output. New sensitivity analyses resample development and evaluation patients, use neighboring landmark times, and apply predefined clinical equations or reduced models to UVA, Curie and OVCAD.

Curie lacks individual age and residual size, and OVCAD records any residual disease rather than a 1 cm dichotomy. OVCAD and UKE share a recruiting institution. Each source is therefore presented separately with its available feature set and neither replaces nor is pooled into the three primary validation cohorts.

Pathology visualizations cannot be used to infer clinical risk from a single section. The TCGA quantitative analysis uses the patient-median clinical biospecimen QC score across valid primary slides. APOLLO-2 is an exploratory clinical transport assessment in a separate pathology-reviewed cohort.

"""Check the core review deliverables and assemble the reproducibility archive."""
from pathlib import Path
import csv
import gzip
import zipfile

import pandas as pd

from download_inputs import digest

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = (
    "Figure_1_survival", "Figure_2_landmark_evidence", "Figure_3_predictor_diagnostics",
    "Figure_4_roc", "Figure_5_calibration_and_decision_curves", "Figure_6_subgroups",
    "Figure_S1_nomogram", "Figure_S2_external_timing", "Figure_S3_sensitivity", "Figure_S4_cohort_flow",
)


def source_inputs():
    paths = []
    with (RESULTS / "sources.csv").open(encoding="utf-8", newline="") as metadata:
        for row in csv.DictReader(metadata):
            path = ROOT / "raw" / row["file"]
            assert digest(path) == row["sha256"], path
            if "series_matrix" in path.name:
                with gzip.open(path, "rb") as stream:
                    while stream.read(2**20):
                        pass
            paths.append(path)
    return paths


def verify_results():
    performance = pd.read_csv(RESULTS / "model_metrics.csv")
    primary = performance.query("variant == 'Primary' and model == 'Status'")
    assert primary.n.sum() == 850 and primary.deaths.sum() == 142
    assert list(primary.cohort) == ["TCGA", "AOCS", "UKE", "ICON7"]
    assert primary.loc[primary.cohort.ne("TCGA"), "oe"].gt(1).all()
    assert performance.query("variant == 'No age' and model == 'Status' and cohort == 'Japan'").n.iloc[0] == 244
    comparisons = pd.read_csv(RESULTS / "model_contrasts.csv").query("variant == 'Primary' and cohort != 'TCGA'")
    status = comparisons.loc[comparisons.contrast.eq("Status vs baseline")]
    timing = comparisons.loc[comparisons.contrast.eq("Timing vs status")]
    assert len(status) == len(timing) == 3
    assert status.lo.gt(0).all() and timing.lo.lt(0).all() and timing.hi.gt(0).all()
    return primary[["cohort", "n", "deaths", "auc"]]


def package():
    inputs = source_inputs()
    results = verify_results()
    names = ["README.md", "Analysis_plan.md", "requirements.txt", "审阅说明.md"]
    paths = [ROOT / name for name in names]
    paths += sorted(p for p in (ROOT / "scripts").iterdir() if p.suffix in {".py", ".txt"})
    paths += sorted(p for p in RESULTS.iterdir() if p.suffix in {".csv", ".json"})
    paths += [ROOT / "figures" / (name + ext) for name in FIGURES for ext in (".pdf", ".png")]
    paths += inputs

    manuscript = ROOT / "Ovarian_Cancer_Landmark_Manuscript.docx"
    with zipfile.ZipFile(manuscript) as docx:
        media = [name for name in docx.namelist() if name.startswith("word/media/")]
        assert len(media) == len(FIGURES) and docx.testzip() is None

    output = ROOT / "Ovarian_Cancer_Research_Package.zip"
    temporary = output.with_suffix(".zip.part")
    with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=5) as archive:
        for path in paths:
            method = zipfile.ZIP_STORED if path in inputs or path.suffix in {".pdf", ".png"} else zipfile.ZIP_DEFLATED
            archive.write(path, path.relative_to(ROOT), compress_type=method)
    with zipfile.ZipFile(temporary) as archive:
        assert archive.testzip() is None
        print("Verified archive entries:", len(archive.namelist()))
    temporary.replace(output)
    print(results.round(3).to_string(index=False))
    print(output, output.stat().st_size)


if __name__ == "__main__":
    package()

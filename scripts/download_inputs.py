"""Download the six immutable public data inputs and verify their checksums."""
from pathlib import Path
from urllib.request import Request, urlopen
import csv
import hashlib
import shutil
import time

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw"
RESULTS = ROOT / "results"
URLS = {
    "TCGA-CDR-SupplementalTableS1.xlsx":
        "https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81",
    "TCGA_firehose.tar.gz":
        "https://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/OV/20160128/"
        "gdac.broadinstitute.org_OV.Merge_Clinical.Level_1.2016012800.0.0.tar.gz",
    "tothillClinical.RData":
        "https://bioinformatics.mdanderson.org/Supplements/ResidualDisease/RDataObjects/tothillClinical.RData",
    "UKE_supplement.zip":
        "https://static-content.springer.com/esm/art%3A10.1186%2Fs13073-026-01727-6/"
        "MediaObjects/13073_2026_1727_MOESM1_ESM.zip",
    "GSE140082_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140082/matrix/GSE140082_series_matrix.txt.gz",
    "GSE32062-GPL6480_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32062/matrix/GSE32062-GPL6480_series_matrix.txt.gz",
}


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def run():
    RAW.mkdir(exist_ok=True)
    with (RESULTS / "sources.csv").open(encoding="utf-8", newline="") as handle:
        checksums = {row["file"]: row["sha256"] for row in csv.DictReader(handle)}
    for name, url in URLS.items():
        path = RAW / name
        checksum = checksums[name]
        if path.is_file() and digest(path) == checksum:
            print(f"Verified {name}")
            continue
        tmp = path.with_name(name + ".part")
        for attempt in range(3):
            try:
                req = Request(url, headers={"User-Agent": "ovarian-landmark-reproducibility/1.0"})
                with urlopen(req, timeout=90) as incoming, tmp.open("wb") as out:
                    shutil.copyfileobj(incoming, out)
                if digest(tmp) != checksum:
                    raise ValueError(f"Downloaded content changed for {name}")
                tmp.replace(path)
                print(f"Downloaded and verified {name}")
                break
            except (OSError, ValueError):
                tmp.unlink(missing_ok=True)
                if attempt == 2:
                    raise
                time.sleep(attempt + 2)


if __name__ == "__main__":
    run()

"""Download the additional public inputs used by the second revision."""
from pathlib import Path
import csv
import hashlib
import json

import requests

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
ZENODO = "https://zenodo.org/records/23125873/files/{}?download=1"
FILES = {
    "GSE9891_eset.rda": ZENODO.format("GSE9891_eset.rda"),
    "GSE32062.GPL6480_eset.rda": ZENODO.format("GSE32062.GPL6480_eset.rda"),
    "GPL14951_full.txt": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL14951&targ=self&form=text&view=data",
    "hgnc_complete_set.txt": "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
    "h.all.v2023.2.Hs.symbols.gmt":
        "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt",
    "TCGA.HRD_withSampleID.txt": "https://api.gdc.cancer.gov/data/66dd07d7-6366-4774-83c3-5ad1e22b177e",
    "GSE165897_cellInfo_HGSOC.tsv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165897/suppl/GSE165897_cellInfo_HGSOC.tsv.gz",
    "GSE165897_UMIcounts_HGSOC.tsv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165897/suppl/GSE165897_UMIcounts_HGSOC.tsv.gz",
}
CBIO = "https://www.cbioportal.org/api"


def fetch(name, url):
    path = DATA / name
    if not path.exists():
        with requests.get(url, stream=True, timeout=300, headers={"User-Agent": "Mozilla/5.0"}) as r:
            r.raise_for_status()
            with path.open("wb") as handle:
                for chunk in r.iter_content(1024 ** 2):
                    handle.write(chunk)
    with path.open("rb") as handle:
        return {"file": name, "url": url, "bytes": path.stat().st_size,
                "sha256": hashlib.file_digest(handle, "sha256").hexdigest()}


def brca():
    """Germline and somatic BRCA1/2 calls of the TCGA 2011 marker paper via cBioPortal."""
    path = DATA / "tcga_ov_pub_brca.json"
    if not path.exists():
        study = "ov_tcga_pub"
        muts = requests.post(f"{CBIO}/molecular-profiles/{study}_mutations/mutations/fetch?projection=DETAILED",
                             json={"sampleListId": f"{study}_sequenced", "entrezGeneIds": [672, 675]}, timeout=120)
        muts.raise_for_status()
        samples = requests.get(f"{CBIO}/sample-lists/{study}_sequenced", timeout=120)
        samples.raise_for_status()
        path.write_text(json.dumps({"mutations": muts.json(), "sequenced": samples.json()["sampleIds"]}),
                        encoding="utf-8")
    with path.open("rb") as handle:
        return {"file": path.name, "url": f"{CBIO} (ov_tcga_pub, BRCA1/BRCA2, sequenced sample list)",
                "bytes": path.stat().st_size, "sha256": hashlib.file_digest(handle, "sha256").hexdigest()}


if __name__ == "__main__":
    DATA.mkdir(exist_ok=True)
    rows = []
    for name, url in FILES.items():
        try:
            rows.append(fetch(name, url))
            print("ok", name, flush=True)
        except Exception as exc:  # report and continue so that one mirror failure is visible
            print("FAILED", name, repr(exc), flush=True)
    rows.append(brca())
    with (HERE / "results" / "v2_source_manifest.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["file", "url", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(rows)

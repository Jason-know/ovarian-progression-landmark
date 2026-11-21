"""Discover available clinical annotations without loading expression matrices."""
from pathlib import Path
from collections import defaultdict
import csv
import gzip
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw"


def parse_matrix(path):
    records = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith(("!Sample_geo_accession", "!Sample_title", "!Sample_characteristics_ch1")):
                continue
            key, *values = next(csv.reader([line], delimiter="\t"))
            if records is None:
                records = [{} for _ in values]
            for record, value in zip(records, values, strict=True):
                if key == "!Sample_characteristics_ch1" and ": " in value:
                    field, value = value.split(": ", 1)
                    record[field.strip().lower()] = value.strip()
                elif key != "!Sample_characteristics_ch1":
                    record[key.removeprefix("!Sample_")] = value
    return pd.DataFrame(records)


def parse_soft(path):
    records = {}
    rec = None
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("^SAMPLE = "):
                accession = line.split(" = ", 1)[1].strip()
                rec = records.setdefault(accession, {"geo_accession": accession})
            elif rec is not None and line.startswith("!Sample_characteristics_ch1 = ") and ": " in line:
                key, value = line.split(" = ", 1)[1].strip().split(": ", 1)
                rec[key.strip().lower()] = value
    return pd.DataFrame(records.values())


if __name__ == "__main__":
    for path in sorted(RAW.glob("*matrix.txt.gz")) + sorted(RAW.glob("*soft.txt")):
        try:
            data = parse_matrix(path) if path.suffix == ".gz" else parse_soft(path)
            print(path.name, data.shape)
            for key in data.columns[2:]:
                print(key, data[key].nunique(), data[key].value_counts(dropna=False).head(12).to_dict())
        except (OSError, ValueError, EOFError) as exc:
            print(path.name, exc)
    print("UKE:")
    import zipfile
    for path in RAW.glob("UKE*.zip"):
        with zipfile.ZipFile(path) as z:
            print(z.namelist())
    p = RAW / "TCGA-CDR-SupplementalTableS1.xlsx"
    print("TCGA CDR:", pd.ExcelFile(p).sheet_names)
    f = pd.read_excel(p, sheet_name=0)
    print(f.shape, list(f.columns), f.head(2).to_string(index=False))
    p = RAW / "tcga_ov_api.json"
    if p.exists() and p.stat().st_size:
        data = json.loads(p.read_text(encoding="utf-8"))
        print("TCGA cBio attributes:", sorted({r["clinicalAttributeId"] for r in data}))

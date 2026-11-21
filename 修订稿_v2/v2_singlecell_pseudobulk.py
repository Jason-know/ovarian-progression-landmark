"""Pseudobulk the GSE165897 single-cell counts by patient, treatment phase and cell class (streamed)."""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA, OUT = HERE / "data", HERE / "results"
CLASS = {"T-cells": "T/NK", "NK": "T/NK", "ILC": "T/NK", "B-cells": "B/plasma", "Plasma-cells": "B/plasma",
         "Macrophages": "Myeloid", "DC-1": "Myeloid", "DC-2": "Myeloid", "pDC": "Myeloid", "Mast-cells": "Myeloid",
         "CAF-1": "Fibroblast", "CAF-2": "Fibroblast", "CAF-3": "Fibroblast", "Mesothelial": "Mesothelial",
         "Endothelial": "Endothelial"}

if __name__ == "__main__":
    info = pd.read_csv(DATA / "GSE165897_cellInfo_HGSOC.tsv.gz", sep="\t")
    info["cell_class"] = np.where(info.cell_type.eq("EOC"), "Malignant", info.cell_subtype.map(CLASS))
    assert info.cell_class.notna().all()
    info["group"] = info.patient_id + "|" + info.treatment_phase + "|" + info.cell_class
    genes, sums = [], []
    with gzip.open(DATA / "GSE165897_UMIcounts_HGSOC.tsv.gz", "rt") as handle:
        cells = handle.readline().rstrip("\n").split("\t")[1:]
        assert cells == list(info.cell)
        codes, labels = pd.factorize(info.group)
        for n, line in enumerate(handle):
            gene, rest = line.split("\t", 1)
            values = np.fromstring(rest, dtype=np.float32, sep="\t")
            genes.append(gene)
            sums.append(np.bincount(codes, weights=values, minlength=len(labels)))
            if n % 5000 == 0:
                print(n, flush=True)
    pd.DataFrame(np.vstack(sums), index=genes, columns=labels).to_csv(OUT / "sc_pseudobulk_counts.csv.gz")
    info.groupby(["patient_id", "treatment_phase", "cell_class"]).size().rename("cells").reset_index().to_csv(
        OUT / "sc_cell_counts.csv", index=False)
    print("done", len(genes))

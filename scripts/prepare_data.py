"""Harmonize public, de-identified clinical records and audit eligibility."""
from pathlib import Path
import hashlib
import io
import tarfile
import warnings
import zipfile

import numpy as np
import pandas as pd
import rdata

from audit_sources import parse_matrix

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw"
OUT = ROOT / "results"
MONTH = 365.25 / 12
SOURCES = {
    "TCGA": "TCGA-CDR-SupplementalTableS1.xlsx",
    "AOCS": "tothillClinical.RData",
    "UKE": "UKE_supplement.zip",
    "ICON7": "GSE140082_series_matrix.txt.gz",
    "Japan": "GSE32062-GPL6480_series_matrix.txt.gz",
}
RESIDUAL_SOURCE = "TCGA_firehose.tar.gz"


def numeric(s):
    return pd.to_numeric(s, errors="coerce").astype(float)


def common(cohort, sid, age, stage, rd, ot, oe, pt, pe, serous, treatment, grade=None):
    stage = stage.astype("string").str.upper().str.replace("STAGE ", "", regex=False)
    return pd.DataFrame({
        "cohort": cohort, "id": sid.astype(str), "age": numeric(age),
        "stage": stage, "iv": stage.eq("IV").astype(float),
        "advanced": stage.str.fullmatch(r"III[A-C]?|IV").fillna(False),
        "rd": numeric(rd), "os": numeric(ot), "death": numeric(oe),
        "pfs": numeric(pt), "progression_event": numeric(pe),
        "serous": serous, "treatment": treatment, "grade": grade,
    })


def aocs():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        s = rdata.read_rda(RAW / SOURCES["AOCS"])["tothillClinical"]
    return common(
        "AOCS", s["GEO.ID"], s.Age, s.Stage, s.ResidDisease.map({"nil": 0, "<1": 0, ">1": 1}),
        s.MosToDeath, s.Status.map({"D": 1, "D*": 1, "PF": 0, "R": 0}),
        s.MosToRelapse, s.Status.map({"D": 1, "D*": 1, "R": 1, "PF": 0}),
        s.ClinicalType.eq("MAL") & s.HistologicSubtype.eq("Ser"),
        pd.Series(np.where(s.Pltx.eq("Y") & s.Neo.eq("N"), "Platinum", "Other"), index=s.index),
        s.Grade,
    )


def uke():
    with zipfile.ZipFile(RAW / SOURCES["UKE"]) as outer:
        nested = next(n for n in outer.namelist() if n.endswith("/Table_S1.zip"))
        with zipfile.ZipFile(io.BytesIO(outer.read(nested))) as inner:
            with inner.open("Tables_S1/Table_S1b.UKE_metadata.tsv") as data:
                s = pd.read_csv(data, sep="\t")
    s = s.loc[s.patient_id.notna() & s.recurrent_tumor.eq("no")].copy()
    clinical = ["age", "OS_time", "OS_event", "PFS_time", "PFS_event", "therapy",
                "FIGO_stage", "operation_result", "histological_type"]
    conflicting = s.groupby("patient_id")[clinical].nunique().gt(1).any(axis=1)
    if conflicting.any():
        raise ValueError("Conflicting clinical endpoints across UKE samples of one patient")
    s = s.sort_values(["patient_id", "tumor_sample"]).groupby("patient_id", as_index=False).first()
    return common(
        "UKE", s.patient_id.astype(int), s.age, s.FIGO_stage,
        s.operation_result.map({"no tumor": 0, "residual tumor <1cm": 0, "residual tumor >1cm": 1}),
        s.OS_time, s.OS_event, s.PFS_time, s.PFS_event, s.histological_type.eq("serous"),
        s.therapy.map({"standard platinum": "Platinum", "platinum+bevacizumab": "Bevacizumab"}).fillna("Other"),
    )


def icon7():
    s = parse_matrix(RAW / SOURCES["ICON7"])
    return common(
        "ICON7", s.geo_accession, s.age, s.figo_stage,
        s.debulking_status.map({"OPTIMAL": 0, "SUB-OPTIMAL": 1}),
        numeric(s.final_ostm) / MONTH, s.final_osid, numeric(s.final_pfstm) / MONTH, s.final_pfsid,
        s["histology.serous"].eq("serous"),
        s.treatment.map({"standard": "Platinum", "bevacizumab": "Bevacizumab"}), s.newgrade,
    )


def japan():
    s = parse_matrix(RAW / SOURCES["Japan"])
    return common(
        "Japan", s.geo_accession, pd.Series(np.nan, index=s.index), s.stage,
        s["surgery status"].str.lower().map({"optimal": 0, "suboptimal": 1}),
        s["os (m)"], s["death (1)"], s["pfs (m)"], s["rec (1)"],
        s.tissue.str.contains("serous"), pd.Series("Platinum", index=s.index), s.grading,
    )


def tcga():
    s = pd.read_excel(RAW / SOURCES["TCGA"], sheet_name="TCGA-CDR")
    s = s.loc[s.type.eq("OV") & s.Redaction.isna()].copy()
    with tarfile.open(RAW / RESIDUAL_SOURCE) as archive:
        table = next(n for n in archive.getnames() if n.endswith("/OV.clin.merged.txt"))
        raw = pd.read_csv(archive.extractfile(table), sep="\t", index_col=0, header=None, dtype=str)
    residual = raw.loc[["patient.bcr_patient_barcode", "patient.tumor_residual_disease"]].T
    residual.columns = ["id", "size"]
    residual["id"] = residual.id.str.upper()
    assert not residual.id.dropna().duplicated().any()
    sized = residual.dropna(subset=["id"]).set_index("id")["size"].str.lower()
    size = s.bcr_patient_barcode.map(sized)
    rd = size.map({"no macroscopic disease": 0, "1-10 mm": 0, "11-20 mm": 1, ">20 mm": 1})
    print("TCGA source residual size categories:", size.value_counts(dropna=False).to_dict())
    # PFI events in CDR can include death with tumor. An individual known to be
    # alive at the landmark cannot have such a death at or before the landmark.
    result = common(
        "TCGA", s.bcr_patient_barcode, s.age_at_initial_pathologic_diagnosis,
        s.clinical_stage, rd,
        numeric(s["OS.time"]) / MONTH, s.OS, numeric(s["PFI.time"]) / MONTH, s.PFI,
        s.histological_type.eq("Serous Cystadenocarcinoma"),
        pd.Series("Unavailable", index=s.index), s.histological_grade,
    )
    result["documented_event"] = s.new_tumor_event_type.isin(["Recurrence", "Progression of Disease"])
    return result


def eligible(frame, landmark=12, require=("age", "rd"), inclusive=True):
    masks = {
        "Available source patients": pd.Series(True, index=frame.index),
        "Stage III–IV invasive serous": frame.advanced & frame.serous,
        "Treatment eligible when reported": ~frame.treatment.eq("Other"),
        "Valid baseline and endpoint fields": (
            frame[list(require)].notna().all(axis=1) & frame.os.ge(0) & frame.pfs.ge(0)
            & frame.death.isin([0, 1]) & frame.progression_event.isin([0, 1])),
        "Alive and under OS observation after landmark": frame.os.gt(landmark),
        "Progression status established at landmark": (
            frame.pfs.ge(landmark) | (frame.progression_event.eq(1) & frame.pfs.lt(landmark))),
    }
    keep = pd.Series(True, index=frame.index)
    counts = {}
    for name, mask in masks.items():
        keep &= mask
        counts[name] = int(keep.sum())
    data = frame.loc[keep].copy()
    if not inclusive:
        data = data.loc[~(data.progression_event.eq(1) & data.pfs.eq(landmark))].copy()
    predates = data.pfs.le(landmark) if inclusive else data.pfs.lt(landmark)
    data["early"] = (predates & data.progression_event.eq(1)).astype(int)
    data["duration"] = (landmark - data.pfs).where(data.early.eq(1), 0)
    data["time"] = data.os - landmark
    assert not data.id.duplicated().any()
    assert data.time.gt(0).all() and data.duration.ge(0).all()
    return data.reset_index(drop=True), counts


def load_all():
    return {f.__name__.upper().replace("JAPAN", "Japan"): f() for f in (aocs, uke, icon7, japan, tcga)}


if __name__ == "__main__":
    manifests = []
    for cohort, name in list(SOURCES.items()) + [("TCGA", RESIDUAL_SOURCE)]:
        path = RAW / name
        with path.open("rb") as handle:
            manifests.append({"cohort": cohort, "file": name, "bytes": path.stat().st_size,
                              "sha256": hashlib.file_digest(handle, "sha256").hexdigest()})
    pd.DataFrame(manifests).to_csv(OUT / "sources.csv", index=False)
    for name, frame in load_all().items():
        needed = {"Japan": ("rd",)}.get(name, ("age", "rd"))
        d, counts = eligible(frame, require=needed)
        print(name, counts, "early:", int(d.early.sum()), "deaths by 12 and 24:",
              [int((d.death.eq(1) & d.time.le(t)).sum()) for t in (12, 24)],
              "observable controls:", [int(d.time.gt(t).sum()) for t in (12, 24)],
              "PFS after OS >1 month:", int(d.pfs.gt(d.os + 1).sum()),
              "progression at landmark:", int((d.early.eq(1) & d.duration.eq(0)).sum()))

"""Portable source ingestion. No credentials, downloads, or access registry."""
from __future__ import annotations

import csv
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import tarfile

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
MOVER_COLUMNS = {
    "patient_information.csv": {"LOG_ID", "MRN", "BIRTH_DATE", "AN_START_DATETIME",
        "AN_STOP_DATETIME", "PRIMARY_ANES_TYPE_NM", "PRIMARY_PROCEDURE_NM",
        "HOSP_ADMSN_TIME", "HOSP_DISCH_TIME", "SEX", "ASA_RATING_C"},
    "patient_labs.csv": {"LOG_ID", "MRN", "Lab Code", "Measurement Units",
        "Observation Value", "Collection Datetime"},
    "patient_visit.csv": {"LOG_ID", "mrn", "diagnosis_code", "dx_name"},
}
VITAL_COLUMNS = {
    "clinical_data.csv": {"caseid", "subjectid", "opname", "age", "ane_type",
        "anestart", "aneend", "adm", "dis", "asa", "sex", "dx"},
    "lab_data.csv": {"caseid", "name", "result", "dt"},
}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def initialize_output(output, *, demo=False):
    """Fresh output only; never mix a real-data run into the distributable repo."""
    target = Path(output).expanduser().resolve()
    if target.exists():
        raise ValueError("Output already exists; select a fresh directory")
    if not demo and target.is_relative_to(REPO):
        raise ValueError("Real-data output must be outside the code repository")
    os.umask(0o077)
    target.mkdir(parents=True, mode=0o700)
    for child in (".private", "results", "evidence"):
        (target / child).mkdir(mode=0o700)
    return target


def validate_header(stream, required):
    """Reject duplicate or missing headers without exposing raw row text."""
    position = stream.tell()
    text = io.TextIOWrapper(stream, encoding="utf-8-sig", newline="")
    try:
        names = next(csv.reader(text), [])
    finally:
        text.detach()
        stream.seek(position)
    if len(names) != len(set(names)) or not required.issubset(names):
        raise ValueError("Source CSV schema does not match the frozen release")


def extract_mover(archive_path, output):
    """Read three regular tar members with extractfile; never extract paths."""
    destinations = {"patient_information.csv": "mover_operations.parquet",
                    "patient_labs.csv": "mover_creatinine.parquet",
                    "patient_visit.csv": "mover_visit.parquet"}
    receipt = {"archive_sha256": sha(archive_path), "tables": {}}
    with tarfile.open(archive_path, "r:*") as archive:
        selected = {}
        for member in archive.getmembers():
            name = PurePosixPath(member.name).name
            if name not in destinations:
                continue
            if not member.isfile() or name in selected:
                raise ValueError("MOVER source has a duplicate or nonregular required member")
            selected[name] = member
        if set(selected) != set(destinations):
            raise ValueError("MOVER archive lacks one of the three required CSV tables")
        for name, member in selected.items():
            stream = archive.extractfile(member)
            if stream is None:
                raise ValueError("Cannot read required archive member")
            with stream:
                validate_header(stream, MOVER_COLUMNS[name])
                retained, scanned = [], 0
                for chunk in pd.read_csv(stream, dtype=str, keep_default_na=False,
                                         encoding="utf-8-sig", chunksize=500000):
                    scanned += len(chunk)
                    if name == "patient_labs.csv":
                        chunk = chunk[chunk["Lab Code"].eq("2160-0")
                                      & chunk["Measurement Units"].str.lower().eq("mg/dl")]
                    retained.append(chunk)
                if not retained:
                    raise ValueError("Required source CSV contains no records")
                frame = pd.concat(retained, ignore_index=True)
            target = output / ".private" / destinations[name]
            frame.to_parquet(target, index=False)
            target.chmod(0o600)
            receipt["tables"][name] = {"rows_scanned": scanned, "rows_retained": len(frame),
                                      "derived_sha256": sha(target)}
    return receipt


def validate_vital(directory):
    receipt = {}
    for name, columns in VITAL_COLUMNS.items():
        path = Path(directory) / name
        with path.open("rb") as stream:
            validate_header(stream, columns)
        receipt[name] = sha(path)
    return receipt


def verify_source_checksums(mover_archive, vitaldb_directory, *, expected=None):
    """Bind formal execution to source-snapshot bytes using published hashes."""
    if expected is None:
        expected = json.loads((REPO / "config/source-checksums.json").read_text())
    sources = {"EPIC_EMR.tar.gz": Path(mover_archive)}
    for name in ("clinical_data.csv", "lab_data.csv", "clinical_parameters.csv", "lab_parameters.csv"):
        sources[name] = Path(vitaldb_directory) / name
    if set(expected) != set(sources):
        raise ValueError("Source-checksum manifest does not cover the fixed five source files")
    actual = {name: sha(path) for name, path in sources.items()}
    if actual != expected:
        raise ValueError("Source bytes differ from the verified study snapshot")
    return actual


def build_inputs(mover_archive, vitaldb_directory, output):
    import build_frozen_inputs as frozen
    frozen.P = output
    receipt = {"mover": extract_mover(mover_archive, output),
               "vitaldb": validate_vital(vitaldb_directory),
               "source_version_required": "MOVER EPIC snapshot used by the study; VitalDB clinical/labs 1.0.0"}
    pm = pd.read_csv(REPO / "config/procedure-map.csv", dtype=str, keep_default_na=False)
    dm = pd.read_csv(REPO / "config/diagnosis-map.csv", dtype=str, keep_default_na=False)
    for mapping in (pm, dm):
        if not set(mapping.decision) <= {"include", "exclude", "uncertain"}:
            raise ValueError("Invalid frozen mapping decision")
    frozen.run_dataset("mover", *frozen.build_mover(), pm, dm)
    frozen.run_dataset("vitaldb", *frozen.build_vital(Path(vitaldb_directory)), pm, dm)
    receipt["algorithm_config_sha256"] = {name: sha(REPO / "config" / name)
                                          for name in ("procedure-map.csv", "diagnosis-map.csv")}
    dump(output / "evidence/source-provenance.json", receipt)
    return receipt

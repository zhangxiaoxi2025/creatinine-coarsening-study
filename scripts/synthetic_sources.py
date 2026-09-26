"""Generate wholly fictional source-format fixtures; no study records."""
import csv
from datetime import datetime, timedelta
import io
from pathlib import Path
import tarfile

import pandas as pd

from portable_io import REPO


EXPECTED = {"main": 10, "parallel7": 11, "S1_no_diagnosis": 11,
            "S2_target_after5min": 10, "S3_exclude_extremes": 9,
            "S4_exclude_over24h": 9, "S5_exclude_conflicts": 9, "S6_all_postop_E48": 9}


def csv_bytes(rows):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def create(directory):
    directory = Path(directory)
    directory.mkdir(parents=True, mode=0o700)
    pm = pd.read_csv(REPO / "config/procedure-map.csv", dtype=str, keep_default_na=False)
    dm = pd.read_csv(REPO / "config/diagnosis-map.csv", dtype=str, keep_default_na=False)
    clinical, labs, visits, vital_clinical, vital_labs = [], [], [], [], []
    origin = datetime(2020, 1, 10, 12)
    def stamp(seconds):
        return (origin + timedelta(seconds=seconds)).strftime("%Y-%m-%d %H:%M:%S")
    for center in ("mover", "vitaldb"):
        procedure = pm.loc[pm.dataset.eq(center) & pm.decision.eq("include"), "term"].iloc[0]
        exclusion = dm.loc[dm.dataset.eq(center) & dm.decision.eq("exclude")].iloc[0]
        for i in range(13):
            key = f"SYNTHETIC_{center}_{i:03}"
            person = f"FICTIONAL_PERSON_{center}_{i:03}"
            end = 90000 if i == 8 else 3600
            discharge = end + 604800
            baseline, current = [("1.00", "1.30"), ("0.30", "0.45"), ("1.00", "1.10"), ("1.00", "1.50")][i % 4]
            if i == 7:
                baseline = "0.05"
            records = [(-200000 if i == 12 else -3600, baseline), (end + 600, current)]
            if i == 4:
                records.append((end + 300000, "1.40"))
            if i == 6:
                records.extend([(end + 600, "2.10"), (end + 1200, "1.25")])
            if i == 9:
                records.append((end + 60, "0.50"))
            if i == 0:
                records.append((-3600, baseline))  # Exact duplicate is removed.
                records.append((end + 7200, "not_numeric"))  # Raw quality exclusion.
            if center == "mover":
                row = {"LOG_ID": key, "MRN": person, "BIRTH_DATE": "90" if i == 11 else "50",
                       "AN_START_DATETIME": stamp(0), "AN_STOP_DATETIME": stamp(end),
                       "PRIMARY_ANES_TYPE_NM": "General", "PRIMARY_PROCEDURE_NM": procedure,
                       "HOSP_ADMSN_TIME": stamp(-604800), "HOSP_DISCH_TIME": stamp(discharge),
                       "SEX": "F", "ASA_RATING_C": "2"}
                clinical.append(row)
                if i == 10:
                    clinical.append(dict(row, LOG_ID=key + "_second"))
                for time, value in records:
                    labs.append({"LOG_ID": key, "MRN": person, "Lab Code": "2160-0",
                                 "Measurement Units": "mg/dL", "Observation Value": value,
                                 "Collection Datetime": stamp(time)})
                visits.append({"LOG_ID": key, "mrn": person,
                               "diagnosis_code": exclusion.code if i == 5 else "",
                               "dx_name": exclusion.term if i == 5 else "SYNTHETIC NONRENAL DIAGNOSIS"})
            else:
                row = {"caseid": key, "subjectid": person, "opname": procedure,
                       "age": ">89" if i == 11 else "50", "ane_type": "General",
                       "anestart": "0", "aneend": str(end), "adm": "-604800", "dis": str(discharge),
                       "asa": "2", "sex": "F", "dx": exclusion.term if i == 5 else "SYNTHETIC NONRENAL DIAGNOSIS"}
                vital_clinical.append(row)
                if i == 10:
                    vital_clinical.append(dict(row, caseid=key + "_second"))
                vital_labs.extend({"caseid": key, "name": "cr", "result": value, "dt": str(time)} for time, value in records)
    # Same value in other analytes/units must never enter the MOVER creatinine pool.
    labs.extend([dict(labs[0], **{"Lab Code": "OTHER"}),
                 dict(labs[0], **{"Measurement Units": "umol/L"})])
    archive = directory / "SYNTHETIC_EPIC_EMR.tar.gz"
    with tarfile.open(archive, "w:gz") as target:
        for name, rows in (("patient_information.csv", clinical), ("patient_labs.csv", labs), ("patient_visit.csv", visits)):
            content = csv_bytes(rows)
            member = tarfile.TarInfo("synthetic/" + name)
            member.size, member.mode, member.mtime = len(content), 0o600, 0
            target.addfile(member, io.BytesIO(content))
    archive.chmod(0o600)
    vital = directory / "vitaldb"
    vital.mkdir(mode=0o700)
    for name, rows in (("clinical_data.csv", vital_clinical), ("lab_data.csv", vital_labs)):
        path = vital / name
        path.write_bytes(csv_bytes(rows))
        path.chmod(0o600)
    return archive, vital

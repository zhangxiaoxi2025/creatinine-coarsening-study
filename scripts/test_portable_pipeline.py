"""Portable I/O and complete CLI tests using fictional source records only."""
from contextlib import redirect_stdout
import io
import json
import os
import subprocess
import sys
from pathlib import Path
import tarfile
import tempfile
import unittest

import numpy as np

from portable_io import REPO, extract_mover, initialize_output, validate_header, verify_source_checksums, sha
from run_pipeline import run


class PortablePipelineTests(unittest.TestCase):
    def setUp(self):
        work = REPO / ".synthetic-test-work"
        work.mkdir(exist_ok=True)
        self.temporary = tempfile.TemporaryDirectory(prefix="synthetic-unit-", dir=work)
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def test_entire_pipeline_denominators_privacy_and_repeatability(self):
        with redirect_stdout(io.StringIO()):
            first = run(output=self.root / "first", demo=True, repetitions=5)
            second = run(output=self.root / "second", demo=True, repetitions=5)
        self.assertEqual(first, second)
        self.assertEqual(first["status"], "SYNTHETIC_END_TO_END_PASS")
        firstdir, seconddir = self.root / "first", self.root / "second"
        self.assertEqual((firstdir / ".private").stat().st_mode & 0o777, 0o700)
        for path in (firstdir / ".private").rglob("*"):
            if path.is_file():
                self.assertEqual(path.stat().st_mode & 0o777, 0o600)
        receipt = json.loads((firstdir / "evidence/source-provenance.json").read_text())
        labs = receipt["mover"]["tables"]["patient_labs.csv"]
        self.assertEqual(labs["rows_scanned"] - labs["rows_retained"], 2)
        files = list((firstdir / "results").glob("*-analysis.json"))
        self.assertEqual(len(files), 16)
        for file in files:
            one = json.loads(file.read_text())
            two = json.loads((seconddir / "results" / file.name).read_text())
            for key in ("point_estimates", "intervals", "mapping_summary", "rng_initial_state", "rng_final_state"):
                self.assertEqual(one[key], two[key])
            prefix = file.name.removesuffix("-analysis.json")
            a = np.load(firstdir / ".private" / f"{prefix}-bootstrap-counts.npz")["replicate_counts"]
            b = np.load(seconddir / ".private" / f"{prefix}-bootstrap-counts.npz")["replicate_counts"]
            np.testing.assert_array_equal(a, b)
            self.assertEqual(a.shape, (5, 4, 4))
        for center in ("mover", "vitaldb"):
            extension = json.loads((firstdir / "results" / f"{center}-extension.json").read_text())
            self.assertEqual(extension["N"], 10)
            self.assertTrue(extension["no_new_interval_computed"])
            self.assertEqual(sum(sum(row) for row in extension["criterion_transition_4x4"]), 10)
            self.assertEqual(extension["fixed_mapping_strata"]["some_observed_postop_outside_E48"]["n"], 1)
        for file in (firstdir / "results").glob("*.json"):
            self.assertNotIn("FICTIONAL_PERSON_", file.read_text())
            self.assertNotIn("SYNTHETIC_mover_", file.read_text())
            self.assertNotIn("SYNTHETIC_vitaldb_", file.read_text())

    def test_missing_and_duplicate_headers_fail_closed(self):
        for raw in (b"caseid,result\na,1\n", b"caseid,name,name,result,dt\na,cr,cr,1,0\n"):
            with self.assertRaises(ValueError):
                validate_header(io.BytesIO(raw), {"caseid", "name", "result", "dt"})
        raw = io.BytesIO(b"\xef\xbb\xbfcaseid,name,result,dt\na,cr,1,0\n")
        validate_header(raw, {"caseid", "name", "result", "dt"})
        self.assertEqual(raw.tell(), 0)

    def test_duplicate_and_symlink_tar_members_rejected_without_extraction(self):
        for malicious in ("duplicate", "symlink"):
            archive = self.root / (malicious + ".tar.gz")
            with tarfile.open(archive, "w:gz") as stream:
                names = ["patient_information.csv", "patient_labs.csv", "patient_visit.csv"]
                if malicious == "duplicate":
                    names.append("other/patient_labs.csv")
                for name in names:
                    member = tarfile.TarInfo(name)
                    if malicious == "symlink" and name == "patient_labs.csv":
                        member.type, member.linkname = tarfile.SYMTYPE, "../outside.csv"
                        stream.addfile(member)
                    else:
                        member.size = 0
                        stream.addfile(member, io.BytesIO(b""))
            out = initialize_output(self.root / (malicious + "-output"), demo=True)
            with self.assertRaises(ValueError):
                extract_mover(archive, out)
            self.assertEqual(list((out / ".private").iterdir()), [])

    def test_source_checksum_contract_uses_all_five_files_and_rejects_changed_bytes(self):
        from synthetic_sources import create
        archive, vital = create(self.root / "checksum-synthetic")
        for name in ("clinical_parameters.csv", "lab_parameters.csv"):
            (vital / name).write_text("parameter,unit\nsynthetic,synthetic\n")
        expected = {"EPIC_EMR.tar.gz": sha(archive)}
        expected.update({p.name: sha(p) for p in vital.glob("*.csv")})
        self.assertEqual(verify_source_checksums(archive, vital, expected=expected), expected)
        with self.assertRaises(ValueError):
            verify_source_checksums(archive, vital, expected={"EPIC_EMR.tar.gz": expected["EPIC_EMR.tar.gz"]})
        (vital / "lab_parameters.csv").write_text("changed synthetic dictionary\n")
        with self.assertRaises(ValueError):
            verify_source_checksums(archive, vital, expected=expected)

    def test_optimized_python_rejected_before_any_input_or_output_access(self):
        output = self.root / "should-not-exist"
        command = [sys.executable, "-B", str(REPO / "scripts/run_pipeline.py"), "run",
                   "--mover-archive", str(self.root / "absent.tar.gz"),
                   "--vitaldb-dir", str(self.root / "absent-csv"), "--output", str(output)]
        process = subprocess.run(command, env=dict(os.environ, PYTHONOPTIMIZE="1"),
                                 capture_output=True, text=True, check=False)
        self.assertEqual(process.returncode, 1)
        self.assertFalse(output.exists())
        self.assertEqual(json.loads(process.stdout)["error_type"], "RuntimeError")

    def test_existing_output_and_real_output_in_repository_rejected(self):
        with self.assertRaises(ValueError):
            initialize_output(self.root, demo=True)
        with self.assertRaises(ValueError):
            initialize_output(self.root / "real-output", demo=False)


if __name__ == "__main__":
    unittest.main(verbosity=2)

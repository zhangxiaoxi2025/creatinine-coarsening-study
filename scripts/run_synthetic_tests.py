"""Run all source-preserving synthetic suites without any real data input."""
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
suite = unittest.defaultTestLoader.discover(str(ROOT / "scripts"), pattern="test_*.py")
result = unittest.TextTestRunner(verbosity=2).run(suite)
report = {"scope": "SYNTHETIC_ONLY", "tests_run": result.testsRun,
          "failures": len(result.failures), "errors": len(result.errors),
          "skipped": len(result.skipped), "status": "PASS" if result.wasSuccessful() else "FAIL"}
(ROOT / "synthetic-test-report.json").write_text(json.dumps(report, indent=2) + "\n")
raise SystemExit(0 if result.wasSuccessful() else 1)

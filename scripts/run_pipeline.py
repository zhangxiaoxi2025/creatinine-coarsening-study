"""Explicit local source-to-results CLI; defaults reproduce the frozen design."""
import argparse
import json
from pathlib import Path
import sys

from portable_io import REPO, build_inputs, dump, initialize_output, sha, verify_source_checksums


def verify_cores():
    preserved = json.loads((REPO / "docs/core-preservation.json").read_text())
    for name, expected in preserved["unchanged_core_sha256"].items():
        if sha(REPO / "scripts" / name) != expected:
            raise ValueError("Preserved scientific core checksum changed")


def run(*, output, mover_archive=None, vitaldb_dir=None, demo=False, repetitions=5000):
    if not __debug__:
        raise RuntimeError("Optimized Python is forbidden: frozen assertion contracts must remain active")
    verify_cores()
    target = initialize_output(output, demo=demo)
    try:
        if demo:
            from synthetic_sources import create
            mover_archive, vitaldb_dir = create(target / ".private/synthetic_sources")
        elif repetitions != 5000:
            raise ValueError("Formal replication requires 5000 repetitions")
        if not demo:
            source_hashes = verify_source_checksums(mover_archive, vitaldb_dir)
            dump(target / "evidence/source-identity-verification.json", {"matches_study_snapshot": True, "sha256": source_hashes})
        build_inputs(mover_archive, vitaldb_dir, target)
        import run_frozen_analysis
        run_frozen_analysis.execute(target, repetitions=repetitions, enforce_expected=not demo)
        import run_extension
        run_extension.execute(target)
        checks = {"synthetic_only": bool(demo), "formal_resamples_per_mode": repetitions,
                  "centers": 2, "analyses_per_center": 8,
                  "status": "SYNTHETIC_END_TO_END_PASS" if demo else "COMPUTED_PENDING_INDEPENDENT_VERIFICATION"}
        if demo:
            from synthetic_sources import EXPECTED
            for center in ("mover", "vitaldb"):
                for analysis, expected in EXPECTED.items():
                    result = json.loads((target / "results" / f"{center}-{analysis}-analysis.json").read_text())
                    if result["N"] != expected:
                        raise ValueError("Synthetic expected denominator mismatch")
            checks["synthetic_expected_denominators"] = EXPECTED
        dump(target / "evidence/pipeline-complete.json", checks)
        return checks
    except Exception as exc:
        # Data exceptions may contain raw values or identifiers: never echo them.
        dump(target / "evidence/pipeline-failure.json", {"status": "STOPPED",
            "error_type": type(exc).__name__, "raw_error_text_omitted": True})
        raise RuntimeError("Pipeline stopped. Review local evidence; no raw error text is emitted.") from None


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    formal = commands.add_parser("run", help="Run on independently obtained authorized study sources")
    formal.add_argument("--mover-archive", type=Path, required=True)
    formal.add_argument("--vitaldb-dir", type=Path, required=True)
    formal.add_argument("--output", type=Path, required=True)
    demo = commands.add_parser("demo", help="Generate and run fictional source-format fixtures")
    demo.add_argument("--output", type=Path, required=True)
    demo.add_argument("--repetitions", type=int, default=20)
    args = parser.parse_args(argv)
    if args.command == "demo" and not 1 <= args.repetitions <= 5000:
        parser.error("Demo repetitions must be between 1 and 5000")
    try:
        report = run(output=args.output, demo=args.command == "demo",
                     repetitions=getattr(args, "repetitions", 5000),
                     mover_archive=getattr(args, "mover_archive", None),
                     vitaldb_dir=getattr(args, "vitaldb_dir", None))
    except Exception as exc:
        print(json.dumps({"status": "STOPPED", "error_type": type(exc).__name__,
                          "detail": "Check schemas, permissions, frozen denominator contract, and local evidence."}))
        return 1
    print(json.dumps(report))
    return 0


if __name__ == "__main__":
    sys.exit(main())

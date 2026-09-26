# Numerical coarsening and creatinine-based AKI research labels

Version 1.0.0. Research software accompanying **Effects of numerical coarsening on creatinine-based acute kidney injury research labels: a paired empirical study**.

This repository includes the complete source-to-results pipeline: source archive/CSV ingestion, frozen cohort construction, exact observation windows and numerical transformations, all 16 center-by-analysis resampling experiments, and the post-result exploratory E1–E3 supplement. It contains no patient data, credentials, access approvals, or data-use agreements.

The original study execution underwent independent verification. The portable source I/O and orchestration released here have been validated with synthetic source-format fixtures only; they have not been rerun on the original patient files as part of packaging. Identical numerical cores and unchanged scientific cohort-function ASTs are documented in `docs/core-preservation.json`. A successful synthetic demonstration is not a new validation of clinical accuracy or of real-data replication.

## Install and test

Use Python 3.11 or newer without `-O` or `PYTHONOPTIMIZE` (the CLI rejects optimized execution before opening inputs) in an environment containing NumPy, pandas, and PyArrow. Tested versions are listed in `requirements.txt`.

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python -B scripts/run_synthetic_tests.py
```

The test command runs the original scalar, independent oracle, vector, adapter, RNG, and exploratory suites, plus synthetic portable-I/O checks. It writes `synthetic-test-report.json`. Tests do not download or discover study data.

## Complete synthetic demonstration

Choose a fresh output directory:

```sh
python -B scripts/run_pipeline.py demo --output ../aki-synthetic-demo --repetitions 20
```

The demo creates explicitly fictional MOVER-format tar members and VitalDB CSVs; reads them through the same ingestion functions; constructs both cohorts; runs all eight analyses for each center; and produces E1, E2, and E3. Cases exercise repeated operations, diagnoses, conflicting timestamps, duplicate records, extreme values, long anesthesia, early targets, and E7-only observation. Hand-declared synthetic denominators are checked. The portable end-to-end test checks all 16 denominators, repeated-run numeric agreement and E1/E2 internal count identities; it is not an independent expected-results oracle for every generated statistic. Reduced resampling is available only in the demo; it must not be presented as the study's inferential run.

`python -B scripts/run_frozen_analysis.py` accepts the same CLI as `run_pipeline.py`.

## Reproduce from independently obtained source files

Obtain the study-compatible MOVER EPIC source archive and VitalDB clinical/laboratory release through their data custodians, and meet all applicable access and reuse conditions independently. The software neither grants data access nor downloads files.

```sh
python -B scripts/run_pipeline.py run \
  --mover-archive /authorized-data/EPIC_EMR.tar.gz \
  --vitaldb-dir /authorized-data/vitaldb-clinical-labs \
  --output /private-results/aki-replication
```

The MOVER archive must contain exactly one regular member for each basename `patient_information.csv`, `patient_labs.csv`, and `patient_visit.csv`; parent archive directories are allowed. Members are read with `tar.extractfile`, never extracted into the filesystem. Only LOINC `2160-0` with case-insensitive exact units `mg/dl` enters the MOVER creatinine table. VitalDB requires `clinical_data.csv`, `lab_data.csv`, `clinical_parameters.csv`, and `lab_parameters.csv`. The parameter dictionaries are included in source identity verification; clinical and laboratory tables supply the numerical records; only `name == cr` enters its creatinine table. CSV strings are preserved. Required headers, including duplicate-header rejection, are checked before parsing.

The study used MOVER's EPIC snapshot dated 2026-07-20 and VitalDB clinical/labs release 1.0.0. MOVER's `BIRTH_DATE` field in that snapshot is interpreted as the supplied age encoding, not parsed as a calendar birth date. Other releases or schemas are not automatically supported. Formal execution compares all five input-file SHA-256 values with `config/source-checksums.json` before source ingestion. This manifest contains only basenames and hashes transcribed from existing verified study-snapshot checksum records, without reopening patient source files during packaging. Cohort denominator agreement alone does not authenticate source identity or establish access rights. The repository does not contain protected source paths or access receipts. Operators must verify that their legally obtained source versions match the study.

Formal `run` always uses 5,000 resamples per center/analysis, the frozen seed 20260921, the two-level PCG64 stream hierarchy, both centers, and all eight analyses. It checks all 16 frozen cohort denominators. A denominator or numeric contract mismatch stops execution; there is no switch to silently accept a changed cohort. Do not edit rules to force a matching result.

## Outputs and handling

Output directories must be new. Formal output must lie outside this repository. The output root and `.private/` are created with mode 0700; generated patient-level files use 0600. These POSIX permissions are defense in depth; users remain responsible for the security of their environment.

- `results/`: cohort-flow summaries, 16 analysis JSON files, all paired comparison cells and intervals, mapping lookups, support/changed-record diagnostics, observation summaries, two extension reports, patient-equal mapping lookups, and resampling diagnostics.
- `.private/`: extracted source tables, frozen patient inputs, patient labels, extension states, and bootstrap count arrays. Keep this directory local and access-controlled.
- `evidence/`: source/code hashes, execution metadata, and completion or sanitized failure reports. Raw exception text is omitted because it may contain source values or identifiers.

The code does not automatically publish even aggregate outputs. Review any intended release under the applicable data terms, including mappings with source support values. Do not upload a run output directory into this repository. `.gitignore` excludes known run artifacts, but it is not an approval mechanism.

A formal completion status is `COMPUTED_PENDING_INDEPENDENT_VERIFICATION`: the portable run is not automatically promoted to the original study's independently verified status.

## Scientific scope

This is a paired numerical-representation experiment. Original released values are the computational reference, not a clinical gold standard. E48 and E7 are observed eligible target sets, not complete surveillance. The algorithm excludes urine-output criteria and does not implement complete clinical AKI diagnosis, staging, or treatment recommendations. **Do not use this research software for clinical decisions.**

Read `METHODS.md` for frozen rules and interpretation; `CHANGELOG.md` for the pre-result implementation changes and post-result exploratory boundary. The preserved `analysis-spec-v1.0.json` is the historical design specification; its cent-path wording and statement that patient-equal analyses were outside v1.0 are superseded only by the dated changes described there. The production adapter uses exact 10^-15 mg/dL units, while the legacy cent kernel is retained for equivalence testing.

Configuration contains 2,010 procedure and 1,623 diagnosis terminology decisions with source frequencies and source row indices removed. Procedure terms not in the frozen map stop the run. Diagnosis screening applies the documented exact code/term exclusions; unmapped diagnoses are not assumed to establish renal health. These maps are study operational rules, not adjudicated clinical diagnoses.

## Figures

The four CSVs in `aggregate-figures/` are frozen aggregate results used by the submission figures, separate from the fictional demo outputs. To rebuild the figures without patient records:

```sh
Rscript scripts/plot_public_figures.R aggregate-figures ../aki-public-figures
```

The script uses R `grid`/`grDevices` and requires Cairo PDF and PNG support. It was tested with R 4.6.0 on macOS. It writes selectable-text PDF and 300 dpi PNG. A 600 dpi LZW TIFF is produced with native TIFF support, or on macOS through Quartz plus system `sips`; otherwise this optional format is omitted. SVG is optional when `svglite` is installed. All four PNG/TIFF outputs were checked against the submission figures. Center-level counts and estimates are unchanged.

## License and citation

Software is provided under the MIT license in `LICENSE`. Source datasets have their own terms and are not covered by that software license. `CITATION.md` gives the software title and version without inventing author names or a repository DOI. Project home: https://github.com/zhangxiaoxi2025/creatinine-coarsening-study . The versioned release is https://github.com/zhangxiaoxi2025/creatinine-coarsening-study/releases/tag/v1.0.0 .

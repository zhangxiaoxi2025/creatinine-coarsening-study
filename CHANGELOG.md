# Change history

## 1.0.0 — 2026-09-25, portable code package

Added explicit local source/output CLI, MOVER archive member ingestion, VitalDB CSV validation, portable output handling, synthetic source fixtures, end-to-end tests, documentation, and terminology-only algorithm configurations. No real patient file was read or rerun while building this package. Five numerical/adapter files remain byte-identical to the analyzed implementations; cohort construction scientific functions retain their original ASTs. Source parsing/orchestration changes are synthetic-tested portability changes, not a newly independently verified real-data execution.

## Exploratory supplement — 2026-09-24

After the original results were known, fixed E1 criterion-state decomposition, E2 original-map coverage strata, and one E3 patient-equal reference distribution were locked before extension output. These are post-result exploratory analyses. No additional resampling intervals or model search were added. Original main outputs and the frozen original plan were retained.

## Exact numerical implementation amendment 1.0.1 — 2026-09-23

Before real outcome computation, source precision checks identified legitimate literals beyond two decimal places. The proposed cent-only production path was replaced with exact 10^-15 mg/dL integer units for both centers and every analysis. Decimal/Fraction conversion must be exact; no rounding or outcome-driven exclusion is permitted. Population, thresholds, observation indices, transformation definitions, and resampling design did not change. The original cent kernel remains here only for historical equivalence tests.

## Frozen design/RNG implementation review — 2026-09-23

The original design was proposed on 2026-09-21 and finalized on 2026-09-23 before outcomes. A pre-result review repaired the RNG interface to implement the specified two-level center/analysis SeedSequence hierarchy. `study_rngs()` and explicit PCG64 checks replaced reliance on a generic center-only helper for study execution. CSV BOM handling in the independent input verifier was also repaired during pre-result validation; this did not change the scientific rules. These implementation repairs were documented rather than relabeled as new prespecified scientific analyses.

The retained `analysis-spec-v1.0.json` records the historical frozen design. Later exact-unit and exploratory changes are described explicitly here and in `METHODS.md`.

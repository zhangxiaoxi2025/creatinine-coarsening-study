# Scientific rules and reproducibility boundary

## Frozen population

Use one operation per patient in the complete source, adults receiving general anesthesia, valid anesthesia and admission/discharge ordering, and the frozen noncardiac procedure map. Conflicting MOVER operation identities or critical metadata exclude the operation. MOVER laboratory and same-visit diagnosis linkage requires both case and patient identifiers to match. VitalDB case identifiers must be unique. The full procedure vocabulary is retained as algorithm configuration, with counts and source row indices removed.

Keep only procedure decisions `include`; `exclude` and `uncertain` are excluded. Same-visit documented ESRD, chronic dialysis dependence, renal transplant recipient state, and uncertain RRT state are screened using the exact frozen diagnosis mappings. Missing or unmatched diagnosis evidence is not proof of absence of disease. Ordinary CKD and recorded AKI are not blanket exclusions. The primary main denominators are MOVER 4,321 and VitalDB 2,059; the E7 parallel denominators are 4,433 and 2,540. All sensitivity denominators are in `config/expected-denominators.json`.

## Values and fixed observations

MOVER uses LOINC 2160-0 in mg/dL; VitalDB uses the released `cr` laboratory variable. Source numerical strings are preserved. The raw inclusion screen is finite 0 < creatinine < 100 mg/dL with a valid timestamp; this is not a physiological-validity guarantee. The global reference window is [anesthesia start − 604800 seconds, min(anesthesia end + 604800 seconds, discharge)]. Exact duplicate patient/time/value rows are collapsed; entire timestamps containing conflicting values are removed before observation sets are constructed.

Require at least one reference record in [anesthesia start − 7 days, anesthesia start). Postoperative targets are in (anesthesia end, min(anesthesia end + 7 days, discharge)]. A strict earlier 48-hour or 7-day reference is required for the respective E48 or E7 target, with inclusive left boundaries and exclusive current time. Minima are calculated within those fixed prior windows. Reference records outside E48, including early postoperative observations, can still supply later priors.

All representation versions use the same patient, timestamp, target, and prior-index sets. At eligible E48 targets the main rule is current − min48 >= 0.3 mg/dL OR 2 × current >= 3 × min7. The patient label is OR over eligible targets. E7 parallel analysis uses only the ratio rule at E7 targets, with its own population and fitted map. Missing eligibility is not interpreted as a negative label.

## Exact representation and four transforms

Source text is converted using Decimal/Fraction to integer units of 10^-15 mg/dL. The conversion must be exact; values requiring finer precision stop execution. Integer arithmetic and explicit bounds protect threshold calculations. Values are never rounded during input conversion.

Fit the reference distribution directly to all retained global-window observations from the selected population, including preoperative, intraoperative, and eligible/ineligible postoperative records. Original main fitting is measurement-equal, using patient multiplicities during bootstrap. Type-1 inverse-CDF quantiles are used.

- I: original released values.
- T: clip to Q(0.025) and Q(0.975).
- B: lower and upper tail values map to the respective endpoints; other values use their original-distribution midrank to select the nearest node at 0.05, 0.10, ..., 0.95 and output that node's quantile. Equal-distance ties choose the lower node.
- G: max(0.2, 0.2 × floor(x / 0.2 + 0.5)). This is half-up rounding with an explicit positive floor.

No post-transform quality screen or cohort reselection is allowed, even when a transformed value equals 100 mg/dL. Monotonicity permits fixed original minima to be transformed directly.

## Paired outcomes and resampling

For each center report n00, n01, n10, n11; the first label is I and the second is the transformed version. D = (n01 + n10) / N. Directional rates are n01 / N and n10 / N; net change is (n01 − n10) / N. Opposite changes can cancel in the net difference.

Run main, parallel7, S1_no_diagnosis, S2_target_after5min, S3_exclude_extremes, S4_exclude_over24h, S5_exclude_conflicts, and S6_all_postop_E48. Each mode refits T/B for its selected cohort. S2 changes targets only, retaining early records as possible priors. S3 excludes whole patients with any retained raw value <0.1 or >30; S4 excludes durations >24 hours; S5 excludes patients with an original within-window conflicting timestamp. S6 retains patients whose retained observed postoperative points all belong to E48 and refits its maps.

There are 16 experiments: 2 centers × 8 analyses × 5,000 whole-patient resamples, totaling 80,000 resamples. One SeedSequence(20260921) spawns center streams in MOVER, VitalDB order; each center spawns eight streams in the fixed analysis order, with PCG64 generators. Complete reference and target clusters receive their patient's sampled multiplicity. Four representations share each replicate's patient weights; T/B are refitted in every replicate. Patients are ordered lexicographically by source case-id strings.

Report empirical type-1 2.5th and 97.5th percentile ranges as patient-resampling stability ranges. Population coverage for this discrete fitted-threshold statistic is not established; these are not automatic 95% confidence guarantees or equivalence tests. A zero central range is distinct from every replicate being zero. Distinct D values describe the latter property directly. Fitted map/quantile-combination counts are center-by-analysis diagnostics shared across representation rows; identity and fixed-grid representations are not being refitted. Quantile-combination counts include both tail endpoints plus the 19 interior nodes. Do not sum repeated diagnostics over the four representations.

## Post-result exploratory E1–E3

The extension was specified on 2026-09-24 after original results were known and before extension calculations. It retains the original main patients, global references, and E48 indices.

E1 reports a 4×4 transition table with state 2*A + R: neither, ratio only, absolute only, both. Each criterion is an OR over all eligible E48 times; both need not occur at the same time. Folding states 1/2/3 positive must recover the original main paired table. Criterion-specific changes overlap and cannot be added as disjoint patient counts.

E2 partitions the original main population by whether all retained observed postoperative points belong to E48, using the original full-cohort B mapping in both strata. All counts must add back to the main table. This differs from S6, which changes the population and refits its map. Coverage of observed records does not reconstruct unobserved measurements.

E3 fits one patient-equal alternative: every retained reference record from patient i gets exact Fraction weight 1/m_i; each patient contributes total weight 1. Normalize by N for inverse-CDF quantiles and midranks. The same tail probabilities, nodes, lower-tie rule, timestamps, thresholds, and patient OR remain fixed. This retains each record rather than averaging a patient's trajectory. Only descriptive point estimates are computed; original B intervals must not be copied to the new mapping. No alternative weights or cutoffs are searched.

## Portable validation

The original study analysis and extension were independently verified in their original controlled execution. Portable packaging preserves five numerical/adapter files byte-for-byte and the scientific cohort functions by AST, while replacing environment-specific access gates, paths, and command-line orchestration. Original synthetic test functions and classes are preserved apart from filesystem lookup/output changes.

The distributed source I/O is tested only on synthetic archives/CSVs: schema and unit filtering, malformed/duplicate member rejection, both full center pipelines, all denominator contracts, repeated-run numerical consistency, E1–E3 folding, private output permissions, and absence of patient identifiers in aggregate outputs. The real source archive has not been loaded in this packaging task. The formal CLI requires the five source files to match `config/source-checksums.json`, which was transcribed from existing verified source-snapshot checksum records. Matching bytes does not establish permission. A new operator must independently verify their access rights and generated results.

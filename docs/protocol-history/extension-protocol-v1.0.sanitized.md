> Sanitized historical copy prepared on 2026-09-29; distributed with software v1.0.1. The versioned release record supplies the actual publication timestamp. Original protocol date: 2026-09-24. The main results were known when this exploratory protocol was written; its internal lock preceded extension computation. Internal workflow and authorization wording has been removed or generalized, and computational-verification terminology clarified. Scientific definitions, expected count checks and stopping conditions are retained. Original and edited hashes are recorded separately in provenance.json.

# Post-result exploratory supplement protocol v1.0

Status: defined and locked before new extension results; original results are already known. Date: 2026-09-24. This document does not amend the frozen original plan retroactively.

## Population and values
Use the original MOVER 4321 and VitalDB 2059 primary patients, exact original retained records, global reference windows, E48 target times and 48-hour/seven-day prior indices. Reconstruct with the unchanged validated main adapter. Retain the pre-outcome exact 10^-15 mg/dL representation. Original main outputs and manuscript v1.0 remain immutable. No alternative cohort, thresholds, nodes, subgroup cutoffs or models will be searched.

## E1 Criterion state decomposition
For I and the fixed original main B mapping, compute A_f as any E48 target satisfying the absolute0.3mg/dL rule and R_f as any E48 target satisfying the1.5 ratio rule. State is 2*A_f+R_f:0 neither,1 ratio only,2 absolute only,3 both across any eligible times. Both need not occur at the same time. Produce one4x4 patient transition table per center and derived directional summaries. Folding states1/2/3 positive must exactly recover original n00,n01,n10,n11 and181/73 discordant patients. Do not add individual criterion-flip counts as if disjoint or make causal/clinical-error claims.

## E2 Fixed-mapping observation coverage strata
Partition each main cohort by whether all retained observed postoperative points belong to E48. Keep the original main B mapping fixed for both strata. Report N, four cells, D, upward/downward counts, original/transformed positive counts and net difference. Each cell and denominator must sum exactly to the main table. No p-values, between-stratum hypothesis tests or new bootstrap. Existing S6 is a separate changed-population/refitted-map analysis; report the distinction explicitly.

## E3 One patient-equal mapping
Each retained reference observation from patient i receives exact rational weight1/m_i, where m_i is that patient's number of retained global-window reference records. Normalize total weight by N. Fit weighted inverse-CDF quantiles Q(p)=min{x:F_w(x)>=p}, weighted midrank M(x)=(weight below x+0.5*weight equal x)/N, tail endpoints Q(.025)/Q(.975), same interior nodes.05,.10,...,.95, nearest-node ties lower, same direct fit to original values. Retain original timestamps, per-patient OR label and all scientific thresholds. Use exact Fraction arithmetic; no floating-point boundary decisions. Do not summarize a patient's trajectory into one value.

This is a descriptive full-cohort point-estimate extension. Report original-versus-new-B2x2, D and directions, change inD relative to the original fit, and mapping endpoints/support. No additional resampling interval is computed; do not attach originalB intervals to the new mapping. This decision is locked before extension outputs and avoids implying unassessed inferential coverage. Retain this one alternative regardless of result direction.

## E4 Existing aggregate diagnostics and writing
Reformat original mappings' support sizes and changed-record counts, extra grid-floor counts, E7 eligibility loss8/2 and existing resampling diagnostics. Do not rerun original80,000 replicates. Report D percentile-interval degeneracy separately from direction-rate degeneracy; distinguish all sampled values equal from a central range with equal endpoints wherever existing aggregates support the distinction, otherwise mark not assessed. No new average-distortion analysis is included in this bounded extension.

## Verification and failure rules
Source/version and documented permitted use must be verified before any patient-data reuse. Tests must exercise exact weighted-CDF ties, weighted midrank ties, unequal measurement counts, constant support, permutation invariance, threshold equality and OR-state collapse.

A separately implemented computational check must reconstruct timing and reference minima directly from original frozen retained records, compute I/B/B_patient states separately, and compare every patient state, fixed-coverage stratum, mapping support/output and all aggregate tables. Reuse original saved B lookup as fixed-map evidence. Verify new weighted map independently with rational weights. Reports expose aggregate results only; patient-level verification vectors remain access-controlled. This computational check does not establish a second human analyst’s complete independent rerun.

STOP on source/version/access mismatch, patient/time/unit ambiguity, differing patient or target sets, non-deterministic or non-monotone mappings, overflow/exactness failures, inability to fold E1 or sum E2 to the original table, or independent disagreement. Preserve failed evidence. Do not change methods because a result is inconvenient.

## Output and interpretation
Produce local exploratory results and verification report; revised English manuscript and supplement v1.1; Chinese change guide. Clearly distinguish original prespecified results from post-result exploratory extensions. This was a local post-result exploratory work package; it did not publicly preregister the original or exploratory analysis.

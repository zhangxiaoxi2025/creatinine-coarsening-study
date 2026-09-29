> Sanitized historical copy prepared on 2026-09-29; distributed with software v1.0.1. The versioned release record supplies the actual publication timestamp. Original amendment date: 2026-09-23. Internal authorization/recovery accounting was omitted and computational-verification terminology clarified. Scientific numerical requirements and aggregate precision-audit counts are retained. This implementation amendment number 1.0.1 is distinct from the later software documentation release 1.0.1. Original and edited hashes are recorded separately in provenance.json.

# Implementation amendment v1.0.1 — exact representation before outcome computation

Date: 2026-09-23. Frozen scientific plan remains analysis-plan-v1.0.md with unchanged population, reference pool, windows, thresholds, transforms, primary comparison, sensitivities and RNG hierarchy.

## Trigger and outcome awareness
A separately implemented computational check stopped BEFORE producing any real patient AKI labels, mapping results or effect estimates: legitimate raw MOVER strings in selected populations are not all exact multiples of 0.01 mg/dL. Aggregate audit: main 13 records/9 patients; S1 23 records/14 patients; normalized fractional precision 15 places. VitalDB selected modes contain no such records. No result-driven exclusion or numeric rounding is permitted.

## Exact implementation repair
It preserves every original numeric literal and the scientific calculation. This is an explicit implementation amendment to the proposed cent-only vector path, not a retrospectively prespecified new analysis.

Use integer units of 10^-15 mg/dL (unit_scale=10^15) for ALL selected observations in both centers and all modes. Require exact conversion using Decimal/Fraction: value*scale must be an integer with no rounding, finite and 0<raw<100. Raw maximum implies units<10^17, leaving room for factors2/3 and rounding increments in signed int64; actual source maxima and exact-conversion checks must also pass. If any selected value is not exactly representable or arithmetic bounds fail, STOP again.

Absolute threshold is 3*scale/10; ratio is 2*current>=3*prior7. Positive-floor fixed-grid step is 2*scale/10, half-step=scale/10. Other fitted quantile/midrank calculations use unchanged integer observation counts and output the corresponding exact source units. No raw-value normalization, removal or new post-transform screening occurs.

The frozen Fraction/Decimal scalar definition already supports these values. A separate scripts/paired_resampling_exact.py implements the generalized integer path; D remains unchanged. Synthetic checks compare scale100 against the old kernel and scale10^15 against the scalar definition, including near-threshold1e-15 values, grid ties, zero-floor and >99.99 transformed outputs. The runner/adapter may proceed only after these checks pass. Comparisons against a separately implemented scalar calculation for patient labels, tables, first bootstrap samples and range extraction remain required. These are computational checks within the same AI-assisted workflow, not a second human analyst’s complete real-data rerun.

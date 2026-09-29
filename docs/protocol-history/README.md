# Internal protocol history and sanitized publication copies

Prepared: 2026-09-29. These records are distributed with software v1.0.1. The versioned release record supplies the actual public timestamp; the preparation date is not a public preregistration date.

The design was drafted on 2026-09-21 and internally frozen on 2026-09-23 before label/effect calculations. The numerical representation amendment was recorded before those calculations on 2026-09-23. The exploratory protocol was internally locked on 2026-09-24 after the original results were known and before the extension calculations. This sequence is supported by local records. It was not publicly preregistered or independently timestamped.

The first software release, [v1.0.0](https://github.com/zhangxiaoxi2025/creatinine-coarsening-study/releases/tag/v1.0.0), was published on 2026-09-26 after the main analysis. It contained the historical JSON analysis specification and methods/change history, but not these original protocol copies or freeze/lock summaries. The package-preparation date 2026-09-25 must not be used as its public-release date. Publishing these documents later improves traceability; it does not retrospectively create preregistration.

## Files

- `analysis-plan-english-summary.md`: retrospective English summary prepared on 2026-09-29, after results; it is not an original pre-analysis registration.

- `analysis-plan-v1.0.sanitized.zh.md`: scientific content of the original Chinese design plan, with internal administrative and access details removed.
- `exact-implementation-amendment-20260923.sanitized.md`: lossless numerical-precision amendment; this historical amendment number is independent of the later software release numbering.
- `extension-protocol-v1.0.sanitized.md`: post-result exploratory E1–E3 protocol and existing-diagnostics reporting plan.
- `freeze-record.sanitized-summary.json`: internal freeze timestamp and hashes of the original plan and specification.
- `extension-lock.sanitized-summary.json`: internal exploratory lock timestamp, result-awareness state and original protocol hash.
- `provenance.json`: original-document hashes and separate edited-file hashes, with editing descriptions.
- `SHA256SUMS`: hashes of all other files in this directory.

The root `analysis-spec-v1.0.json` remains byte-identical to the historical specification. Its proposed cent-only vector path was superseded by the dated precision amendment; its exclusion of patient-equal analyses from the original plan is compatible with the later expressly exploratory extension. Historical references to a stability interval describe the same empirical percentile calculation now reported as a **resampling range**, not a calibrated population confidence interval.

## Limits of this evidence

File hashes establish document identity and allow checks for later modification. Internal timestamps and copied hashes do not independently establish when an external party could see a document. The supplied public copies are edited for disclosure safety; their hashes differ from the original source hashes. Sensitive authorization records, machine-specific paths, credentials, individual records and patient-derived artifact fingerprints are not included in these new history files.

References to separate or independent implementations describe computational checks performed within the same AI-assisted analytical workflow. They do not establish that a second human analyst independently rebuilt the complete real-data pipeline. Current author contribution statements, where supplied with the manuscript, are separate from this computational evidence.

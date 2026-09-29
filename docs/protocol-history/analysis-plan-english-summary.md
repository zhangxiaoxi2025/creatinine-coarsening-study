# English summary of the internally frozen analysis plan

Prepared for documentation on 2026-09-29, after the study results were known. The original Chinese plan was drafted on 2026-09-21 and finalized and internally frozen on 2026-09-23, before real-patient labels, fitted cutpoints, transformation effects or resampling ranges were calculated. Those historical dates are unchanged. Original and screened-copy hashes are recorded separately in provenance.json. This English summary is retrospective documentation, not a publicly preregistered protocol.

## Objective and data

The plan defined a paired empirical methods study of creatinine-based research labels in MOVER EPIC and VitalDB clinical/laboratory release 1.0.0. It compared original released values with predefined transformations while fixing patients, retained sampling times and temporal rules. Original values were the computational reference. The primary estimand was within-cohort patient-level label discordance, reported separately for each center without a pooled effect. Prediction, treatment effects and validation of complete clinical AKI diagnoses were outside scope.

## Population and records

Operation counts were established from the complete released clinical source before downstream filtering. Eligible adults had one operation in that source, general anesthesia, valid chronological anesthesia/admission/discharge times, and an included noncardiac procedure. MOVER cases with conflicting critical identity or clinical fields were excluded; laboratory linkage required agreement of case and patient identifiers. VitalDB used unique case identifiers and patient identifiers. Top-coded ages established eligibility without inventing exact ages.

The frozen exact-name procedure map retained included terms and excluded both excluded and uncertain terms, including explicit kidney-transplant or dialysis procedures. Recorded end-stage kidney disease, chronic dialysis dependence, kidney-transplant recipient status/complications, and renal replacement states of uncertain acute/chronic timing triggered diagnosis-based exclusion. Ordinary CKD, CKD stage 5, AKI, ambiguous CKD5-or-ESRD descriptions and donor records did not automatically exclude patients. MOVER diagnosis linkage stayed within the same operation and patient. These restrictions did not establish complete preoperative renal histories.

MOVER used LOINC 2160-0 measurements in mg/dL; VitalDB used its official creatinine variable. Valid original values were finite and strictly between 0 and 100 mg/dL. The global window ran from seven days before anesthesia start through the earlier of discharge and seven days after anesthesia end. Identical patient–time–value records were deduplicated; every record at a timestamp containing conflicting values was removed. At least one preoperative measurement was required.

## Time sets and labels

Postoperative targets were strictly after anesthesia end. Prior windows were [t−48 hours, t) and [t−7 days, t), intersected with the global window; each used its lowest retained value. E48 contained targets with a strictly earlier value within 48 hours. A nonempty E48 was required for the main comparison. A patient was positive if any E48 target satisfied an absolute rise of at least 0.3 mg/dL over the 48-hour minimum or a ratio of at least 1.5 to the seven-day minimum. Missing target eligibility was unobserved, not negative. The parallel E7 analysis required a seven-day prior value and used only the ratio criterion, with a separately selected cohort and refitted mapping.

## Transformations and weights

Each center–analysis combination used all retained global-window measurements from eligible patients as its reference distribution, including intraoperative records and records outside the target set. Records had equal weight; frequent testing therefore contributed more reference mass. Patient-equal quantiles were explicitly not performed in the original plan.

Let Q denote the type-1 empirical inverse-distribution quantile and M(x) the fraction below x plus half the fraction equal to x. I retained original values. T clipped values to Q(0.025) and Q(0.975). The primary transformation B mapped values at or beyond those endpoints to the respective endpoint. Interior values were assigned to the nearest M(x) node among 0.05, 0.10, …, 0.95, with exact ties assigned downward, and replaced by Q at that node. Coincident endpoints mapped everything to that value. G was max(0.2, 0.2 × floor(x/0.2 + 0.5)), in mg/dL. T and B were fitted directly from original values. No post-transformation reselection, matched privacy guarantee or equal output-support constraint was specified.

## Estimates, resampling and sensitivities

For I versus B, the paired table contained n00, n01, n10 and n11, with original status first. The primary statistic was D=(n01+n10)/N. Directional changes, positive-label proportions and net change (n01−n10)/N were explanatory summaries. I versus T, I versus G and E7 results were secondary.

Each combination used 5,000 whole-patient resamples with shared patient multiplicities across transformations and refitted T/B mappings. PCG64 used SeedSequence(20260921), split by MOVER/VitalDB, then main, parallel7 and S1–S6. Empirical type-1 2.5th and 97.5th percentiles summarized resampling variation; population coverage was unestablished. Current reporting calls these resampling ranges.

Six separate sensitivities retained other rules: omit diagnosis exclusions; require targets more than five minutes after anesthesia end while retaining earlier records as potential priors; exclude patients with retained values below 0.1 or above 30 mg/dL; exclude anesthesia durations above 24 hours; exclude patients with in-window conflicting timestamps; and retain patients whose every observed postoperative point belonged to E48. Each refitted its reference mapping. No result-driven combinations were specified.

## Implementation history and stopping rules

The original plan permitted exact decimal/rational calculations and a vectorized integer-cent implementation only when every value multiplied by 100 was integral; otherwise execution had to stop. As documented in Supplement S6, this restriction triggered a stop before real label/effect calculations. A subsequent implementation amendment retained scientific rules but generalized lossless integer representation to 10^-15 mg/dL. This representation unit is not assay resolution.

Execution also stopped for unresolved identity/unit/time issues, empty cohorts, inconsistent implementations, invalid transformations or changing paired membership/times. Failed resamples could not be discarded and replaced. Effect size, significance or differing center results could not justify changing rules. Later changes required a documented amendment and disclosure of whether results were known.

The extension dated 2026-09-24 added criterion-state tables, fixed-mapping coverage strata and one patient-equal mapping after the original results. These remained exploratory and did not retroactively become prespecified analyses.

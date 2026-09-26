"""E1-E3 post-result exploratory supplement on the fixed main cohorts."""
from fractions import Fraction
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis_adapter import build_analysis
from extension_core import criterion_states, four_by_four, paired_table, patient_equal_mapping
from portable_io import dump, sha


def execute(output):
    output = Path(output)
    results = []
    for center in ("mover", "vitaldb"):
        cohort = pd.read_parquet(output / ".private" / f"{center}-cohort.parquet")
        labs = pd.read_parquet(output / ".private" / f"{center}-measurements.parquet")
        adapted = build_analysis(cohort, labs, analysis="main")
        design, keys = adapted.design, adapted.patient_keys
        fit = design.fit_mappings()
        original = criterion_states(design, design.support, design.support)
        fixed = criterion_states(design, design.support, fit.mappings["tail_plus_bins"])
        main = json.loads((output / "results" / f"{center}-main-analysis.json").read_text())
        fixed_table = paired_table(original, fixed)
        for key, value in main["point_estimates"]["tail_plus_bins"].items():
            if fixed_table[key] != value:
                raise ValueError("E1 failed to recover the original main comparison")
        matrix = four_by_four(original, fixed)
        folded = [int(matrix[0, 0]), int(matrix[0, 1:].sum()),
                  int(matrix[1:, 0].sum()), int(matrix[1:, 1:].sum())]
        if folded != [fixed_table[k] for k in ("n00", "n01", "n10", "n11")]:
            raise ValueError("Criterion-state OR collapse failed")
        selected = labs[labs.case_id.isin(keys)]
        incomplete = set(selected.loc[selected.is_post & ~selected.E48, "case_id"])
        complete = np.asarray([key not in incomplete for key in keys], dtype=bool)
        strata = {"all_observed_postop_E48": paired_table(original, fixed, complete),
                  "some_observed_postop_outside_E48": paired_table(original, fixed, ~complete)}
        for key in ("n", "n00", "n01", "n10", "n11", "discordance_count"):
            if sum(table[key] for table in strata.values()) != fixed_table[key]:
                raise ValueError("Fixed-map coverage strata do not recover the main comparison")
        weighted = patient_equal_mapping(design.reference_patient, design.reference_units,
                                         design.n_patients)
        patient_equal = criterion_states(design, weighted["support"], weighted["mapped_units"])
        table = paired_table(original, patient_equal)
        lookup = []
        for raw, transformed, weight, midrank in zip(weighted["support"], weighted["mapped_units"],
                                                      weighted["support_weights"], weighted["weighted_midranks"]):
            lookup.append({"raw_units": int(raw), "B_patient": int(transformed),
                           "unit_scale": design.unit_scale,
                           "weight_numerator": weight.numerator, "weight_denominator": weight.denominator,
                           "midrank_numerator": midrank.numerator, "midrank_denominator": midrank.denominator})
        pd.DataFrame(lookup).to_csv(output / "results" / f"{center}-patient-equal-mapping.csv", index=False)
        path = output / ".private" / f"{center}-extension-labels.parquet"
        pd.DataFrame({"case_id": keys, "I_state": original, "B_state": fixed,
                      "B_patient_state": patient_equal, "complete_E48": complete}).to_parquet(path, index=False)
        path.chmod(0o600)
        result = {"center": center, "N": design.n_patients,
                  "analysis_status": "post-result exploratory descriptive point estimates",
                  "criterion_state_order": ["neither", "ratio_only", "absolute_only", "both"],
                  "criterion_transition_4x4": matrix.tolist(), "fixed_B_table": fixed_table,
                  "fixed_mapping_strata": strata, "patient_equal_table": table,
                  "patient_equal_delta_discordance": table["discordance"] - fixed_table["discordance"],
                  "same_composite_label_but_different_criterion_state": int(np.sum((original != fixed) & ((original > 0) == (fixed > 0)))),
                  "patient_equal_mapping": {"q025_units": weighted["q025"], "q975_units": weighted["q975"],
                      "node_values": weighted["node_values"], "support_size": len(weighted["support"]),
                      "output_support_size": len(np.unique(weighted["mapped_units"])),
                      "weight_total": str(sum(weighted["support_weights"], Fraction(0)))},
                  "no_new_interval_computed": True, "unit_scale": design.unit_scale,
                  "private_output_sha256": sha(path)}
        dump(output / "results" / f"{center}-extension.json", result)
        results.append(result)
    diagnostics = []
    for path in sorted((output / "results").glob("*-analysis.json")):
        analysis = json.loads(path.read_text())
        diagnostic = analysis["diagnostics"]
        for version, summary in diagnostic["versions"].items():
            intervals = analysis["intervals"][version]
            diagnostics.append({"center": analysis["center"], "analysis": analysis["analysis"],
                "version": version, "N": analysis["N"], "repetitions": diagnostic["n_resamples"],
                "shared_tail_and_bin_map_combinations": diagnostic["distinct_tail_and_bin_mapping_combinations"],
                "shared_quantile_endpoint_and_node_combinations": diagnostic["distinct_quantile_node_combinations"],
                "distinct_D_values": summary["distinct_discordance_rates"],
                "all_D_equal": summary["distinct_discordance_rates"] == 1,
                "D_central_interval_degenerate": summary["discordance_interval_degenerate"],
                "all_D_zero": summary["all_resampled_discordance_zero"],
                "D_lower": intervals["discordance"][0], "D_upper": intervals["discordance"][1],
                "upward_lower": intervals["upward"][0], "upward_upper": intervals["upward"][1],
                "all_upward_equal": "not assessed in retained aggregate diagnostics",
                "failed_replicates": diagnostic["failed_replicates"]})
    pd.DataFrame(diagnostics).to_csv(output / "results/resampling-diagnostics.csv", index=False)
    return results

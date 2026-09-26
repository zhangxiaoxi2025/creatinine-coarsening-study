"""Synthetic adapter checks; no participant files are read."""
from decimal import Decimal, localcontext
import unittest

import numpy as np
import pandas as pd

from analysis_adapter import build_analysis, exact_units, UNIT_SCALE
from aki_core import evaluate_patient, fit_transforms, prepare_patient
from paired_resampling_exact import ANALYSIS_ORDER


def synthetic_case(key, records=None, *, start=0, end=100, discharge=700000,
                   diagnosis=False, conflict=False):
    records = records or [(-1, "1.00"), (101, "1.30")]
    ordered = sorted((Decimal(str(t)), str(v)) for t, v in records)
    lower, upper = Decimal(str(start)) - 604800, min(Decimal(str(end)) + 604800, Decimal(str(discharge)))
    labs = []
    for t, value in ordered:
        post = end < t <= upper
        previous = [s for s, _ in ordered if s < t]
        labs.append(dict(case_id=key, patient_id="p" + key, start=start, end=end,
                         discharge=discharge, time=t, value_str=value,
                         is_pre7=lower <= t < start, is_post=post,
                         E48=post and any(t - 172800 <= s < t for s in previous),
                         E7=post and any(t - 604800 <= s < t for s in previous)))
    baseline = any(l["is_pre7"] for l in labs)
    e48 = any(l["E48"] for l in labs)
    e7 = any(l["E7"] for l in labs)
    primary, secondary = baseline and e48 and not diagnosis, baseline and e7 and not diagnosis
    cohort = dict(case_id=key, patient_id="p" + key, start=start, end=end, discharge=discharge,
                  documented_renal_exclusion=diagnosis, primary=primary, secondary7=secondary,
                  procedural_only_primary=baseline and e48, conflict_in_window=conflict,
                  extreme_in_window=any(Decimal(v) < Decimal(".1") or Decimal(v) > 30 for _, v in ordered),
                  duration_over24h=end-start > 86400,
                  primary_after5min=primary and any(l["E48"] and l["time"]-end > 300 for l in labs),
                  secondary7_after5min=secondary and any(l["E7"] and l["time"]-end > 300 for l in labs))
    return cohort, labs


def frames(*cases):
    return pd.DataFrame([c for c, _ in cases]), pd.DataFrame([l for _, labs in cases for l in labs])


def build(cases, mode="main"):
    return build_analysis(*frames(*cases), analysis=mode)


class AdapterTests(unittest.TestCase):
    def test_lexicographic_original_strings_and_reference_order(self):
        actual = build([synthetic_case("2"), synthetic_case("10"), synthetic_case("01")])
        self.assertEqual(actual.patient_keys, ("01", "10", "2"))
        np.testing.assert_array_equal(actual.design.reference_patient, [0, 0, 1, 1, 2, 2])
        self.assertEqual(actual.diagnostics["target_points"], 3)

    def test_reference_includes_preop_intraop_and_E_outside(self):
        c = synthetic_case("1", [(-1, "0.9"), (100, "1.0"), (200000, "1.1"), (200001, "1.4")])
        a = build([c])
        self.assertEqual(a.diagnostics["reference_measurements"], 4)
        self.assertEqual(a.diagnostics["target_points"], 1)
        self.assertEqual(a.diagnostics["postop_points_outside_E48"], 1)
        np.testing.assert_array_equal(a.design.support[a.design._min48_code], [11 * UNIT_SCALE // 10])

    def test_min7_is_distinct_from_min48(self):
        c = synthetic_case("1", [(-200000, ".7"), (-1, "1.0"), (101, "1.2")])
        a = build([c])
        np.testing.assert_array_equal(a.design.support[a.design._min7_code], [7 * UNIT_SCALE // 10])
        np.testing.assert_array_equal(a.design.support[a.design._min48_code], [UNIT_SCALE])

    def test_prior_48h_lower_endpoint_included(self):
        a = build([synthetic_case("1", [(-172699, "1.0"), (101, "1.3")])])
        self.assertEqual(a.diagnostics["target_points"], 1)

    def test_parallel7_includes_patient_without_E48(self):
        both = synthetic_case("a")
        only7 = synthetic_case("b", [(-200000, "1.0"), (101, "1.5")])
        a = build([both, only7], "parallel7")
        self.assertEqual(a.patient_keys, ("a", "b"))
        self.assertEqual(a.design.mode, "ratio")
        self.assertIsNone(a.design._min48_code)
        self.assertEqual(build([both, only7]).patient_keys, ("a",))

    def test_S1_reinstates_diagnosis_but_no_others(self):
        included, diagnosed = synthetic_case("a"), synthetic_case("b", diagnosis=True)
        self.assertEqual(build([included, diagnosed]).patient_keys, ("a",))
        self.assertEqual(build([included, diagnosed], "S1_no_diagnosis").patient_keys, ("a", "b"))

    def test_S2_keeps_first5min_as_reference_and_prior(self):
        c = synthetic_case("1", [(-1, "2"), (101, "1"), (400, ".8"), (401, "1.1")])
        a = build([c], "S2_target_after5min")
        self.assertEqual(a.diagnostics["target_points"], 1)
        self.assertEqual(a.diagnostics["first5min_target_points"], 0)
        self.assertEqual(a.diagnostics["first5min_reference_measurements_retained"], 2)
        self.assertEqual(a.diagnostics["reference_measurements"], 4)
        np.testing.assert_array_equal(a.design.support[a.design._min48_code], [8 * UNIT_SCALE // 10])

    def test_S2_excludes_patient_with_only_early_targets(self):
        a = build([synthetic_case("early"), synthetic_case("late", [(-1, "1"), (401, "1.3")])], "S2_target_after5min")
        self.assertEqual(a.patient_keys, ("late",))
        self.assertEqual(a.diagnostics["mode_excluded_patients"], 1)

    def test_S3_excludes_whole_patient_and_preserves_endpoints(self):
        cases = [synthetic_case(str(i), [(-1, x), (101, "1")]) for i, x in enumerate([".09", ".1", "30", "30.01"])]
        a = build(cases, "S3_exclude_extremes")
        self.assertEqual(a.patient_keys, ("1", "2"))
        self.assertEqual(a.diagnostics["reference_measurements"], 4)

    def test_S4_preserves_exact_24hours(self):
        cases = [synthetic_case("a", [(-1, "1"), (86401, "1.3")], end=86400),
                 synthetic_case("b", [(-1, "1"), (86402, "1.3")], end=86401)]
        self.assertEqual(build(cases, "S4_exclude_over24h").patient_keys, ("a",))

    def test_S5_uses_original_conflict_flag(self):
        cases = [synthetic_case("a"), synthetic_case("b", conflict=True)]
        self.assertEqual(build(cases, "S5_exclude_conflicts").patient_keys, ("a",))

    def test_S6_requires_all_observed_postop_points_E48(self):
        cases = [synthetic_case("a"), synthetic_case("b", [(-1, "1"), (200000, "1"), (200001, "1.3")])]
        a = build(cases, "S6_all_postop_E48")
        self.assertEqual(a.patient_keys, ("a",))
        self.assertEqual(a.diagnostics["postop_points_outside_E48"], 0)

    def test_raw_99p99_value_not_lost_after_rounding(self):
        a = build([synthetic_case("1", [(-1, "99.99"), (101, "99.99")])])
        self.assertIn(9999 * UNIT_SCALE // 100, a.design.reference_units)
        self.assertIn(100 * UNIT_SCALE, a.design.fit_mappings().mappings["lower_bounded_round_0p2"])

    def test_extra_floor_diagnostic_distinct_from_output_at_floor(self):
        a = build([synthetic_case("1", [(-3, ".09"), (-2, ".1"), (-1, ".29"), (101, ".3")])])
        self.assertEqual(a.diagnostics["grid_extra_positive_floor_reference_measurements"], 1)
        self.assertEqual(a.diagnostics["grid_extra_positive_floor_patients"], 1)
        self.assertEqual(a.diagnostics["grid_output_equals_0p2_reference_measurements"], 3)

    def test_exact_units_rejects_hidden_precision_without_decimal_context_rounding(self):
        with localcontext() as context:
            context.prec = 3
            self.assertEqual(exact_units("99.99"), 9999 * UNIT_SCALE // 100)
            for value in ["1.0000000000000001", "1.000000000000000000000000000000000000001"]:
                with self.assertRaisesRegex(ValueError, "not exactly"):
                    exact_units(value)

    def test_mode_excluded_precision_does_not_constrain_selected_pool(self):
        common = synthetic_case("included", [(-1, "1"), (401, "1.3")])
        excluded = {
            "S2_target_after5min": synthetic_case("excluded", [(-1, "1.0000000000000001"), (101, "1.3")]),
            "S3_exclude_extremes": synthetic_case("excluded", [(-1, ".0900000000000001"), (101, "1.3")]),
            "S4_exclude_over24h": synthetic_case("excluded", [(-1, "1.0000000000000001"), (86402, "1.3")], end=86401),
            "S5_exclude_conflicts": synthetic_case("excluded", [(-1, "1.0000000000000001"), (101, "1.3")], conflict=True),
            "S6_all_postop_E48": synthetic_case("excluded", [(-1, "1.0000000000000001"), (200000, "1.3"), (200001, "1.4")]),
        }
        for mode, excluded_case in excluded.items():
            with self.subTest(mode=mode):
                self.assertEqual(build([common, excluded_case], mode).patient_keys, ("included",))
                with self.assertRaisesRegex(ValueError, "not exactly"):
                    build([common, excluded_case])

    def test_exact_15digit_values_are_preserved(self):
        self.assertEqual(exact_units("1.299999999999999"), 1299999999999999)
        self.assertEqual(exact_units("1.300000000000001"), 1300000000000001)
        with localcontext() as context:
            context.prec = 3
            self.assertEqual(exact_units("0.333333333333333"), 333333333333333)

    def test_15digit_near_absolute_threshold_keeps_distinct_labels(self):
        cases = [synthetic_case("below", [(-1, "1"), (101, "1.299999999999999")]),
                 synthetic_case("equal", [(-1, "1"), (101, "1.3")]),
                 synthetic_case("above", [(-1, "1"), (101, "1.300000000000001")])]
        a = build(cases)
        self.assertEqual(a.patient_keys, ("above", "below", "equal"))
        np.testing.assert_array_equal(a.design.evaluate().labels["identity"], [True, False, True])

    def test_15digit_near_ratio_threshold_keeps_distinct_labels(self):
        cases = [synthetic_case("below", [(-200000, "1"), (101, "1.499999999999999")]),
                 synthetic_case("equal", [(-200000, "1"), (101, "1.5")]),
                 synthetic_case("above", [(-200000, "1"), (101, "1.500000000000001")])]
        a = build(cases, "parallel7")
        np.testing.assert_array_equal(a.design.evaluate().labels["identity"], [True, False, True])

    def test_float_value_str_forbidden(self):
        with self.assertRaisesRegex(ValueError, "original decimal text"):
            exact_units(1.0)

    def test_invalid_original_values_forbidden(self):
        for value in ["0", "-1", "100", "NaN", "Infinity"]:
            with self.subTest(value=value), self.assertRaises(ValueError):
                exact_units(value)

    def test_observation_flag_drift_stops(self):
        c, m = frames(synthetic_case("1"))
        m.loc[1, "E48"] = False
        with self.assertRaisesRegex(ValueError, "observation flags"):
            build_analysis(c, m, analysis="main")

    def test_cohort_flag_drift_stops(self):
        c, m = frames(synthetic_case("1"))
        c.loc[0, "extreme_in_window"] = True
        with self.assertRaisesRegex(ValueError, "cohort flags"):
            build_analysis(c, m, analysis="main")

    def test_duplicate_and_conflicting_frozen_timestamp_stop(self):
        c, m = frames(synthetic_case("1"))
        for altered in (False, True):
            more = m.iloc[[0]].copy()
            if altered:
                more["value_str"] = "2.00"
            with self.subTest(altered=altered), self.assertRaisesRegex(ValueError, "deduplicated"):
                build_analysis(c, pd.concat([m, more]), analysis="main")

    def test_identity_linkage_and_metadata_drift_stop(self):
        c, original = frames(synthetic_case("1"))
        for column, value in (("patient_id", "other"), ("start", 1), ("end", 101), ("discharge", 700001)):
            m = original.copy()
            m.loc[0, column] = value
            with self.subTest(column=column), self.assertRaises(ValueError):
                build_analysis(c, m, analysis="main")

    def test_global_window_is_asserted_not_recropped(self):
        c, m = frames(synthetic_case("1"))
        m.loc[0, "time"] = -604801
        with self.assertRaisesRegex(ValueError, "outside the global"):
            build_analysis(c, m, analysis="main")

    def test_empty_mode_unknown_mode_and_nonboolean_flags_stop(self):
        c, m = frames(synthetic_case("1"))
        with self.assertRaisesRegex(ValueError, "No eligible"):
            build_analysis(c, m, analysis="S2_target_after5min")
        with self.assertRaisesRegex(ValueError, "Unknown analysis"):
            build_analysis(c, m, analysis="invented")
        c["primary"] = "False"
        with self.assertRaisesRegex(ValueError, "actual booleans"):
            build_analysis(c, m, analysis="main")

    def test_unique_original_string_identity_contract(self):
        c, m = frames(synthetic_case("1"))
        with self.assertRaisesRegex(ValueError, "one unique"):
            build_analysis(pd.concat([c, c]), m, analysis="main")
        c["case_id"] = 1
        with self.assertRaisesRegex(ValueError, "original nonempty strings"):
            build_analysis(c, m, analysis="main")

    def test_dataframes_not_mutated(self):
        c, m = frames(synthetic_case("1"))
        c_before, m_before = c.copy(deep=True), m.copy(deep=True)
        build_analysis(c, m, analysis="main")
        pd.testing.assert_frame_equal(c, c_before)
        pd.testing.assert_frame_equal(m, m_before)

    def test_scalar_core_agreement_across_all_transforms_main_and_parallel(self):
        cases = [synthetic_case("a", [(-200000, ".8"), (-1, "1"), (101, "1.2"), (400, "1.3")]),
                 synthetic_case("b", [(-1, ".1"), (101, ".2"), (604900, ".4")]),
                 synthetic_case("c", [(-1, "1.9"), (101, "2.1"), (102, "2.2")])]
        for mode in ("main", "parallel7"):
            adapted = build(cases, mode)
            reference = [value for _, rows in cases for value in (r["value_str"] for r in rows)]
            transforms = fit_transforms(reference)
            actual = adapted.design.evaluate()
            for version, transform in transforms.versions().items():
                expected = []
                for cohort, labs in cases:
                    plan = prepare_patient([(r["time"], r["value_str"]) for r in labs],
                                           aneend=cohort["end"], discharge=cohort["discharge"])
                    result = evaluate_patient(plan, transform)
                    expected.append(result.secondary_ratio_label if mode == "parallel7" else result.primary_label)
                with self.subTest(mode=mode, version=version):
                    np.testing.assert_array_equal(actual.labels[version], expected)

    def test_all_eight_frozen_modes_supported(self):
        case = synthetic_case("1", [(-1, "1"), (401, "1.3")])
        self.assertEqual(len(ANALYSIS_ORDER), 8)
        for mode in ANALYSIS_ORDER:
            self.assertEqual(build([case], mode).patient_keys, ("1",))



if __name__ == "__main__":
    unittest.main(verbosity=2)

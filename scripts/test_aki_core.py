"""Synthetic boundary tests only: no patient files are opened by this suite."""

import hashlib
import json
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
import sys
import unittest

from aki_core import (
    clean_series, decimal, evaluate_patient, fit_transforms,
    lower_bounded_round_0p2, prepare_patient, type1_quantile,
)

D = Decimal


def result(records, end=0, discharge=604800):
    return evaluate_patient(prepare_patient(records, aneend=end, discharge=discharge))


class WindowAndLabelTests(unittest.TestCase):
    def test_absolute_exact_boundary_float_trap(self):
        self.assertLess(1.4 - 1.1, .3)
        r = result([(-1, 1.1), (1, 1.4)])
        self.assertTrue(r.primary_label)
        self.assertTrue(r.primary_points[0].absolute_positive)
        self.assertFalse(r.primary_points[0].ratio_positive)

    def test_just_below_absolute_threshold(self):
        r = result([(-1, '1.1'), (1, '1.399999999999999999999999999999')])
        self.assertFalse(r.primary_label)

    def test_ratio_exact_boundary(self):
        r = result([(-1, '.2'), (1, '.3')])
        self.assertTrue(r.primary_label)
        self.assertFalse(r.primary_points[0].absolute_positive)
        self.assertTrue(r.primary_points[0].ratio_positive)

    def test_just_below_ratio_threshold(self):
        self.assertFalse(result([(-1, '.2'), (1, '.299999999999999999999999999999')]).primary_label)

    def test_48h_inclusive_left(self):
        p = prepare_patient([(1 - 172800, '1'), (1, '1.3')], aneend=0, discharge=2)
        self.assertEqual(p.primary_points[0].prior48, (0,))
        self.assertTrue(evaluate_patient(p).primary_label)

    def test_48h_one_second_outside(self):
        r = result([(1 - 172801, '1'), (1, '1.5')])
        self.assertIsNone(r.primary_label)
        self.assertTrue(r.secondary_ratio_label)

    def test_7d_inclusive_left(self):
        r = result([(1 - 604800, '1'), (1, '1.5')])
        self.assertIsNone(r.primary_label)
        self.assertTrue(r.secondary_ratio_label)

    def test_7d_one_second_outside(self):
        r = result([(1 - 604801, '1'), (1, '2')])
        self.assertIsNone(r.primary_label)
        self.assertIsNone(r.secondary_ratio_label)

    def test_window_minimum_not_latest(self):
        r = result([(-2, '1'), (-1, '1.25'), (1, '1.35')])
        self.assertEqual(r.primary_points[0].min48, D('1'))
        self.assertTrue(r.primary_label)

    def test_separate_window_minima(self):
        r = result([(1 - 604800, '.2'), (1 - 172800, '.5'), (1, '.3')])
        p = r.primary_points[0]
        self.assertEqual((p.min48, p.min7), (D('.5'), D('.2')))
        self.assertFalse(p.absolute_positive)
        self.assertTrue(p.ratio_positive)

    def test_only_E_primary_points_contribute_to_primary_or(self):
        r = result([(0, '1'), (1, '1'), (200000, '2')])
        self.assertEqual([p.time for p in r.primary_points], [D(1)])
        self.assertFalse(r.primary_label)
        self.assertTrue(r.secondary_ratio_label)

    def test_postoperative_and_discharge_boundaries(self):
        p = prepare_patient([(-1, '1'), (0, '5'), (1, '1'), (2, '1'), (3, '9')], aneend=0, discharge=2)
        self.assertEqual(tuple(p.cleaned.measurements[i].time for i in p.postoperative_indices), (D(1), D(2)))
        self.assertFalse(evaluate_patient(p).primary_label)

    def test_postop_7d_right_inclusive_one_second_after_excluded(self):
        p = prepare_patient([(604799, '1'), (604800, '1.3'), (604801, '9')], aneend=0, discharge=900000)
        self.assertEqual(tuple(p.cleaned.measurements[i].time for i in p.postoperative_indices), (D(604799), D(604800)))
        self.assertTrue(evaluate_patient(p).primary_label)

    def test_no_prior_is_none_not_negative(self):
        r = result([(1, '4')])
        self.assertIsNone(r.primary_label)
        self.assertIsNone(r.secondary_ratio_label)

    def test_empty_observation_is_none(self):
        self.assertIsNone(result([]).primary_label)
        self.assertIsNone(result([(-1, '1'), (0, '2')]).primary_label)

    def test_nonzero_time_origin_and_fractional_seconds(self):
        p = prepare_patient([('-173099.875', '1'), ('-299.875', '1.3')], aneend='-300', discharge='-200')
        self.assertTrue(evaluate_patient(p).primary_label)
        self.assertEqual(p.primary_points[0].prior48, (0,))

    def test_transformation_never_recomputes_E(self):
        p = prepare_patient([(0, '1'), (1, '1.3'), (200000, '2')], aneend=0, discharge=300000)
        tr = fit_transforms(['.2', '.4', '1', '1.3', '2'])
        original_times = tuple(v.time for v in evaluate_patient(p).primary_points)
        for name, fn in tr.versions().items():
            self.assertEqual(tuple(v.time for v in evaluate_patient(p, fn).primary_points), original_times, name)

    def test_transformed_values_can_change_label_but_not_sample_times(self):
        p = prepare_patient([(-1, '1.1'), (1, '1.4')], aneend=0, discharge=100)
        t = fit_transforms(['1.1', '1.4'])
        original = evaluate_patient(p, t.identity)
        rounded = evaluate_patient(p, lower_bounded_round_0p2)
        self.assertTrue(original.primary_label)
        self.assertFalse(rounded.primary_label)
        self.assertEqual([x.time for x in original.primary_points], [x.time for x in rounded.primary_points])

    def test_conversion_does_not_resolve_an_original_time_conflict(self):
        p = prepare_patient([(-1, '1.11'), (-1, '1.12'), (1, '1.4')], aneend=0, discharge=100)
        self.assertEqual(lower_bounded_round_0p2('1.11'), lower_bounded_round_0p2('1.12'))
        self.assertIsNone(evaluate_patient(p, lower_bounded_round_0p2).primary_label)

    def test_arbitrary_nonpositive_transform_fails_instead_of_rescreening(self):
        p = prepare_patient([(-1, '1'), (1, '1.3')], aneend=0, discharge=100)
        with self.assertRaises(ValueError):
            evaluate_patient(p, lambda x: D(0))

    def test_decimal_context_does_not_change_boundaries(self):
        with localcontext() as ctx:
            ctx.prec = 3
            r = result([(-1, '1.1'), (1, '1.399999999999999999999999999999')])
            self.assertFalse(r.primary_label)
            p = prepare_patient([
                ('99999999999999999999827201.001', '1'),
                ('100000000000000000000000001.001', '1.3'),
            ], aneend='100000000000000000000000000', discharge='100000000000000000000000002')
            self.assertEqual(p.primary_points[0].prior48, (0,))
            self.assertTrue(evaluate_patient(p).primary_label)


class CleaningTests(unittest.TestCase):
    def test_exact_duplicates_then_entire_conflict_time_excluded(self):
        c = clean_series([(1, '1'), (1, '1.0'), (1, '2'), (2, '1.3'), (2, '1.30')])
        self.assertEqual(c.exact_duplicates_removed, 2)
        self.assertEqual(c.conflicting_times, (D(1),))
        self.assertEqual(c.conflicting_distinct_records_removed, 2)
        self.assertEqual([(m.time, m.value) for m in c.measurements], [(D(2), D('1.3'))])

    def test_conflicted_prior_does_not_supply_observability(self):
        r = result([(0, '1'), (0, '2'), (1, '3')])
        self.assertIsNone(r.primary_label)

    def test_conflicted_current_does_not_supply_an_outcome(self):
        r = result([(-1, '1'), (1, '1.3'), (1, '2')])
        self.assertIsNone(r.primary_label)

    def test_unordered_records_sort_by_exact_time(self):
        c = clean_series([(2, '1.2'), (-1, '1'), (1, '1.1')])
        self.assertEqual([m.time for m in c.measurements], [D(-1), D(1), D(2)])

    def test_invalid_values_and_endpoints_raise(self):
        for value in ['nan', 'inf', '-inf', '0', '-1', True]:
            with self.subTest(value=value), self.assertRaises(ValueError):
                clean_series([(0, value)])
        with self.assertRaises(ValueError):
            clean_series([('nan', '1')])
        with self.assertRaises(ValueError):
            prepare_patient([], aneend=1, discharge=0)

    def test_float_contract_is_explicit_not_heuristic_rounding(self):
        self.assertEqual(decimal(.1), D('.1'))
        self.assertEqual(decimal(.1 + .2), D('.30000000000000004'))


class NumericTransformationTests(unittest.TestCase):
    def test_type1_integer_and_fractional_rank(self):
        values = tuple(D(i) for i in range(1, 41))
        self.assertEqual(type1_quantile(values, Fraction(1, 40)), D(1))
        self.assertEqual(type1_quantile(values, Fraction(39, 40)), D(39))
        self.assertEqual(type1_quantile(values, Fraction(51, 1000)), D(3))
        self.assertEqual(type1_quantile(values, Fraction(0)), D(1))
        self.assertEqual(type1_quantile(values, Fraction(1)), D(40))

    def test_identity_is_exact(self):
        t = fit_transforms(['.01', '.3', '1.4'])
        for value in ['.00000001', '.30000000000000001', '10000000000000000.11111111']:
            self.assertEqual(t.identity(value), D(value))

    def test_strict_tails_and_equal_endpoints(self):
        t = fit_transforms([str(i) for i in range(1, 41)])
        for fn in [t.tail_only, t.tail_plus_bins]:
            self.assertEqual(fn('.1'), D(1))
            self.assertEqual(fn('1'), D(1))
            self.assertEqual(fn('39'), D(39))
            self.assertEqual(fn('40'), D(39))
        self.assertEqual(t.tail_only('1.5'), D('1.5'))
        self.assertEqual(t.tail_plus_bins('1.5'), D(2))

    def test_midrank_nearest_node_tie_goes_lower(self):
        # x=2: (one lower + half one equal)/20=.075; .05 wins tie with .10.
        t = fit_transforms([str(i) for i in range(1, 21)])
        self.assertEqual(t.tail_plus_bins('2'), D(1))
        self.assertEqual(t.tail_plus_bins('3'), D(2))

    def test_parallel_values_use_original_reference_midrank(self):
        # Five observations equal 2: rank=(1+2.5)/10=.35 -> type1 q35 = 2.
        t = fit_transforms(['1', '2', '2', '2', '2', '2', '4', '5', '6', '7'])
        self.assertEqual(t.tail_plus_bins('2'), D(2))
        self.assertEqual(len(t.reference), 10)

    def test_singleton_and_all_tied_reference(self):
        for values in [['.8'], ['.8'] * 40]:
            t = fit_transforms(values)
            self.assertEqual(t.tail_plus_bins('.8'), D('.8'))
            self.assertEqual(t.tail_plus_bins('.1'), D('.8'))
            self.assertEqual(t.tail_plus_bins('2'), D('.8'))

    def test_two_value_reference_preserves_positive_monotonic_results(self):
        t = fit_transforms(['.01', '2'])
        output = [t.tail_plus_bins(x) for x in ['.001', '.01', '.1', '1', '2', '3']]
        self.assertEqual(output, sorted(output))
        self.assertTrue(all(x > 0 for x in output))

    def test_lower_bounded_half_up_control(self):
        expected = {
            '.00001': '.2', '.02': '.2', '.099999': '.2', '.1': '.2',
            '.2': '.2', '.29999999999999999999999999': '.2', '.3': '.4',
            '.5': '.6', '1.1': '1.2',
        }
        for value, target in expected.items():
            self.assertEqual(lower_bounded_round_0p2(value), D(target))

    def test_rounding_independent_of_decimal_context(self):
        with localcontext() as ctx:
            ctx.prec = 3
            self.assertEqual(lower_bounded_round_0p2('123456789.3'), D('123456789.4'))

    def test_all_versions_monotone_positive_for_many_ties_and_small_pools(self):
        pools = [['.02'], ['.02', '.03'], ['.02', '.1', '.1', '.2', '.7', '1', '1', '2', '30']]
        pools.append([str(D(i) / 100) for i in range(1, 301)])
        grid = [D(i) / 100 for i in range(1, 3101)]
        for pool in pools:
            t = fit_transforms(pool)
            for name, fn in t.versions().items():
                output = [fn(x) for x in grid]
                self.assertTrue(all(x > 0 for x in output), name)
                self.assertEqual(output, sorted(output), name)

    def test_empty_or_nonpositive_reference_rejected(self):
        for values in [[], [0], ['nan'], [-1]]:
            with self.assertRaises(ValueError):
                fit_transforms(values)


class RecordingResult(unittest.TextTestResult):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.records = []

    def addSuccess(self, test):
        super().addSuccess(test)
        self.records.append({'test': test.id(), 'status': 'PASS'})

    def addFailure(self, test, err):
        super().addFailure(test, err)
        self.records.append({'test': test.id(), 'status': 'FAIL'})

    def addError(self, test, err):
        super().addError(test, err)
        self.records.append({'test': test.id(), 'status': 'ERROR'})



if __name__ == "__main__":
    unittest.main(verbosity=2)

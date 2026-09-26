"""Synthetic-only tests; no participant data or result files are accessed."""
from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

import numpy as np

from extension_core import (MAX_SAFE_UNITS, criterion_states, four_by_four,
                            paired_table, patient_equal_mapping)


ROOT = Path(__file__).resolve().parents[1]
ORIGINAL_CORE = ROOT / "scripts/paired_resampling_exact.py"
spec = importlib.util.spec_from_file_location("original_exact_core", ORIGINAL_CORE)
original_core = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = original_core
spec.loader.exec_module(original_core)


def independent_weighted_oracle(patient, values, n):
    """Dictionary mass, linear CDF scan and all-node distance comparison."""
    sizes = Counter(map(int, patient))
    mass = defaultdict(Fraction)
    for person, value in zip(patient, values):
        mass[int(value)] += Fraction(1, sizes[int(person)] * n)
    support = sorted(mass)

    def quantile(p):
        acc = Fraction(0)
        for value in support:
            acc += mass[value]
            if acc >= p:
                return value
        raise AssertionError("incomplete reference mass")

    nodes = [quantile(Fraction(k, 20)) for k in range(1, 20)]
    lower, upper = quantile(Fraction(1, 40)), quantile(Fraction(39, 40))
    below, mapped = Fraction(0), []
    for value in support:
        midrank = below + mass[value] / 2
        node = min(range(1, 20), key=lambda k: (abs(Fraction(k, 20) - midrank), k))
        mapped.append(lower if value <= lower else upper if value >= upper else nodes[node - 1])
        below += mass[value]
    return support, mapped, lower, upper, nodes


def explicit_design(current, minimum48, minimum7, patient=None, scale=10**15):
    current, minimum48, minimum7 = [np.asarray(x, dtype=np.int64)
                                   for x in (current, minimum48, minimum7)]
    patient = np.arange(len(current), dtype=np.int64) if patient is None else np.asarray(patient, dtype=np.int64)
    support = np.unique(np.concatenate((current, minimum48, minimum7)))
    return SimpleNamespace(n_patients=int(patient.max()) + 1, unit_scale=scale,
                           point_patient=patient, current_units=current,
                           min48_units=minimum48, min7_units=minimum7,
                           support=support, mode="composite")


class WeightedMappingTests(unittest.TestCase):
    def test_weighted_cdf_exact_boundary(self):
        fit = patient_equal_mapping(np.array([0, 0, 0, 1]), np.array([10, 20, 30, 40]), 2)
        self.assertEqual(fit["support_weights"], (Fraction(1, 3),) * 3 + (Fraction(1),))
        self.assertEqual(fit["normalized_weights"], (Fraction(1, 6),) * 3 + (Fraction(1, 2),))
        self.assertEqual(fit["node_values"][9], 30)  # F(30) is exactly 0.5.
        self.assertEqual(fit["node_values"][10], 40)
        self.assertEqual(fit["total_weight"], 2)

    def test_exact_midrank_ties_choose_lower_node(self):
        fit = patient_equal_mapping(np.arange(20), np.arange(1, 21), 20)
        self.assertEqual(fit["weighted_midranks"][1], Fraction(3, 40))
        self.assertEqual(fit["node_values"], list(range(1, 20)))
        self.assertEqual(fit["nearest_node_indices"][1], 1)
        self.assertEqual(fit["mapped_units"][1], 1)
        self.assertEqual(fit["nearest_node_indices"][10], 10)
        self.assertEqual(fit["mapped_units"][10], 10)

    def test_tail_cdf_exact_boundaries(self):
        fit = patient_equal_mapping(np.arange(40), np.arange(1, 41), 40)
        self.assertEqual((fit["q025"], fit["q975"]), (1, 39))
        self.assertEqual(fit["mapped_units"][-1], 39)

    def test_constant_support(self):
        fit = patient_equal_mapping(np.array([0, 0, 1, 2, 2, 2]), np.full(6, 777), 3)
        np.testing.assert_array_equal(fit["support"], [777])
        np.testing.assert_array_equal(fit["mapped_units"], [777])
        self.assertEqual(fit["support_weights"], (Fraction(3),))
        self.assertEqual(fit["weighted_midranks"], (Fraction(1, 2),))
        self.assertEqual(fit["node_values"], [777] * 19)

    def test_unequal_counts_preserve_one_total_weight_per_patient(self):
        patient = np.array([0] * 100 + [1])
        values = np.array([10] * 100 + [100])
        fit = patient_equal_mapping(patient, values, 2)
        self.assertEqual(fit["support_weights"], (Fraction(1), Fraction(1)))
        self.assertEqual((fit["q025"], fit["q975"]), (10, 100))
        design = original_core.PairedDesign(n_patients=2, reference_patient=patient,
                   reference_units=values, point_patient=np.array([0, 1]),
                   current_units=np.array([10, 100]), min48_units=np.array([10, 100]),
                   min7_units=np.array([10, 100]), unit_scale=100)
        self.assertEqual(design.fit_mappings().q975, 10)

    def test_equal_record_counts_match_original_uniform_map(self):
        rng = np.random.default_rng(721)
        for _ in range(40):
            patient = np.repeat(np.arange(13), 11)
            values = rng.integers(1, 90, size=len(patient), dtype=np.int64)
            current = values[::11]
            design = original_core.PairedDesign(n_patients=13, reference_patient=patient,
                       reference_units=values, point_patient=np.arange(13),
                       current_units=current, min48_units=current, min7_units=current,
                       unit_scale=100)
            baseline = design.fit_mappings()
            fit = patient_equal_mapping(patient, values, 13)
            np.testing.assert_array_equal(fit["support"], design.support)
            np.testing.assert_array_equal(fit["mapped_units"], baseline.mappings["tail_plus_bins"])
            self.assertEqual((fit["q025"], fit["q975"]), (baseline.q025, baseline.q975))
            self.assertEqual(fit["node_values"], baseline.node_values.tolist())

    def test_unequal_counts_match_independent_rational_oracle(self):
        rng = np.random.default_rng(993)
        for _ in range(75):
            counts = rng.integers(1, 14, size=9)
            patient = np.repeat(np.arange(9), counts)
            values = rng.integers(1, 65, size=len(patient), dtype=np.int64)
            fit = patient_equal_mapping(patient, values, 9)
            support, mapped, lower, upper, nodes = independent_weighted_oracle(patient, values, 9)
            np.testing.assert_array_equal(fit["support"], support)
            np.testing.assert_array_equal(fit["mapped_units"], mapped)
            self.assertEqual((fit["q025"], fit["q975"], fit["node_values"]), (lower, upper, nodes))

    def test_record_and_patient_permutation_invariance(self):
        patient = np.repeat(np.arange(5), [2, 4, 3, 8, 1])
        values = np.array([5, 10, 7, 8, 11, 17, 2, 5, 19, 1, 2, 3, 4, 5, 6, 7, 9, 18])
        fit = patient_equal_mapping(patient, values, 5)
        rng = np.random.default_rng(771)
        for _ in range(15):
            order = rng.permutation(len(patient))
            relabel = rng.permutation(5)
            changed = patient_equal_mapping(relabel[patient[order]], values[order], 5)
            for key in ("support", "mapped_units"):
                np.testing.assert_array_equal(fit[key], changed[key])
            for key in ("support_weights", "normalized_weights", "weighted_midranks", "node_values"):
                self.assertEqual(fit[key], changed[key])

    def test_fit_fail_closed_for_inexact_missing_or_out_of_range_input(self):
        for patient, values, n in (([0, 1], [1.0, 2.0], 2), ([0, 0], [1, 2], 2),
                                   ([0, 2], [1, 2], 2), ([0, 1], [0, 2], 2),
                                   ([0, 1], [1, MAX_SAFE_UNITS + 1], 2),
                                   ([0, 1], [1, 2], True)):
            with self.subTest(patient=patient, values=values, n=n), self.assertRaises(ValueError):
                patient_equal_mapping(np.asarray(patient), np.asarray(values), n)


class CriterionStateTests(unittest.TestCase):
    def test_exact_absolute_and_ratio_threshold_equality(self):
        s = 10**15
        design = explicit_design([13*s//10, 13*s//10-1, 45*s//100, 45*s//100-1, 15*s//10],
                                 [s, s, 3*s//10, 3*s//10, s],
                                 [s, s, 3*s//10, 3*s//10, s])
        states = criterion_states(design, design.support, design.support)
        np.testing.assert_array_equal(states, [2, 0, 1, 0, 3])
        self.assertEqual(states.dtype, np.int8)

    def test_absolute_and_ratio_can_occur_at_different_times(self):
        design = explicit_design([130, 45, 120], [100, 30, 100], [100, 30, 100],
                                 patient=[0, 0, 1], scale=100)
        states = criterion_states(design, design.support, design.support)
        np.testing.assert_array_equal(states, [3, 0])
        simultaneous = ((design.current_units - design.min48_units >= 30)
                        & (2 * design.current_units >= 3 * design.min7_units))
        self.assertFalse(simultaneous.any())

    def test_fixed_point_order_does_not_change_patient_or(self):
        design = explicit_design([130, 45, 120, 190], [100, 30, 100, 170], [100, 30, 100, 170],
                                 patient=[0, 0, 1, 1], scale=100)
        baseline = criterion_states(design, design.support, design.support)
        order = np.array([3, 1, 2, 0])
        for field in ("point_patient", "current_units", "min48_units", "min7_units"):
            setattr(design, field, getattr(design, field)[order])
        np.testing.assert_array_equal(criterion_states(design, design.support, design.support), baseline)

    def test_original_coded_design_and_or_labels_match(self):
        rng = np.random.default_rng(888)
        n, m = 17, 8
        patient = np.repeat(np.arange(n), m)
        values = rng.integers(20, 201, size=len(patient), dtype=np.int64)
        minimum48 = np.minimum(values, values.reshape(n, m).min(axis=1).repeat(m))
        minimum7 = minimum48.copy()
        design = original_core.PairedDesign(n_patients=n, reference_patient=patient,
                   reference_units=values, point_patient=patient, current_units=values,
                   min48_units=minimum48, min7_units=minimum7, unit_scale=100)
        evaluated = design.evaluate()
        for name, mapping in evaluated.mapping_fit.mappings.items():
            states = criterion_states(design, design.support, mapping)
            np.testing.assert_array_equal(states > 0, evaluated.labels[name])
        identity = criterion_states(design, design.support, design.support)
        transformed = criterion_states(design, design.support,
                                       evaluated.mapping_fit.mappings["tail_plus_bins"])
        table = paired_table(identity, transformed)
        for key, value in evaluated.tables()["tail_plus_bins"].items():
            self.assertEqual(table[key], value)

    def test_monotone_transformed_minimum_is_valid(self):
        original = np.array([10, 11, 12, 19, 20])
        mapped = np.array([10, 10, 15, 15, 15])
        for indices in ([1, 3], [0, 2, 4], [2, 3, 4]):
            self.assertEqual(mapped[indices].min(), mapped[np.searchsorted(original, original[indices].min())])

    def test_state_validation_rejects_changed_support_nonmonotone_and_overflow(self):
        design = explicit_design([130, 45], [100, 30], [100, 30], scale=100)
        with self.assertRaises(ValueError):
            criterion_states(design, design.support[::-1], design.support)
        with self.assertRaises(ValueError):
            criterion_states(design, design.support, design.support[::-1])
        with self.assertRaises(ValueError):
            criterion_states(design, design.support[:-1], design.support[:-1])
        with self.assertRaises(ValueError):
            criterion_states(design, design.support, np.full(len(design.support), MAX_SAFE_UNITS + 1))
        design.unit_scale = 101
        with self.assertRaises(ValueError):
            criterion_states(design, design.support, design.support)

    def test_missing_point_patient_or_invalid_nested_minimum_rejected(self):
        design = explicit_design([130, 45], [100, 30], [100, 30], scale=100)
        design.n_patients = 3
        with self.assertRaises(ValueError):
            criterion_states(design, design.support, design.support)
        design.n_patients = 2
        design.min7_units = np.array([130, 30])
        with self.assertRaises(ValueError):
            criterion_states(design, design.support, design.support)


class PairedSummaryTests(unittest.TestCase):
    def test_all_sixteen_transitions_and_or_collapse(self):
        original, transformed = np.repeat(np.arange(4), 4), np.tile(np.arange(4), 4)
        matrix = four_by_four(original, transformed)
        np.testing.assert_array_equal(matrix, np.ones((4, 4), dtype=np.int64))
        table = paired_table(original, transformed)
        self.assertEqual([table[key] for key in ("n00", "n01", "n10", "n11")], [1, 3, 3, 9])
        self.assertEqual(table["discordance_count"], 6)
        self.assertEqual(table["discordance"], 6/16)
        self.assertEqual(table["net_change"], 0)
        self.assertEqual(table["original_positive_count"], 12)
        self.assertEqual(table["transformed_positive_count"], 12)

    def test_fixed_coverage_partition_counts_sum_and_preserve_direction(self):
        original = np.array([0, 0, 1, 2, 3, 0, 2])
        transformed = np.array([0, 3, 0, 1, 3, 2, 0])
        mask = np.array([True, False, True, True, False, True, False])
        total = paired_table(original, transformed)
        first, second = paired_table(original, transformed, mask), paired_table(original, transformed, ~mask)
        for key in ("n", "n00", "n01", "n10", "n11", "discordance_count", "upward_count",
                    "downward_count", "net_change_count", "original_positive_count", "transformed_positive_count"):
            self.assertEqual(first[key] + second[key], total[key])
        self.assertEqual((total["n01"], total["n10"]), (2, 2))

    def test_empty_subset_uses_null_rates(self):
        result = paired_table(np.array([0, 1]), np.array([1, 0]), np.array([False, False]))
        self.assertEqual(result["n"], 0)
        self.assertIsNone(result["discordance"])
        self.assertIsNone(result["net_change"])

    def test_bad_masks_states_and_lengths_rejected(self):
        with self.assertRaises(ValueError):
            paired_table(np.array([0, 1]), np.array([0, 1]), np.array([1, 0]))
        with self.assertRaises(ValueError):
            four_by_four(np.array([0, 4]), np.array([0, 1]))
        with self.assertRaises(ValueError):
            paired_table(np.array([0]), np.array([0, 1]))



if __name__ == "__main__":
    unittest.main(verbosity=2)

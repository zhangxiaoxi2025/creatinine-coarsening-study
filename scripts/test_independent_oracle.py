"""Independent synthetic-only oracle for aki_core; never reads participant data.

Oracle deliberately uses enumerated sets, exact fractions, and CDF counts rather
than the production bisect windows/quantile rank formula. Generated examples are
not clinical estimates and provide no information about either study dataset.
"""
from collections import Counter
from decimal import Decimal, localcontext
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import random
import unittest
import aki_core as core

ROOT = Path(__file__).resolve().parents[1]
COUNTS = {"generated_patient_versions": 0, "generated_mapping_queries": 0,
          "small_tie_reference_pools": 0}

def cstr(cents):
    return f"{cents // 100}.{cents % 100:02d}"

def quantile_by_cdf(ref, p):
    counts = Counter(ref)
    cumulative = 0
    for x in sorted(counts):
        cumulative += counts[x]
        if F(cumulative, len(ref)) >= p:
            return x
    raise AssertionError("unreachable quantile")

def mapping_oracle(ref, x, version):
    if version == "identity":
        return x
    if version == "lower_bounded_round_0p2":
        candidates = [F(k, 5) for k in range(1, int(x * 5) + 3)]
        # Nearest positive .2 lattice point, ties upward: equivalent semantics
        # independently calculated without the production floor expression.
        return min(candidates, key=lambda g: (abs(x - g), -g))
    lo, hi = quantile_by_cdf(ref, F(1, 40)), quantile_by_cdf(ref, F(39, 40))
    if x <= lo:
        return lo
    if x >= hi:
        return hi
    if version == "tail_only":
        return x
    cdf_less = sum(z < x for z in ref)
    equal = sum(z == x for z in ref)
    mid = F(cdf_less, len(ref)) + F(equal, 2 * len(ref))
    node = min([F(k, 20) for k in range(1, 20)], key=lambda p: (abs(p-mid), p))
    return quantile_by_cdf(ref, node)

def patient_oracle(records, end, discharge, tf=lambda x: x):
    by_time = {}
    for t, v in records:
        by_time.setdefault(F(t), set()).add(F(v))
    clean = sorted((t, next(iter(v))) for t, v in by_time.items() if len(v) == 1)
    upper = min(F(end) + 604800, F(discharge))
    primary, secondary = [], []
    for t, raw in clean:
        if not F(end) < t <= upper:
            continue
        p48 = [tf(v) for s, v in clean if t-F(172800) <= s < t]
        p7 = [tf(v) for s, v in clean if t-F(604800) <= s < t]
        if not p7:
            continue
        now = tf(raw)
        m7 = min(p7)
        rp = now >= F(3, 2) * m7
        ap = None if not p48 else now - min(p48) >= F(3, 10)
        row = (t, now, None if not p48 else min(p48), m7, ap, rp,
               None if ap is None else ap or rp)
        secondary.append(row)
        if p48:
            primary.append(row)
    label = any(p[-1] for p in primary) if primary else None
    ratio = any(p[-2] for p in secondary) if secondary else None
    return label, ratio, tuple(primary), tuple(secondary)

def normalize_actual(r):
    def row(p):
        return (F(p.time), F(p.value), None if p.min48 is None else F(p.min48),
                F(p.min7), p.absolute_positive, p.ratio_positive, p.composite_positive)
    return r.primary_label, r.secondary_ratio_label, tuple(map(row, r.primary_points)), tuple(map(row, r.secondary_points))

class IndependentOracle(unittest.TestCase):
    def test_exact_absolute_and_ratio_boundaries(self):
        pairs = [("1.00", "1.30", True), ("1.00", "1.2999999999999999999999999999", False),
                 ("0.20", "0.30", True), ("0.20", "0.2999999999999999999999999999", False)]
        for old, new, expected in pairs:
            records = [("-1", old), ("1", new)]
            actual = core.evaluate_patient(core.prepare_patient(records, aneend=0, discharge=100))
            self.assertEqual(actual.primary_label, expected)
            self.assertEqual(normalize_actual(actual), patient_oracle(records, 0, 100))

    def test_minimum_window_is_not_nearest(self):
        records = [("-100", "1.0"), ("-1", "1.4"), ("1", "1.3")]
        actual = core.evaluate_patient(core.prepare_patient(records, aneend=0, discharge=100))
        self.assertTrue(actual.primary_label)
        self.assertTrue(actual.primary_points[0].absolute_positive)
        self.assertEqual(actual.primary_points[0].min48, Decimal("1.0"))
        self.assertEqual(normalize_actual(actual), patient_oracle(records, 0, 100))

    def test_strict_prior_and_inclusive_left_window(self):
        for width in (172800, 604800):
            for epsilon in ("0", "0.0001", "-0.0001"):
                old_time = Decimal(1-width) - Decimal(epsilon)
                records = [(str(old_time), "1.0"), ("1", "1.5")]
                actual = core.evaluate_patient(core.prepare_patient(records, aneend=0, discharge=2))
                self.assertEqual(normalize_actual(actual), patient_oracle(records, 0, 2))
                if width == 172800:
                    self.assertEqual(actual.primary_label is not None, Decimal(epsilon) <= 0)
                if width == 604800:
                    self.assertEqual(actual.secondary_ratio_label is not None, Decimal(epsilon) <= 0)

    def test_end_exclusive_and_discharge_inclusive(self):
        records = [("-1", "1.0"), ("0", "1.5"), ("100", "1.0"), ("100.01", "2.0")]
        plan = core.prepare_patient(records, aneend=0, discharge=100)
        actual = core.evaluate_patient(plan)
        self.assertEqual([p.time for p in actual.primary_points], [Decimal(100)])
        self.assertEqual(normalize_actual(actual), patient_oracle(records, 0, 100))

    def test_postoperative_day_seven_inclusive(self):
        records = [("604799", "1.0"), ("604800", "1.3"), ("604800.0001", "2.0")]
        actual = core.evaluate_patient(core.prepare_patient(records, aneend=0, discharge=700000))
        self.assertEqual([p.time for p in actual.primary_points], [Decimal(604800)])
        self.assertEqual(normalize_actual(actual), patient_oracle(records, 0, 700000))

    def test_no_prior_is_unknown_and_e7_only_does_not_enter_primary(self):
        records = [("-200000", "1.0"), ("1", "1.5")]
        actual = core.evaluate_patient(core.prepare_patient(records, aneend=0, discharge=10))
        self.assertIsNone(actual.primary_label)
        self.assertTrue(actual.secondary_ratio_label)
        no_prior = core.evaluate_patient(core.prepare_patient([("1", "1.5")], aneend=0, discharge=10))
        self.assertIsNone(no_prior.primary_label)
        self.assertIsNone(no_prior.secondary_ratio_label)

    def test_e7_only_record_can_be_later_prior(self):
        # t=1 has only a 7-day prior. It remains a valid prior for t=2 in E48.
        records = [("-200000", "1.0"), ("1", "1.1"), ("2", "1.4")]
        actual = core.evaluate_patient(core.prepare_patient(records, aneend=0, discharge=10))
        self.assertEqual([p.time for p in actual.primary_points], [Decimal(2)])
        self.assertEqual(actual.primary_points[0].min48, Decimal("1.1"))
        self.assertEqual(normalize_actual(actual), patient_oracle(records, 0, 10))

    def test_duplicates_conflict_entire_timestamp_and_order(self):
        records = [("-1", "1.0"), ("-1", "1.00"), ("1", "1.0"), ("1", "1.3"), ("2", "1.1")]
        for incoming in (records, list(reversed(records))):
            plan = core.prepare_patient(incoming, aneend=0, discharge=10)
            self.assertEqual(plan.cleaned.exact_duplicates_removed, 1)
            self.assertEqual(plan.cleaned.conflicting_times, (Decimal(1),))
            self.assertEqual(normalize_actual(core.evaluate_patient(plan)), patient_oracle(records, 0, 10))

    def test_fixed_set_under_all_transforms_and_output_one_hundred(self):
        records = [("-1", "99.80"), ("1", "99.99")]
        plan = core.prepare_patient(records, aneend=0, discharge=10)
        transforms = core.fit_transforms(["0.02", "0.1", "0.3", "1.0", "99.80", "99.99"])
        points = tuple(p.index for p in plan.primary_points)
        for name, tf in transforms.versions().items():
            result = core.evaluate_patient(plan, tf)
            self.assertEqual(tuple(p.time for p in result.primary_points), (Decimal(1),))
            self.assertEqual(tuple(p.index for p in plan.primary_points), points)
        rounded = core.evaluate_patient(plan, transforms.versions()["lower_bounded_round_0p2"])
        self.assertEqual(rounded.primary_points[0].value, Decimal(100))
        self.assertEqual(core.lower_bounded_round_0p2("0.02"), Decimal("0.2"))
        self.assertEqual(core.lower_bounded_round_0p2("0.3"), Decimal("0.4"))

    def test_global_study_left_edge_is_callers_filter_contract(self):
        # Numeric core deliberately has no anestart argument. The caller must
        # filter the global a-7d boundary before fitting AND preparing plans.
        a, e = -700000, 0
        records = [(str(a-604800-1), "0.1"), (str(a-604800), "0.9"), ("-1", "1.0"), ("1", "1.3")]
        filtered = [(t,v) for t,v in records if F(a-604800) <= F(t) <= 10]
        plan = core.prepare_patient(filtered, aneend=e, discharge=10)
        self.assertTrue(all(F(m.time) >= a-604800 for m in plan.cleaned.measurements))
        self.assertEqual(normalize_actual(core.evaluate_patient(plan)), patient_oracle(filtered, e, 10))

    def test_short_decimal_context_exactness(self):
        with localcontext() as context:
            context.prec = 3
            records = [("-1", "1.00000000000000000000000000000"),
                       ("1", "1.29999999999999999999999999999")]
            actual = core.evaluate_patient(core.prepare_patient(records, aneend=0, discharge=10))
            self.assertFalse(actual.primary_label)
            self.assertEqual(normalize_actual(actual), patient_oracle(records, 0, 10))
            self.assertEqual(core.lower_bounded_round_0p2("99.99"), Decimal("100.0"))

    def test_exhaustive_small_tied_reference_pools(self):
        for n in range(1, 8):
            for raw in itertools.combinations_with_replacement(("0.10", "0.30", "1.00", "1.50"), n):
                ref = [F(v) for v in raw]
                fit = core.fit_transforms(raw)
                COUNTS["small_tie_reference_pools"] += 1
                for query in ("0.02", "0.10", "0.20", "0.30", "0.40", "0.99", "1.00", "1.25", "1.50", "2.00"):
                    for name, tf in fit.versions().items():
                        self.assertEqual(F(tf(query)), mapping_oracle(ref, F(query), name), (raw, query, name))
                        COUNTS["generated_mapping_queries"] += 1

    def test_generated_quantile_endpoints_monotonicity_and_ties(self):
        rng = random.Random(2026092101)
        for _ in range(100):
            ref = [rng.choice([2, 10, 20, 30, 50, 70, 90, 100, 110, 130, 150, 200, 400]) for _ in range(rng.randint(41, 200))]
            raw = list(map(cstr, ref))
            expected_ref = list(map(F, raw))
            fit = core.fit_transforms(raw)
            queries = sorted(set(list(range(1, 501, 7)) + ref))
            for name, tf in fit.versions().items():
                actual = [F(tf(cstr(z))) for z in queries]
                expected = [mapping_oracle(expected_ref, F(cstr(z)), name) for z in queries]
                self.assertEqual(actual, expected, name)
                self.assertEqual(actual, sorted(actual), name)
                self.assertTrue(all(v > 0 for v in actual))
                COUNTS["generated_mapping_queries"] += len(queries)
            self.assertEqual(fit.tail_plus_bins(fit.q025), fit.q025)
            self.assertEqual(fit.tail_plus_bins(fit.q975), fit.q975)

    def test_generated_patient_versions_against_independent_oracle(self):
        rng = random.Random(2026092102)
        time_choices = [-604801,-604800,-604799,-200000,-172800,-172799,-100,-1,0,1,2,100,172800,172801,604799,604800,604801]
        value_choices = [2,10,20,29,30,31,40,70,80,99,100,109,110,120,129,130,131,149,150,151,200,230,300,400,9999]
        for _ in range(300):
            records = [(str(rng.choice(time_choices)), cstr(rng.choice(value_choices))) for _ in range(rng.randint(1, 35))]
            rng.shuffle(records)
            discharge = rng.choice([0,1,100,172800,604800,700000])
            plan = core.prepare_patient(records, aneend=0, discharge=discharge)
            ref_raw = [v for _,v in records]
            fit = core.fit_transforms(ref_raw)
            ref_fractions = list(map(F, ref_raw))
            for name, tf in fit.versions().items():
                expected = patient_oracle(records, 0, discharge, lambda x: mapping_oracle(ref_fractions, x, name))
                actual = core.evaluate_patient(plan, tf)
                self.assertEqual(normalize_actual(actual), expected, (name, records, discharge))
                COUNTS["generated_patient_versions"] += 1

    def test_invalid_numerical_inputs_raise(self):
        for value in (True, False, "NaN", "Infinity", "0", "-1"):
            with self.assertRaises(ValueError):
                core.fit_transforms([value])
            with self.assertRaises(ValueError):
                core.prepare_patient([(1, value)], aneend=0, discharge=10)
        with self.assertRaises(ValueError):
            core.fit_transforms([])
        with self.assertRaises(ValueError):
            core.prepare_patient([], aneend=10, discharge=9)


if __name__ == "__main__":
    unittest.main(verbosity=2)

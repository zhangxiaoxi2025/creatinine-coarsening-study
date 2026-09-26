"""Purely synthetic independent checks for paired_resampling. No patient files."""
from decimal import Decimal
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import random
import time
import unittest

import numpy as np
import aki_core as core
import paired_resampling as fast

ROOT = Path(__file__).resolve().parents[1]
CHECKS = {"weighted_core_comparisons": 0, "monotone_minimum_checks": 0}


def cents_decimal(c):
    return Decimal((0, tuple(map(int, str(int(c)))), -2))


def exact_cents(d):
    f = Fraction(d) * 100
    if f.denominator != 1:
        raise ValueError("synthetic helper refuses to round")
    return int(f)


def design_from_records(patients, mode="composite"):
    plans = [core.prepare_patient(records, aneend=0, discharge=604800) for records in patients]
    rp, rv, pp, cv, m48, m7 = [], [], [], [], [], []
    for i, plan in enumerate(plans):
        rp.extend([i] * len(plan.cleaned.measurements))
        rv.extend(exact_cents(m.value) for m in plan.cleaned.measurements)
        points = plan.primary_points if mode == "composite" else plan.secondary_points
        for point in points:
            pp.append(i)
            cv.append(exact_cents(plan.cleaned.measurements[point.index].value))
            m7.append(min(exact_cents(plan.cleaned.measurements[k].value) for k in point.prior7))
            if mode == "composite":
                m48.append(min(exact_cents(plan.cleaned.measurements[k].value) for k in point.prior48))
    design = fast.PairedDesign(n_patients=len(patients), reference_patient=rp, reference_cents=rv,
                              point_patient=pp, current_cents=cv, min7_cents=m7,
                              min48_cents=m48 if mode == "composite" else None, mode=mode)
    return design, plans


def core_weighted(plans, weights, mode):
    reference = [m.value for i,p in enumerate(plans) for _ in range(int(weights[i])) for m in p.cleaned.measurements]
    transforms = core.fit_transforms(reference)
    labels = {}
    for name, transform in transforms.versions().items():
        result = [core.evaluate_patient(p, transform) for p in plans]
        labels[name] = [r.primary_label if mode == "composite" else r.secondary_ratio_label for r in result]
    counts = []
    for name in fast.VERSIONS:
        cells = [0,0,0,0]
        for i,(old,new) in enumerate(zip(labels["identity"], labels[name])):
            if old is None or new is None:
                raise AssertionError("unexpected empty synthetic observable set")
            cells[2*int(old)+int(new)] += int(weights[i])
        counts.append(cells)
    return np.array(counts), labels, transforms


def constant_design(n=4):
    return fast.PairedDesign(n_patients=n, reference_patient=np.repeat(np.arange(n), 2),
                            reference_cents=np.full(2*n, 100), point_patient=np.arange(n),
                            current_cents=np.full(n,100), min48_cents=np.full(n,100),
                            min7_cents=np.full(n,100))


class PairedResamplingTests(unittest.TestCase):
    def assert_core_match(self, design, plans, weights):
        observed = design.evaluate(weights)
        expected, expected_labels, transforms = core_weighted(plans, weights, design.mode)
        np.testing.assert_array_equal(observed.counts, expected)
        for name in fast.VERSIONS:
            np.testing.assert_array_equal(observed.labels[name], expected_labels[name])
            expected_map = [exact_cents(transforms.versions()[name](cents_decimal(c))) for c in design.support]
            np.testing.assert_array_equal(observed.mapping_fit.mappings[name], expected_map)
            CHECKS["weighted_core_comparisons"] += 1

    def test_handcrafted_weighted_tables_and_mapping_exactness(self):
        patients = [
            [(-3600,"1.0"),(1,"1.3"),(2,"1.1")],
            [(-3600,"0.2"),(1,"0.3"),(2,"0.1")],
            [(-3600,"99.8"),(1,"99.99")],
            [(-3600,"0.02"),(1,"0.1"),(2,"0.31")],
        ]
        design, plans = design_from_records(patients)
        for weights in ([1,1,1,1],[4,0,0,0],[0,2,0,2],[0,0,1,3],[1,2,3,4]):
            self.assert_core_match(design, plans, weights)
        out = design.evaluate().tables()
        self.assertEqual(out["identity"]["discordance_count"], 0)
        self.assertEqual(out["identity"]["net_change"], 0)
        self.assertEqual(out["identity"]["n"],4)
        self.assertIn(10000, design.evaluate().mapping_fit.mappings["lower_bounded_round_0p2"])

    def test_generated_weighted_composite_and_ratio_match_core(self):
        rng = random.Random(2026092103)
        choices = [2,10,20,29,30,31,50,60,70,99,100,110,129,130,131,150,200,300,9999]
        nrng = np.random.default_rng(2026092104)
        for mode in ("composite", "ratio"):
            for _ in range(35):
                n = rng.randint(2,8)
                patients = []
                for i in range(n):
                    count = rng.randint(3,15)
                    baseline_t = -3600 if mode == "composite" else -200000
                    times = [baseline_t] + list(range(1, count))
                    patients.append([(t, cents_decimal(rng.choice(choices))) for t in times])
                design, plans = design_from_records(patients, mode)
                for _ in range(3):
                    weights = nrng.multinomial(n, np.full(n,1/n))
                    self.assert_core_match(design,plans,weights)

    def test_whole_patient_reference_multiplicity_not_lab_sampling(self):
        patients = [[(-10,"1"),(1,"2")],
                    [(-10,"3"),(1,"4"),(2,"5"),(3,"6")]]
        design, _ = design_from_records(patients)
        np.testing.assert_array_equal(design.reference_frequency([2,0]), [2,2,0,0,0,0])
        np.testing.assert_array_equal(design.reference_frequency([0,2]), [0,0,2,2,2,2])
        self.assertEqual(int(design.reference_frequency([2,0]).sum()),4)
        self.assertEqual(int(design.reference_frequency([0,2]).sum()),8)
        # Both resamples still have exactly two patients in every paired table.
        self.assertTrue(np.all(design.evaluate([0,2]).counts.sum(axis=1)==2))

    def test_explicit_copied_patients_equal_frequency_weights(self):
        patients = [[(-10,"0.2"),(1,"0.3")], [(-10,"1"),(1,"1.3"),(2,"1.6")],
                    [(-10,"3"),(1,"3.2"),(2,"2.8")]]
        original, _ = design_from_records(patients)
        for weights in ([3,0,0], [0,1,2], [1,2,3]):
            expanded = [records for records,count in zip(patients,weights) for _ in range(count)]
            copied, _ = design_from_records(expanded)
            original_result = original.evaluate(weights)
            copied_result = copied.evaluate()
            np.testing.assert_array_equal(original_result.counts, copied_result.counts)
            for name in fast.VERSIONS:
                repeated_labels = np.repeat(original_result.labels[name], weights)
                np.testing.assert_array_equal(repeated_labels, copied_result.labels[name])

    def test_zero_weight_support_stays_positive_and_exact(self):
        patients = [[(-10,"0.02"),(1,"99.99")], [(-10,"1.0"),(1,"1.0")]]
        design, plans = design_from_records(patients)
        self.assert_core_match(design, plans, [0,2])
        fit = design.fit_mappings([0,2])
        self.assertEqual(fit.q025,100)
        self.assertEqual(fit.q975,100)
        self.assertTrue(np.all(fit.mappings["tail_plus_bins"]==100))
        self.assertTrue(all(np.all(m>0) for m in fit.mappings.values()))

    def test_mapped_minimum_commutes_with_monotone_transforms(self):
        rng = random.Random(2026092105)
        for _ in range(70):
            reference = [rng.randint(1,300) for _ in range(rng.randint(1,80))]
            mapping = core.fit_transforms([cents_decimal(c) for c in reference])
            for _ in range(15):
                window = [rng.choice(reference) for _ in range(rng.randint(1,15))]
                for name, tf in mapping.versions().items():
                    self.assertEqual(tf(cents_decimal(min(window))), min(tf(cents_decimal(c)) for c in window))
                    CHECKS["monotone_minimum_checks"] += 1

    def test_all_small_count_compositions_exact_quantile_ties(self):
        patients = [[(-1,"0.2"),(1,"0.3")], [(-1,"0.7"),(1,"1")],
                    [(-1,"1.1"),(1,"1.3")], [(-1,"1.5"),(1,"2")]]
        design, plans = design_from_records(patients)
        for weights in itertools.product(range(4), repeat=4):
            if sum(weights)==4:
                self.assert_core_match(design, plans, weights)

    def test_single_patient_default_generator_and_degenerate_intervals(self):
        d=constant_design(1)
        result=fast.run_bootstrap(d,rng=fast.center_rngs()["MOVER"],n_resamples=40)
        self.assertEqual(result.replicate_counts.shape,(40,4,4))
        for name in fast.VERSIONS:
            self.assertEqual(result.intervals[name]["discordance"],(0.0,0.0))
            self.assertTrue(result.diagnostics["versions"][name]["zero_degenerate_warning"])
        self.assertEqual(result.diagnostics["distinct_tail_and_bin_mapping_combinations"],1)

    def test_bootstrap_reproducible_pcg64_children_and_denominators(self):
        patients = [[(-10,str(v/100)), (1,str((v+30)/100))] for v in [2,10,20,30,70,100,150,300]]
        d,_=design_from_records(patients)
        first=fast.center_rngs(); second=fast.center_rngs()
        self.assertEqual(tuple(first), fast.CENTER_ORDER)
        one=fast.run_bootstrap(d,rng=first["MOVER"],n_resamples=79)
        two=fast.run_bootstrap(d,rng=second["MOVER"],n_resamples=79)
        np.testing.assert_array_equal(one.replicate_counts,two.replicate_counts)
        self.assertEqual(one.initial_rng_state,two.initial_rng_state)
        self.assertEqual(one.intervals,two.intervals)
        self.assertNotEqual(first["MOVER"].bit_generator.state,first["VitalDB"].bit_generator.state)
        self.assertEqual(one.diagnostics["bit_generator"],"PCG64")
        self.assertTrue(np.all(one.replicate_counts.sum(axis=-1)==8))
        self.assertTrue(np.all(one.replicate_counts[:,0,1:3]==0))
        self.assertEqual(one.diagnostics["failed_replicates"],0)

    def test_bootstrap_refits_mapping_and_preserves_whole_patients(self):
        patients = [[(-10,cents_decimal(v)), (1,cents_decimal(v+30)),(2,cents_decimal(v+50))]
                    for v in [10,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,390]]
        d,_=design_from_records(patients)
        rng=fast.center_rngs()["VitalDB"]
        observed=fast.run_bootstrap(d,rng=rng,n_resamples=53)
        replay=fast.center_rngs()["VitalDB"]
        expected=[]
        for _ in range(53):
            weights=fast.draw_patient_multiplicities(d.n_patients,replay)
            expected.append(d.evaluate(weights).counts)
        np.testing.assert_array_equal(observed.replicate_counts,np.asarray(expected))
        self.assertGreater(observed.diagnostics["distinct_quantile_node_combinations"],1)
        self.assertGreater(observed.diagnostics["distinct_tail_and_bin_mapping_combinations"],1)

    def test_type1_interval_ranks_and_directional_counts(self):
        data=np.arange(40,dtype=float)
        self.assertEqual(fast._type1_central_95(data),(0.0,38.0))
        self.assertEqual(fast._type1_central_95(np.array([7.])),(7.0,7.0))
        rates=fast._metric_rates(np.array([5,3,1,1]))
        self.assertEqual(rates["discordance"],0.4)
        self.assertEqual(rates["net_change"],0.2)
        self.assertEqual(rates["upward"],0.3)
        self.assertEqual(rates["downward"],0.1)

    def test_fixed_inputs_are_copied_not_changed_by_resampling(self):
        rp=np.array([0,0,1,1]); rv=np.array([100,130,200,230])
        d=fast.PairedDesign(n_patients=2,reference_patient=rp,reference_cents=rv,
                           point_patient=[0,1],current_cents=[130,230],min48_cents=[100,200],min7_cents=[100,200])
        before=d.evaluate().counts.copy()
        rv[:]=1; rp[:]=0
        d.evaluate([2,0])
        np.testing.assert_array_equal(before,d.evaluate().counts)
        self.assertFalse(d.support.flags.writeable)
        self.assertFalse(d.reference_patient.flags.writeable)

    def test_unordered_point_patients_correct_or_aggregation(self):
        d=fast.PairedDesign(n_patients=2,reference_patient=[0,0,0,1,1,1],
                           reference_cents=[100,110,130,200,210,230],point_patient=[1,0,1,0],
                           current_cents=[210,130,230,110],min48_cents=[200,100,200,100],min7_cents=[200,100,200,100])
        np.testing.assert_array_equal(d.evaluate().labels["identity"],[True,True])
        self.assertTrue(np.all(d.evaluate().counts.sum(axis=1)==2))

    def test_invalid_input_and_weight_contracts(self):
        base=dict(n_patients=1,reference_patient=[0,0],reference_cents=[100,130],
                  point_patient=[0],current_cents=[130],min48_cents=[100],min7_cents=[100])
        for key,bad in [("n_patients",True),("reference_cents",[100.0,130.0]),
                        ("reference_cents",[0,130]),("current_cents",[131]),
                        ("reference_patient",[0,1]),("min7_cents",[130]),
                        ("point_patient",[True])]:
            args={**base,key:bad}
            with self.assertRaises(ValueError,msg=key): fast.PairedDesign(**args)
        with self.assertRaises(ValueError): fast.PairedDesign(**{**base,"n_patients":2})
        d=fast.PairedDesign(**base)
        for bad in ([0],[-1],[1.0],[True],[1,0],[2**53]):
            with self.assertRaises(ValueError): d.evaluate(bad)
        with self.assertRaises(ValueError): fast.run_bootstrap(d,rng=None)
        with self.assertRaises(ValueError): fast.run_bootstrap(d,rng=fast.center_rngs()["MOVER"],n_resamples=0)
        with self.assertRaises(ValueError): fast.PairedDesign(**{**base,"mode":"ratio"})


if __name__ == "__main__":
    unittest.main(verbosity=2)

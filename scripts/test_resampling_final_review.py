"""Independent final review checks; reads specification only, never patients."""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import math
import unittest
from copy import deepcopy
import numpy as np
import aki_core as core
import paired_resampling as fast

ROOT = Path(__file__).resolve().parents[1]

def tiny_design(mode="composite"):
    kwargs=dict(n_patients=2, reference_patient=[0,0,1,1], reference_cents=[100,130,200,230],
                point_patient=[0,1], current_cents=[130,230], min7_cents=[100,200], mode=mode)
    if mode == "composite": kwargs["min48_cents"]=[100,200]
    return fast.PairedDesign(**kwargs)

class FinalReviewTests(unittest.TestCase):
    def test_frozen_hierarchy_exactly_matches_specification(self):
        spec=json.loads((ROOT/'analysis-spec-v1.0.json').read_text())
        self.assertEqual(tuple(spec['centers']),fast.STUDY_CENTER_ORDER)
        self.assertEqual(tuple(spec['analysis_order']),fast.ANALYSIS_ORDER)
        actual=fast.study_rngs()
        center_seeds=np.random.SeedSequence(spec['resampling']['seed']).spawn(2)
        states=[]
        for center,seed in zip(spec['centers'],center_seeds):
            self.assertEqual(tuple(actual[center]),tuple(spec['analysis_order']))
            for analysis,child in zip(spec['analysis_order'],seed.spawn(8)):
                expected=np.random.Generator(np.random.PCG64(child))
                self.assertEqual(actual[center][analysis].bit_generator.state,expected.bit_generator.state)
                np.testing.assert_array_equal(actual[center][analysis].integers(0,2**63,size=16),expected.integers(0,2**63,size=16))
                states.append(json.dumps(actual[center][analysis].bit_generator.state,sort_keys=True))
        self.assertEqual(len(set(states)),16)

    def test_analysis_execution_order_does_not_advance_other_streams(self):
        one=fast.study_rngs(); other=fast.study_rngs()
        fast.run_bootstrap(tiny_design(),rng=one['mover']['main'],n_resamples=17)
        self.assertEqual(one['mover']['parallel7'].bit_generator.state,other['mover']['parallel7'].bit_generator.state)
        a=fast.run_bootstrap(tiny_design('ratio'),rng=one['mover']['parallel7'],n_resamples=17)
        b=fast.run_bootstrap(tiny_design('ratio'),rng=other['mover']['parallel7'],n_resamples=17)
        np.testing.assert_array_equal(a.replicate_counts,b.replicate_counts)
        self.assertEqual(one['vitaldb']['main'].bit_generator.state,other['vitaldb']['main'].bit_generator.state)

    def test_wrong_bit_generators_rejected_without_consuming_state(self):
        for kind in (np.random.MT19937,np.random.Philox,np.random.SFC64,np.random.PCG64DXSM):
            rng=np.random.Generator(kind(3))
            # JSON supports scalar-only generator states inconsistently; compare
            # explicit generated sequences with a fresh counterpart after error.
            with self.assertRaisesRegex(ValueError,'PCG64'):
                fast.run_bootstrap(tiny_design(),rng=rng,n_resamples=3)
            np.testing.assert_array_equal(rng.integers(0,100,size=12),np.random.Generator(kind(3)).integers(0,100,size=12))

    def test_main_absolute_rule_does_not_leak_into_parallel_ratio(self):
        np.testing.assert_array_equal(tiny_design().evaluate().labels['identity'],[True,True])
        np.testing.assert_array_equal(tiny_design('ratio').evaluate().labels['identity'],[False,False])
        # A later observation can be E7 but outside E48; scalar temporal core
        # decides point membership, rather than the array kernel inventing it.
        p=core.prepare_patient([(-1,'1'),(300000,'1.5')],aneend=0,discharge=604800)
        self.assertEqual(len(p.primary_points),0)
        self.assertEqual(len(p.secondary_points),1)
        out=core.evaluate_patient(p)
        self.assertIsNone(out.primary_label)
        self.assertTrue(out.secondary_ratio_label)
        d=fast.PairedDesign(n_patients=1,reference_patient=[0,0],reference_cents=[100,150],
                            point_patient=[0],current_cents=[150],min7_cents=[100],mode='ratio')
        self.assertTrue(d.evaluate().labels['identity'][0])

    def test_type1_95_intervals_across_rank_boundaries(self):
        for b in (1,2,3,39,40,41,79,80,81,4999,5000):
            values=np.arange(b,dtype=float)[::-1]
            expected=tuple(float(math.ceil(Fraction(num,40)*b)-1) for num in (1,39))
            self.assertEqual(fast._type1_central_95(values),expected)

    def test_exact_cent_thresholds_rounding_and_floor(self):
        # 1.29 is below +.3; 1.30 is exactly +.3. 0.29 rounds to .2;
        # .30 rounds to .4; .01 has an explicitly positive .2 floor.
        d=fast.PairedDesign(n_patients=5,reference_patient=[0,0,1,1,2,2,3,3,4,4],
                            reference_cents=[100,129,100,130,29,30,1,10,9998,9999],
                            point_patient=[0,1,2,3,4],current_cents=[129,130,30,10,9999],
                            min7_cents=[100,100,29,1,9998],min48_cents=[100,100,29,1,9998])
        result=d.evaluate()
        np.testing.assert_array_equal(result.labels['identity'],[False,True,False,True,False])
        rounded=dict(zip(d.support,result.mapping_fit.mappings['lower_bounded_round_0p2']))
        self.assertEqual((rounded[1],rounded[29],rounded[30],rounded[9999]),(20,20,40,10000))
        with self.assertRaisesRegex(ValueError,'integer'):
            fast.PairedDesign(n_patients=1,reference_patient=[0,0],reference_cents=[100,100.1],
                              point_patient=[0],current_cents=[100],min7_cents=[100],mode='ratio')


if __name__ == "__main__":
    unittest.main(verbosity=2)

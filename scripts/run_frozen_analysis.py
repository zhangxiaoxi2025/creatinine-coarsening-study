"""Execute frozen study; output only aggregate logs, private patient label files."""
from pathlib import Path
from datetime import datetime,timezone
import sys,json,hashlib,platform,time
import numpy as np,pandas as pd
from portable_io import sha
REPO=Path(__file__).resolve().parents[1]
A=D=None
REPETITIONS=5000
ENFORCE_EXPECTED=True
from paired_resampling_exact import study_rngs,run_bootstrap,VERSIONS
from analysis_adapter import build_analysis

def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def descr(v):
 x=pd.to_numeric(pd.Series(v),errors='coerce').dropna()
 return {'n':len(x),'missing':len(v)-len(x),'min':float(x.min()) if len(x) else None,'q25':float(x.quantile(.25)) if len(x) else None,'median':float(x.median()) if len(x) else None,'q75':float(x.quantile(.75)) if len(x) else None,'max':float(x.max()) if len(x) else None}
def main():
 access={"scope":"portable local execution; permission obtained independently by operator"}
 assert not (A/'results/execution-complete.json').exists(),'Run already complete; do not silently overwrite.'
 for f in (A/'.private').glob('*-labels.parquet'):raise AssertionError('Existing partial result: review recovery before rerun.')
 spec=json.loads((REPO/'analysis-spec-v1.0.json').read_text())
 expected_config=json.loads((REPO/'config/expected-denominators.json').read_text())
 expected={(c,m):n for c,values in expected_config.items() for m,n in values.items()}
 source={'scope':'PORTABLE_SYNTHETIC_ONLY' if not ENFORCE_EXPECTED else 'FORMAL_REPLICATION_PENDING_INDEPENDENT_VERIFICATION','python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scripts':{p.name:sha(p) for p in (REPO/'scripts').glob('*.py')},'started_at':datetime.now(timezone.utc).isoformat(),'frozen_spec_sha256':sha(REPO/'analysis-spec-v1.0.json'),'unit_scale':10**15,'resampling_repetitions':REPETITIONS}
 dump(A/'evidence/execution-start.json',source)
 rngs=study_rngs(spec['resampling']['seed'])
 combined=[];obs=[];mappings=[];outputs=[]
 for center in spec['centers']:
  cohort=pd.read_parquet(D/'.private'/f'{center}-cohort.parquet')
  labs=pd.read_parquet(D/'.private'/f'{center}-measurements.parquet')
  for mode in spec['analysis_order']:
   started=time.monotonic()
   adapted=build_analysis(cohort,labs,analysis=mode)
   design=adapted.design;keys=adapted.patient_keys
   assert len(keys)==design.n_patients
   if ENFORCE_EXPECTED:assert len(keys)==expected[(center,mode)], "Frozen cohort denominator mismatch"
   assert tuple(sorted(keys))==keys and len(set(keys))==len(keys)
   boot=run_bootstrap(design,rng=rngs[center][mode],n_resamples=REPETITIONS)
   tables=boot.point_estimate.tables()
   assert boot.replicate_counts.shape==(REPETITIONS,4,4)
   assert np.all(boot.replicate_counts.sum(axis=-1)==len(keys))
   assert tables['identity']['discordance_count']==0
   private=A/'.private'
   labels=pd.DataFrame({'case_id':keys,**boot.point_estimate.labels})
   f=private/f'{center}-{mode}-labels.parquet';labels.to_parquet(f,index=False);f.chmod(0o600)
   g=private/f'{center}-{mode}-bootstrap-counts.npz';np.savez_compressed(g,replicate_counts=boot.replicate_counts);g.chmod(0o600)
   dump(A/'results'/f'{center}-{mode}-point-estimates.json',tables)
   chosen=cohort.set_index('case_id').loc[list(keys)]
   chosen_labs=labs[labs.case_id.isin(keys)].copy()
   baseline=chosen_labs[chosen_labs.is_pre7].groupby('case_id').value.min()
   points=chosen_labs[chosen_labs.E7 if mode=='parallel7' else chosen_labs.E48]
   if mode=='S2_target_after5min':points=points[(points.time-points.end)>300]
   postop=chosen_labs[chosen_labs.is_post]
   observation={'center':center,'analysis':mode,'N':len(keys),'age_not_topcoded':descr(chosen.loc[~chosen.age_topcoded,'age']),'age_topcoded':int(chosen.age_topcoded.sum()),'sex_recorded_counts':chosen.sex.fillna('MISSING').value_counts(dropna=False).to_dict(),'ASA':descr(chosen.asa),'baseline_minimum_creatinine_mg_dl':descr(baseline.reindex(keys)),'anesthesia_hours':descr((chosen.end-chosen.start)/3600),'reference_measurements_per_patient':descr(chosen_labs.groupby('case_id').size().reindex(keys,fill_value=0)),'target_points_per_patient':descr(points.groupby('case_id').size().reindex(keys,fill_value=0)),'observed_postop_points_outside_E48':int((postop.is_post&~postop.E48).sum()),'patients_with_postop_points_outside_E48':int(postop.loc[~postop.E48,'case_id'].nunique()),'last_observed_postop_days':descr(postop.assign(day=(postop.time-postop.end)/86400).groupby('case_id').day.max().reindex(keys))}
   obs.append(observation)
   fit=boot.point_estimate.mapping_fit
   floor_mask=design.reference_units<design.unit_scale//10
   mapping={'center':center,'analysis':mode,'reference_measurements':int(len(design.reference_units)),'original_support_n':len(design.support),'q025_mg_dl':fit.q025/design.unit_scale,'q975_mg_dl':fit.q975/design.unit_scale,'node_values_mg_dl':(fit.node_values/design.unit_scale).tolist(),'rounding_floor_affected_records':int(floor_mask.sum()),'rounding_floor_affected_patients':int(len(np.unique(design.reference_patient[floor_mask]))),'versions':{}}
   lookup=[]
   for v in VERSIONS:
    vals=fit.mappings[v]
    mapping['versions'][v]={'output_support_n':int(len(np.unique(vals))),'reference_changed_records':int(fit.reference_counts[vals!=design.support].sum())}
   for j,z in enumerate(design.support):
    lookup.append({'unit_scale':design.unit_scale,'raw_units':int(z),'reference_count':int(fit.reference_counts[j]),**{v:int(fit.mappings[v][j]) for v in VERSIONS}})
   pd.DataFrame(lookup).to_csv(A/'results'/f'{center}-{mode}-mapping-lookup-exact-units.csv',index=False)
   mappings.append(mapping)
   out={'center':center,'analysis':mode,'N':len(keys),'reference_measurements':len(design.reference_units),'target_points':len(design.point_patient),'adapter_diagnostics':adapted.diagnostics,'point_estimates':tables,'intervals':boot.intervals,'diagnostics':boot.diagnostics,'mapping_summary':mapping,'rng_initial_state':boot.initial_rng_state,'rng_final_state':boot.final_rng_state,'input_hashes':{f'{center}-cohort.parquet':sha(D/'.private'/f'{center}-cohort.parquet'),f'{center}-measurements.parquet':sha(D/'.private'/f'{center}-measurements.parquet')},'private_output_hashes':{f.name:sha(f),g.name:sha(g)},'elapsed_seconds':time.monotonic()-started}
   dump(A/'results'/f'{center}-{mode}-analysis.json',out)
   for v,t in tables.items():
    row={'center':center,'analysis':mode,'version':v,**t}
    for metric,bounds in boot.intervals[v].items():row[metric+'_lower']=bounds[0];row[metric+'_upper']=bounds[1]
    combined.append(row)
   outputs.append({'center':center,'analysis':mode,'N':len(keys),'D_B':tables['tail_plus_bins']['discordance'],'up_B':tables['tail_plus_bins']['upward'],'down_B':tables['tail_plus_bins']['downward'],'elapsed_seconds':out['elapsed_seconds']})
   print(json.dumps(outputs[-1]),flush=True)
 pd.DataFrame(combined).to_csv(A/'results/all-comparisons.csv',index=False)
 dump(A/'results/patient-observation-summary.json',obs)
 dump(A/'results/mapping-diagnostics.json',mappings)
 dump(A/'results/execution-complete.json',{'completed_at':datetime.now(timezone.utc).isoformat(),'status':'COMPUTED_PENDING_INDEPENDENT_VERIFICATION','modes_completed':len(outputs),'resamples_per_mode':REPETITIONS,'total_resamples':len(outputs)*REPETITIONS,'center_modes':outputs,'original_source_cores_verified_before_run':True,'portable_result_requires_independent_verification':True})
 print('MAIN_RUN_COMPLETE; independent verification pending',flush=True)
def execute(output, *, repetitions=5000, enforce_expected=True):
 global A,D,REPETITIONS,ENFORCE_EXPECTED
 A=D=Path(output);REPETITIONS=int(repetitions);ENFORCE_EXPECTED=bool(enforce_expected)
 if REPETITIONS<1 or (ENFORCE_EXPECTED and REPETITIONS!=5000):
  raise ValueError("Formal replication requires exactly 5000 resamples")
 main()


if __name__ == '__main__':
 from run_pipeline import main as portable_cli
 raise SystemExit(portable_cli())

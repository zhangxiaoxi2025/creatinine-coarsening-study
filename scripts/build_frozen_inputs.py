"""Frozen cohort rules. Portable caller sets P to a fresh local output directory.
Scientific function bodies preserved from the frozen implementation.
"""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
P = None

def dump(path,obj):path.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n')

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def parse_dates(s):
 out=pd.Series(pd.NaT,index=s.index,dtype='datetime64[ns]')
 for fmt in ['%Y-%m-%d %H:%M:%S','%Y-%m-%d %H:%M','%m/%d/%y %H:%M','%m/%d/%Y %H:%M','%m/%d/%Y %H:%M:%S','%m/%d/%y %H:%M:%S']:
  m=out.isna()&s.ne('');out.loc[m]=pd.to_datetime(s.loc[m],format=fmt,errors='coerce')
 return out

def distribution(s):
 s=pd.to_numeric(s,errors='coerce').dropna()
 return {'n':len(s),**{n:float(s.quantile(q)) if len(s) else None for n,q in [('min',0),('q25',.25),('median',.5),('q75',.75),('max',1)]}}

def build_mover():
 raw=pd.read_parquet(P/'.private/mover_operations.parquet'); labs=pd.read_parquet(P/'.private/mover_creatinine.parquet');dx=pd.read_parquet(P/'.private/mover_visit.parquet')
 counts=raw.loc[raw.MRN.ne('')&raw.LOG_ID.ne('')].groupby('MRN').LOG_ID.nunique()
 critical=['MRN','BIRTH_DATE','AN_START_DATETIME','AN_STOP_DATETIME','PRIMARY_ANES_TYPE_NM','PRIMARY_PROCEDURE_NM','HOSP_ADMSN_TIME','HOSP_DISCH_TIME']
 conflicts=raw.groupby('LOG_ID')[critical].nunique(dropna=False).gt(1).any(axis=1);badids=set(conflicts[conflicts].index)
 x=raw[~raw.LOG_ID.isin(badids)&raw.LOG_ID.ne('')&raw.MRN.ne('')].drop_duplicates('LOG_ID').copy()
 for col in ['AN_START_DATETIME','AN_STOP_DATETIME','HOSP_ADMSN_TIME','HOSP_DISCH_TIME']:x[col+'_p']=parse_dates(x[col])
 x['case_id']=x.LOG_ID;x['patient_id']=x.MRN;x['procedure']=x.PRIMARY_PROCEDURE_NM;x['age']=pd.to_numeric(x.BIRTH_DATE,errors='coerce');x['age_topcoded']=x.age.eq(90)
 x['adult']=x.age.ge(18)&x.age.le(90);x['general']=x.PRIMARY_ANES_TYPE_NM.str.strip().eq('General');x['single']=x.MRN.map(counts).eq(1)
 # Seconds relative to anesthesia start; preserve exact source second resolution.
 x['start']=0.;x['end']=(x.AN_STOP_DATETIME_p-x.AN_START_DATETIME_p).dt.total_seconds()
 x['admission']=(x.HOSP_ADMSN_TIME_p-x.AN_START_DATETIME_p).dt.total_seconds();x['discharge']=(x.HOSP_DISCH_TIME_p-x.AN_START_DATETIME_p).dt.total_seconds()
 x['valid_time']=x.AN_START_DATETIME_p.notna()&x.end.gt(0)&np.isfinite(x.end)
 x['valid_admission_discharge']=x.admission.notna()&x.admission.le(0)&x.discharge.notna()&x.discharge.ge(x.end)
 x['sex']=x.SEX;x['asa']=pd.to_numeric(x.ASA_RATING_C,errors='coerce')
 lab=labs.copy();lab['case_id']=lab.LOG_ID;lab['patient_id']=lab.MRN;lab['value_str']=lab['Observation Value'];lab['value']=pd.to_numeric(lab.value_str,errors='coerce');lab['collection']=parse_dates(lab['Collection Datetime'])
 identity=x.set_index('case_id').patient_id;native=lab.case_id.map(identity)
 qa={'raw_operation_rows':len(raw),'unique_caseids':int(raw.LOG_ID.nunique()),'critical_conflicting_caseids':len(badids),'raw_serum_creatinine_rows':len(lab),'unlinked_or_conflicted_lab_rows':int(native.isna().sum()),'identity_mismatch_lab_rows':int((native.notna()&native.ne(lab.patient_id)).sum()),'invalid_value_or_time_rows':int((~np.isfinite(lab.value)|~lab.value.gt(0)|~lab.value.lt(100)|lab.collection.isna()).sum())}
 lab=lab[native.notna()&native.eq(lab.patient_id)&np.isfinite(lab.value)&lab.value.gt(0)&lab.value.lt(100)&lab.collection.notna()].copy()
 starts=lab.case_id.map(x.set_index('case_id').AN_START_DATETIME_p);lab['time']=(lab.collection-starts).dt.total_seconds()
 d=dx.rename(columns={'LOG_ID':'case_id','mrn':'patient_id','diagnosis_code':'code','dx_name':'term'}).copy();dmatch=d.case_id.map(identity)
 qa['visit_diagnosis_rows']=len(d);qa['visit_diagnosis_identity_mismatch_rows']=int((dmatch.notna()&dmatch.ne(d.patient_id)).sum());qa['visit_diagnosis_unlinked_rows']=int(dmatch.isna().sum())
 d=d[dmatch.notna()&dmatch.eq(d.patient_id)][['case_id','patient_id','code','term']].drop_duplicates()
 return x,lab[['case_id','patient_id','time','value','value_str']],d,qa

def build_vital(path):
 raw=pd.read_csv(path/'clinical_data.csv',dtype=str,keep_default_na=False);labs=pd.read_csv(path/'lab_data.csv',dtype=str,keep_default_na=False)
 assert raw.caseid.ne('').all() and not raw.caseid.duplicated().any()
 x=raw.copy();x['case_id']=x.caseid;x['patient_id']=x.subjectid;x['procedure']=x.opname
 counts=x.groupby('patient_id').case_id.nunique();x['single']=x.patient_id.map(counts).eq(1)
 x['age_topcoded']=x.age.eq('>89');x['age']=pd.to_numeric(x.age,errors='coerce');x['adult']=x.age.ge(18)|x.age_topcoded
 x['general']=x.ane_type.eq('General')
 for to,fr in [('start','anestart'),('end','aneend'),('admission','adm'),('discharge','dis')]:x[to]=pd.to_numeric(x[fr],errors='coerce')
 x['valid_time']=np.isfinite(x.start)&np.isfinite(x.end)&x.end.gt(x.start)
 x['valid_admission_discharge']=np.isfinite(x.admission)&x.admission.le(x.start)&np.isfinite(x.discharge)&x.discharge.ge(x.end)
 x['asa']=pd.to_numeric(x.asa,errors='coerce')
 lab=labs[labs.name.eq('cr')].copy();lab['case_id']=lab.caseid;lab['patient_id']=lab.case_id.map(x.set_index('case_id').patient_id);lab['value_str']=lab.result;lab['value']=pd.to_numeric(lab.result,errors='coerce');lab['time']=pd.to_numeric(lab.dt,errors='coerce')
 qa={'raw_operation_rows':len(raw),'unique_caseids':int(raw.caseid.nunique()),'critical_conflicting_caseids':0,'raw_serum_creatinine_rows':len(lab),'unlinked_or_conflicted_lab_rows':int(lab.patient_id.isna().sum()),'identity_mismatch_lab_rows':0,'invalid_value_or_time_rows':int((~np.isfinite(lab.value)|~lab.value.gt(0)|~lab.value.lt(100)|~np.isfinite(lab.time)).sum())}
 lab=lab[lab.patient_id.notna()&np.isfinite(lab.value)&lab.value.gt(0)&lab.value.lt(100)&np.isfinite(lab.time)].copy()
 d=x[['case_id','patient_id','dx']].rename(columns={'dx':'term'}).copy();d['code']='';d=d[['case_id','patient_id','code','term']]
 return x,lab[['case_id','patient_id','time','value','value_str']],d,qa

def run_dataset(ds,x,lab,dx,qa,procmap,dxmap):
 mapping=procmap[procmap.dataset.eq(ds)].copy();assert not mapping.term.duplicated().any()
 assert set(x.procedure)<=set(mapping.term),f'Unmapped procedure in {ds}'
 x=x.merge(mapping[['term','decision','reason']].rename(columns={'term':'procedure','decision':'procedure_decision','reason':'procedure_reason'}),on='procedure',how='left',validate='many_to_one')
 mappingdx=dxmap[dxmap.dataset.eq(ds)].copy();assert not mappingdx.duplicated(['code','term']).any()
 dx=dx.merge(mappingdx[['code','term','decision','reason']],on=['code','term'],how='left',validate='many_to_one')
 bad_dx=dx[dx.decision.isin(['exclude','uncertain'])]
 dxflag=set(bad_dx.case_id);x['documented_renal_exclusion']=x.case_id.isin(dxflag)
 flow=[]
 def keep(label,mask):
  nonlocal x
  old=len(x);x=x.loc[mask].copy();flow.append({'stage':label,'excluded':old-len(x),'remaining':len(x)})
 flow.append({'stage':'nonconflicting_identity_cases','excluded':None,'remaining':len(x)})
 keep('one_operation_in_complete_source',x.single)
 keep('adult_general_anesthesia',x.adult&x.general)
 keep('valid_anesthesia_times',x.valid_time)
 keep('valid_admission_discharge_order',x.valid_admission_discharge)
 beforeproc=x.copy();exproc=x[~x.procedure_decision.eq('include')]
 procreasons=exproc.groupby(['procedure_decision','procedure_reason']).size().rename('case_count').reset_index().to_dict('records')
 keep('noncardiac_without_recorded_excluded_or_uncertain_procedure',x.procedure_decision.eq('include'))
 procedural=x.copy();documented_counts=bad_dx[bad_dx.case_id.isin(x.case_id)].groupby(['decision','reason']).case_id.nunique().rename('case_count').reset_index().to_dict('records')
 keep('without_documented_renal_replacement_transplant_or_ESRD_screen',~x.documented_renal_exclusion)
 # All labs processed for procedural cohort, so the no-diagnosis-exclusion sensitivity is reproducible.
 lab=lab[lab.case_id.isin(procedural.case_id)].merge(procedural[['case_id','start','end','discharge']],on='case_id',validate='many_to_one')
 lab=lab[lab.time.ge(lab.start-604800)&lab.time.le(np.minimum(lab.end+604800,lab.discharge))].copy()
 duplicate_n=int(lab.duplicated(['case_id','time','value']).sum());lab=lab.drop_duplicates(['case_id','time','value']).copy()
 conflicts=lab.groupby(['case_id','time']).value.nunique().gt(1);bad_times=conflicts[conflicts].index;bad_cases=set(bad_times.get_level_values(0))
 lab=lab[~pd.MultiIndex.from_frame(lab[['case_id','time']]).isin(bad_times)].sort_values(['case_id','time']).copy()
 lab['previous']=lab.groupby('case_id').time.shift(1);lab['gap']=lab.time-lab.previous
 lab['is_pre7']=lab.time.lt(lab.start)&lab.time.ge(lab.start-604800)
 lab['is_post']=lab.time.gt(lab.end)&lab.time.le(np.minimum(lab.end+604800,lab.discharge))
 lab['E48']=lab.is_post&lab.gap.gt(0)&lab.gap.le(172800)
 lab['E7']=lab.is_post&lab.gap.gt(0)&lab.gap.le(604800)
 b=set(lab.loc[lab.is_pre7,'case_id']);post=set(lab.loc[lab.is_post,'case_id']);e48=set(lab.loc[lab.E48,'case_id']);e7=set(lab.loc[lab.E7,'case_id'])
 primary=set(x.case_id)&b&e48;secondary=set(x.case_id)&b&e7;core=set(x.case_id)&b&post
 flow.append({'stage':'pre7_and_observed_inhospital_post7','excluded':len(x)-len(core),'remaining':len(core)})
 flow.append({'stage':'E48_nonempty_primary','excluded':len(core)-len(primary),'remaining':len(primary)})
 procedural_primary=set(procedural.case_id)&b&e48
 x=procedural.copy();x['primary']=x.case_id.isin(primary);x['secondary7']=x.case_id.isin(secondary);x['procedural_only_primary']=x.case_id.isin(procedural_primary)
 x['conflict_in_window']=x.case_id.isin(bad_cases)
 extreme=set(lab.loc[(lab.value.lt(.1)|lab.value.gt(30)),'case_id']);x['extreme_in_window']=x.case_id.isin(extreme)
 x['duration_over24h']=(x.end-x.start).gt(86400)
 # Time-boundary sensitivity changes postoperative evaluation points only; all eligible prior/reference records stay fixed.
 e48late=set(lab.loc[lab.E48&(lab.time-lab.end).gt(300),'case_id']);e7late=set(lab.loc[lab.E7&(lab.time-lab.end).gt(300),'case_id'])
 x['primary_after5min']=x.case_id.isin(primary&e48late)
 x['secondary7_after5min']=x.case_id.isin(secondary&e7late)
 lab['primary']=lab.case_id.isin(primary);lab['secondary7']=lab.case_id.isin(secondary);lab['procedural_only_primary']=lab.case_id.isin(procedural_primary)
 # Private inputs retain source numeric strings, independent of float-based quality screening.
 columns=['case_id','patient_id','procedure','age','age_topcoded','sex','asa','start','end','admission','discharge','documented_renal_exclusion','primary','secondary7','procedural_only_primary','conflict_in_window','extreme_in_window','duration_over24h','primary_after5min','secondary7_after5min']
 f=P/'.private'/f'{ds}-cohort.parquet';x[columns].to_parquet(f,index=False);f.chmod(0o600)
 g=P/'.private'/f'{ds}-measurements.parquet';lab.drop(columns=['previous','gap']).to_parquet(g,index=False);g.chmod(0o600)
 pri=x[x.primary];pl=lab[lab.primary];sl=lab[lab.secondary7]
 out={'dataset':ds,'flow':flow,'exclusions_by_procedure':procreasons,'documented_renal_screen_categories_overlapping':documented_counts,'quality':qa|{'duplicate_rows_in_procedural_window_removed':duplicate_n,'conflicting_times_in_procedural_window_removed':len(bad_times)},'denominators':{'procedural_eligible':len(procedural),'after_documented_renal_screen':int((~procedural.documented_renal_exclusion).sum()),'baseline7_post7_coverage':len(core),'primary_E48_patients':len(primary),'parallel_E7_patients':len(secondary),'primary_E48_points':int(pl.E48.sum()),'parallel_E7_points':int(sl.E7.sum()),'primary_reference_measurements':len(pl),'parallel_reference_measurements':len(sl)},'sensitivity_denominators':{'no_diagnosis_screen_primary':len(procedural_primary),'exclude_any_conflicting_time_primary':int((~pri.conflict_in_window).sum()),'exclude_extreme_values_patient_primary':int((~pri.extreme_in_window).sum()),'exclude_duration_over24h_primary':int((~pri.duration_over24h).sum()),'exclude_first5min_eval_primary':int(pri.primary_after5min.sum()),'all_quality_restrictions_and_after5min_primary':int((~pri.conflict_in_window&~pri.extreme_in_window&~pri.duration_over24h&pri.primary_after5min).sum())},'primary_observation':{'reference_measurements_per_patient':distribution(pl.groupby('case_id').size()),'E48_points_per_patient':distribution(pl[pl.E48].groupby('case_id').size()),'last_postop_days':distribution(pl[pl.is_post].assign(days=lambda a:(a.time-a.end)/86400).groupby('case_id').days.max()),'anesthesia_hours':distribution((pri.end-pri.start)/3600),'age_numeric_not_topcoded':distribution(pri.loc[~pri.age_topcoded,'age']),'age_topcoded_n':int(pri.age_topcoded.sum())},'private_input_hashes':{f.name:sha(f),g.name:sha(g)},'scope':'No patient AKI label or coarsening effect computed; only eligibility, numeric data quality and timestamps.'}
 assert len(primary)<=len(core) and primary<=secondary and pl.E48.sum()>=len(primary)
 assert not x.patient_id.duplicated().any() and not lab.duplicated(['case_id','time']).any()
 assert lab.time.ge(lab.start-604800).all() and lab.time.le(np.minimum(lab.end+604800,lab.discharge)).all()
 dump(P/'results'/f'{ds}-frozen-cohort.json',out)
 pd.DataFrame(flow).to_csv(P/'results'/f'{ds}-cohort-flow.csv',index=False)
 print(json.dumps({'dataset':ds,'denominators':out['denominators'],'sensitivities':out['sensitivity_denominators']},ensure_ascii=False),flush=True)
 return out

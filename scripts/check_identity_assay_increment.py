"""Independent reconstruction of group constraints, fitted states and assay increments."""
from pathlib import Path
import os
os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import csv,json,math,hashlib,statistics
import numpy as np
import joblib
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'PUBLIC_RELEASE.json').is_file())
RUN=ROOT/'runs/identity_assay_increment_v1_20261031';ADEQ=ROOT/'runs/model_adequacy_v1_20261031'
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def load(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,rows):
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf8')
CATS=['host_system','probe_oat1','probe_oat3'];NUMS=['probe_oat1_uM','probe_oat3_uM','uptake_time_min','uptake_temperature_C','pH'];FIELDS=CATS+NUMS
def transform(v,k):return math.log10(v) if v is not None and k in NUMS[:3] else v
def assay(ids,raw,e):
    out=[]
    for oid in ids:
        vals=[];r=raw[oid]
        for k in NUMS:
            n=e['numeric'][k];v=transform(r[k],k);v=n['median'] if v is None else v;vals.append((v-n['center'])/n['scale'])
        for k in CATS:vals.extend(float(r[k]==level) for level in e['categorical_levels'][k])
        vals.extend(float(r[k] is None) for k in FIELDS);out.append(vals)
    return np.array(out)/math.sqrt(8)
def build():
    checks=[]
    def check(name,v,detail=''):
        checks.append({'check':name,'passed':bool(v),'detail':str(detail)})
        if not v:raise AssertionError(name+': '+str(detail))
    spec=load(RUN/'analysis_spec_v0_26.json');summary=load(RUN/'summary.json');check('all_inputs_and_current_spec_hashes_match',all(sha(ROOT/r['path_from_project_root'])==r['sha256'] for r in read(RUN/'INPUT_SHA256.csv')))
    by={r['observation_id']:r for r in read(ROOT/'runs/exploratory_baselines_v1_20261030/derived/exploratory_outcomes_49.csv') if r['in_main']=='True'};raw=load(ROOT/'runs/assay_ridge_v1_20261030/raw_assay_metadata_49.json');y={k:float(r['observed_assay_conditional_delta']) for k,r in by.items()}
    names=load(ADEQ/'analysis_spec_v0_23.json')['descriptors'];desc={r['observation_id']:np.array([float(r[k]) for k in names]) for r in read(ADEQ/'derived/descriptor_bank_46.csv')};kd=load(ADEQ/'derived/tanimoto_matrix_46.json');K=np.array(kd['matrix']);index={x:i for i,x in enumerate(kd['ids'])}
    outer=read(RUN/'derived/outer_partitions_17.csv');inner=read(RUN/'derived/inner_partitions_128.csv');mem=read(RUN/'derived/group_membership_75.csv');states=load(RUN/'models/outer_states_187.json');pred=read(RUN/'derived/outer_predictions_975.csv');pm={(r['view'],r['family'],r['observation_id']):r for r in pred};gmap={(r['view'],r['observation_id']):r['group_id'] for r in mem}
    check('two_views17outer128inner187states975predictions',len(outer)==17 and len(inner)==128 and len(states)==187 and len(pred)==975)
    graph={x for x in by if by[x]['reported_graph_qualified']=='True'};check('graph29_view_is_exact_reported_graph_subset',len(graph)==29 and {r['observation_id'] for r in mem if r['view']==spec['views'][1]}==graph)
    for view in spec['views']:
        ids=sorted(r['observation_id'] for r in mem if r['view']==view);parent={x:x for x in ids}
        def find(x):
            while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
            return x
        for i,a in enumerate(ids):
            for b in ids[i+1:]:
                if by[a]['study_key']==by[b]['study_key'] or by[a]['constitution_smiles']==by[b]['constitution_smiles']:parent[find(a)]=find(b)
        assert all((find(a)==find(b))==(gmap[(view,a)]==gmap[(view,b)]) for a in ids for b in ids)
        tests=[x for r in outer if r['view']==view for x in r['test_ids'].split('|')];assert sorted(tests)==ids
        for f in spec['models']:assert sorted(r['observation_id'] for r in pred if r['view']==view and r['family']==f)==ids
    check('connected_groups_recomputed_by_independent_union_find',True)
    for r in outer:
        a=r['fit_ids'].split('|');b=r['test_ids'].split('|');assert set(a).isdisjoint(b)
        for k in ['study_key','constitution_smiles']:assert {by[x][k] for x in a}.isdisjoint({by[x][k] for x in b})
        if r['view']==spec['views'][1]:assert set(a+b)<=graph
    for r in inner:
        a=r['fit_ids'].split('|');b=r['validation_ids'].split('|');test=r['outer_test_ids'].split('|');assert set(a).isdisjoint(b) and set(a+b).isdisjoint(test)
        for k in ['study_key','constitution_smiles']:assert {by[x][k] for x in a}.isdisjoint({by[x][k] for x in b})
        if r['view']==spec['views'][1]:assert set(a+b+test)<=graph
    check('all_inner_outer_source_constitution_and_view_exclusions',True)
    errs=[];svd=[];encoder_diffs=[];foresterr=[];replayrows=[]
    for s in states:
        view=s['view'];f=s['family'];fit=s['fit_ids'];query=s['test_ids'];yy=np.array([y[x] for x in fit]);kind=s['state_kind']
        if kind=='BASELINE':
            if f=='TRAIN_MEAN':p=np.full(len(query),yy.mean())
            elif f=='TRAIN_MEDIAN':p=np.full(len(query),np.median(yy))
            else:
                p=[]
                for x in query:
                    sims=[K[index[x],index[k]] for k in fit];best=max(sims);p.append(statistics.mean(y[k] for k,v in zip(fit,sims) if abs(v-best)<=1e-12))
        elif kind=='LINEAR':
            dt=np.array([desc[x] for x in fit]);dq=np.array([desc[x] for x in query]);mu=dt.mean(0);sd=dt.std(0);sd=np.where(sd>1e-12,sd,1.);da=(dt-mu)/sd/math.sqrt(10);db=(dq-mu)/sd/math.sqrt(10)
            de=s['descriptor_encoder'];assert de['fit_ids']==fit;encoder_diffs.extend(np.abs(mu-np.array(de['center'])).tolist()+np.abs(sd-np.array(de['scale'])).tolist())
            e=s['assay_encoder'];assert sorted(e['fit_observation_ids'])==sorted(fit)
            for k in CATS:assert e['categorical_levels'][k]==sorted({raw[x][k] for x in fit if raw[x][k] is not None})
            for k in NUMS:
                vals=[transform(raw[x][k],k) for x in fit if raw[x][k] is not None];med=statistics.median(vals) if vals else 0.;filled=np.array([transform(raw[x][k],k) if raw[x][k] is not None else med for x in fit]);scale=float(filled.std());scale=scale if scale>1e-12 else 1.
                encoder_diffs.extend([abs(med-e['numeric'][k]['median']),abs(float(filled.mean())-e['numeric'][k]['center']),abs(scale-e['numeric'][k]['scale'])])
            aa,ab=assay(fit,raw,e),assay(query,raw,e);ma=np.array([[float(raw[x][k] is None) for k in FIELDS] for x in fit])/math.sqrt(8);mb=np.array([[float(raw[x][k] is None) for k in FIELDS] for x in query])/math.sqrt(8)
            if f=='DESCRIPTOR_RIDGE':a,b=da,db
            elif f=='DESCRIPTOR_HALF_BLOCK_RIDGE':a,b=da/math.sqrt(2),db/math.sqrt(2)
            elif f=='DESCRIPTOR_MISSING_RIDGE':a,b=np.hstack([da,ma])/math.sqrt(2),np.hstack([db,mb])/math.sqrt(2)
            elif f=='DESCRIPTOR_ASSAY_RIDGE':a,b=np.hstack([da,aa])/math.sqrt(2),np.hstack([db,ab])/math.sqrt(2)
            elif f=='ASSAY_RIDGE_EXPANDED':a,b=aa,ab
            else:a,b=ma,mb
            st=s['model'];p=b@np.array(st['coefficients'])+st['intercept'];ac=a-a.mean(0);u,v,vh=np.linalg.svd(ac,full_matrices=False);coef=vh.T@((v/(v*v+len(fit)*s['config']['lambda']))*(u.T@(yy-yy.mean())));ind=yy.mean()+(b-a.mean(0))@coef;svd.extend(np.abs(ind-p).tolist())
        elif kind=='KERNEL':
            ka=K[np.ix_([index[x] for x in fit],[index[x] for x in fit])];kb=K[np.ix_([index[x] for x in query],[index[x] for x in fit])];col=ka.mean(0);grand=ka.mean();a=ka-ka.mean(1)[:,None]-col+grand;b=kb-kb.mean(1)[:,None]-col+grand;assert np.allclose(col,s['kernel_encoder']['fit_column_mean'],rtol=0,atol=1e-15)
            p=s['target_mean']+b@np.array(s['dual']);ev,q=np.linalg.eigh(a);du=q@((q.T@(yy-yy.mean()))/(ev+len(fit)*s['config']['lambda']));ind=yy.mean()+b@du;svd.extend(np.abs(ind-p).tolist())
        else:
            fp=ROOT/s['model_path'];assert sha(fp)==s['model_sha256'];rf=joblib.load(fp);z=np.array([desc[x] for x in query],dtype=np.float32);p=[]
            for row in z:
                vv=[]
                for tree in rf.estimators_:
                    t=tree.tree_;n=0
                    while t.children_left[n]!=-1:n=t.children_left[n] if row[t.feature[n]]<=t.threshold[n] else t.children_right[n]
                    vv.append(float(t.value[n,0,0]))
                p.append(math.fsum(vv)/len(vv))
            foresterr.extend(np.abs(np.array(p)-rf.predict(z)).tolist())
        for oid,v in zip(query,p):
            err=abs(float(v)-float(pm[(view,f,oid)]['predicted_delta']));errs.append(err);replayrows.append({'view':view,'family':f,'observation_id':oid,'replay_difference':err})
    check('all_fit_only_encoder_statistics_match',max(encoder_diffs)<1e-12,max(encoder_diffs))
    check('825_saved_model_and_baseline_predictions_replayed',len(errs)==825 and max(errs)<1e-10,max(errs))
    check('independent_SVD_kernel_replay',max(svd)<1e-7,max(svd));check('all75_manual_forest_traversals_match',len(foresterr)==75 and max(foresterr)<1e-12,max(foresterr))
    old={(r['family'],r['observation_id']):r for r in read(ADEQ/'derived/outer_predictions_460.csv')};legacy=[]
    for f in ['TRAIN_MEAN','TRAIN_MEDIAN','MORGAN_1NN_TIE_MEAN','DESCRIPTOR_RIDGE','TANIMOTO_KRR','DESCRIPTOR_RF']:
        for oid in sorted(by):legacy.append(abs(float(pm[(spec['views'][0],f,oid)]['predicted_delta'])-float(old[(f,oid)]['predicted_delta'])))
    check('all276_matching_MAIN46_predictions_reproduce_v023',len(legacy)==276 and max(legacy)<1e-12,max(legacy))
    scores=read(RUN/'derived/inner_scores_11264.csv');curves=read(RUN/'derived/objective_curves_1496.csv');selected=read(RUN/'derived/selected_hyperparameters_187.csv');selectors=read(RUN/'derived/selector_choices_34.csv')
    check('11264scores1496curves187choices34selectors',len(scores)==11264 and len(curves)==1496 and len(selected)==187 and len(selectors)==34)
    ix={}
    for r in scores:ix.setdefault((r['view'],r['outer_fold'],r['family'],r['config_key']),[]).append(float(r['inner_MAE']))
    for r in curves:
        vals=ix[(r['view'],r['outer_fold'],r['family'],r['config_key'])];assert len(vals)==int(r['inner_groups']) and abs(statistics.mean(vals)-float(r['equal_inner_component_MAE']))<1e-12
    for s in selected:
        rr=[r for r in curves if all(r[k]==s[k] for k in ['view','outer_fold','family'])];m=min(float(r['equal_inner_component_MAE']) for r in rr);assert abs(float(s['inner_selection_MAE'])-m)<=1e-12;ties=[json.loads(r['config_key']) for r in rr if abs(float(r['equal_inner_component_MAE'])-m)<=1e-12];conf=json.loads(s['config_key'])
        if 'lambda' in conf:assert conf['lambda']==max(x['lambda'] for x in ties)
        elif s['family']=='DESCRIPTOR_RF':assert (conf['min_samples_leaf'],conf['max_features']=='sqrt')==max((x['min_samples_leaf'],x['max_features']=='sqrt') for x in ties)
    check('all_inner_objectives_and_selection_ties_recomputed',True)
    for r in selectors:
        order=spec['chem_selector_tie_order'] if r['selector']=='NESTED_CHEM_FAMILY' else spec['all_selector_tie_order'];ss={x['family']:x for x in selected if x['view']==r['view'] and x['outer_fold']==r['outer_fold']};best=min(float(ss[f]['inner_selection_MAE']) for f in order);family=next(f for f in order if abs(float(ss[f]['inner_selection_MAE'])-best)<=1e-12);assert r['selected_family']==family
        for p in [p for p in pred if p['view']==r['view'] and p['family']==r['selector'] and p['outer_fold']==r['outer_fold']]:assert p['predicted_delta']==pm[(p['view'],family,p['observation_id'])]['predicted_delta']
    check('all150_nested_selector_predictions_follow_inner_only_choices',len([r for r in pred if r['family'].startswith('NESTED_')])==150)
    for r in pred:assert math.isfinite(float(r['predicted_delta'])) and abs(float(r['observed_delta'])-y[r['observation_id']])<1e-12 and abs(abs(float(r['predicted_delta'])-y[r['observation_id']])-float(r['absolute_error']))<1e-12
    for r in read(RUN/'derived/overall_metrics_26.csv'):
        z=[x for x in pred if x['view']==r['view'] and x['family']==r['family']];assert len(z)==int(r['queries']) and abs(statistics.mean(float(x['absolute_error']) for x in z)-float(r['MAE']))<1e-12
    check('all975_responses_errors_and26overall_metrics_recomputed',True)
    for r in read(RUN/'derived/paired_assay_comparisons_10.csv'):
        ids=[x['observation_id'] for x in mem if x['view']==r['view']];dd=[float(pm[(r['view'],r['candidate'],x)]['absolute_error'])-float(pm[(r['view'],r['comparator'],x)]['absolute_error']) for x in ids];assert abs(statistics.mean(dd)-float(r['paired_MAE_difference']))<1e-12
    check('all10_fixed_assay_increment_comparisons_same_query_replay',True)
    for r in read(RUN/'derived/common_graph29_scoring_comparison_26.csv'):
        z=[x for x in pred if x['view']==r['training_view'] and x['family']==r['family'] and x['observation_id'] in graph];assert len(z)==29 and abs(statistics.mean(float(x['absolute_error']) for x in z)-float(r['MAE']))<1e-12
    check('both_training_views_compare_on_identical_graph29_queryIDs',True)
    cache=read(RUN/'reproduction/RF_fit_cache_manifest.csv');seen=set()
    for r in cache:
        calculated=hashlib.sha256((r['fit_ids']+'|'+r['config_key']).encode()).hexdigest();assert calculated==r['cache_key'] and r['cache_hit']==str(calculated in seen);seen.add(calculated)
    check('RF_cache_key_and_hit_manifest_recomputed',len(seen)==summary['fresh_RF_fits'] and sum(r['cache_hit']=='True' for r in cache)==summary['cached_RF_reuses'])
    f=load(ROOT/'runs/independent_validation_setup_v1_20261031/models/model_freeze_manifest.json');check('external_validation_candidate_unchanged',sha(ROOT/f['bundle_path'])==f['bundle_sha256']==summary['frozen_validation_candidate_sha256'])
    check('no_new_biology_external_test_or_calibration',summary['new_measurements']==0 and not summary['independent_external_validation'] and not summary['calibration_performed'])
    write(RUN/'reproduction/verification_checks_numeric.csv',checks);write(RUN/'reproduction/saved_prediction_replay_825.csv',replayrows)
    dump(RUN/'reproduction/independent_replay_summary.json',{'checks_passed':len(checks),'saved_prediction_max_diff':max(errs),'SVD_eigen_max_diff':max(svd),'forest_manual_max_diff':max(foresterr),'legacy276prediction_max_diff':max(legacy),'independent_human_review':False})
    print(json.dumps({'checks_passed':len(checks),'saved_prediction_max_diff':max(errs),'SVD_eigen_max_diff':max(svd),'MAIN46_legacy_max_diff':max(legacy)}))
if __name__=='__main__':build()

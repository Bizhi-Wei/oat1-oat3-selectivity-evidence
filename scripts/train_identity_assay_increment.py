"""Nested graph29 retraining and controlled descriptor/assay increment comparisons."""
from pathlib import Path
from collections import defaultdict
import os
os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import csv,json,hashlib,math,statistics,sys,platform
import numpy as np
from sklearn.ensemble import RandomForestRegressor
import sklearn,joblib
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'PUBLIC_RELEASE.json').is_file())
sys.path.insert(0,str(ROOT/'scripts'))
from train_assay_ridge import fit_encoder,encode,path_predict,ridge_fit
from train_model_adequacy import desc_encode,kernel_blocks,baseline,cfgkey
RUN=ROOT/'runs/identity_assay_increment_v1_20261031';ADEQ=ROOT/'runs/model_adequacy_v1_20261031';PREV=ROOT/'runs/prediction_robustness_v1_20261031'
BASES=['TRAIN_MEAN','TRAIN_MEDIAN','MORGAN_1NN_TIE_MEAN']
LINEAR=['DESCRIPTOR_RIDGE','DESCRIPTOR_HALF_BLOCK_RIDGE','DESCRIPTOR_MISSING_RIDGE','DESCRIPTOR_ASSAY_RIDGE','ASSAY_RIDGE_EXPANDED','MISSINGNESS_RIDGE_EXPANDED']
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(p,rows):
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
def load(p):return json.loads(p.read_text('utf-8-sig'))
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def groups_for(rows):
    adj=defaultdict(set)
    for r in rows:
        a='S:'+r['study_key'];b='C:'+r['constitution_smiles'];adj[a].add(b);adj[b].add(a)
    assigned={};n=0
    for node in sorted(adj):
        if node in assigned:continue
        n+=1;stack=[node]
        while stack:
            q=stack.pop()
            if q in assigned:continue
            assigned[q]=f'G{n:02d}';stack.extend(adj[q]-assigned.keys())
    return {r['observation_id']:assigned['S:'+r['study_key']] for r in rows}
def matrices(train,query,raw,desc):
    a,b,de=desc_encode(train,query,desc);enc=fit_encoder(train,raw);aa,ab=encode(train,raw,enc),encode(query,raw,enc)
    fields=['host_system','probe_oat1','probe_oat3','probe_oat1_uM','probe_oat3_uM','uptake_time_min','uptake_temperature_C','pH']
    ma=np.array([[float(raw[x][k] is None) for k in fields] for x in train])/math.sqrt(8);mb=np.array([[float(raw[x][k] is None) for k in fields] for x in query])/math.sqrt(8)
    return {'DESCRIPTOR_RIDGE':(a,b),'DESCRIPTOR_HALF_BLOCK_RIDGE':(a/math.sqrt(2),b/math.sqrt(2)),
        'DESCRIPTOR_MISSING_RIDGE':(np.hstack([a,ma])/math.sqrt(2),np.hstack([b,mb])/math.sqrt(2)),
        'DESCRIPTOR_ASSAY_RIDGE':(np.hstack([a,aa])/math.sqrt(2),np.hstack([b,ab])/math.sqrt(2)),
        'ASSAY_RIDGE_EXPANDED':(aa,ab),'MISSINGNESS_RIDGE_EXPANDED':(ma,mb)},enc,de
def choose(options,family):
    best=min(r['score'] for r in options);ties=[r for r in options if abs(r['score']-best)<=1e-12]
    if family in LINEAR+['TANIMOTO_KRR']:return max(ties,key=lambda r:r['config']['lambda'])
    if family=='DESCRIPTOR_RF':return max(ties,key=lambda r:(r['config']['min_samples_leaf'],r['config']['max_features']=='sqrt'))
    return ties[0]
def build():
    for d in ['derived','models','reproduction','figures']:(RUN/d).mkdir(parents=True,exist_ok=True)
    spec=load(RUN/'analysis_spec_v0_26.json');families=spec['models'];trainable=BASES+LINEAR+['TANIMOTO_KRR','DESCRIPTOR_RF'];assert families==trainable+['NESTED_CHEM_FAMILY','NESTED_ALL_FAMILY']
    inputs=[RUN/'analysis_spec_v0_26.json',PREV/'ARTIFACT_SHA256.csv',ROOT/'runs/exploratory_baselines_v1_20261030/derived/exploratory_outcomes_49.csv',ADEQ/'derived/descriptor_bank_46.csv',ADEQ/'derived/tanimoto_matrix_46.json',
        ROOT/'runs/assay_ridge_v1_20261030/raw_assay_metadata_49.json',ADEQ/'derived/outer_partition_manifest.csv',ADEQ/'derived/outer_predictions_460.csv',ROOT/'runs/bridge_structure_review_v1_20261030/derived/cumulative_observations_90_structure_overlay.csv',
        ROOT/'runs/independent_validation_setup_v1_20261031/models/frozen_candidate.json',ROOT/'scripts/train_assay_ridge.py',ROOT/'scripts/train_model_adequacy.py',ROOT/'scripts/train_identity_assay_increment.py']
    # The public release verifies its own input manifest rather than unrelated historical publication artifacts.
    write(RUN/'INPUT_SHA256.csv',[{'path_from_project_root':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in inputs])
    rows=[r for r in read(inputs[2]) if r['in_main']=='True'];by={r['observation_id']:r for r in rows};allids=sorted(by);ledger={r['observation_id']:r for r in read(inputs[8])};assert len(rows)==46
    assert all(r['reported_graph_qualified']==ledger[r['observation_id']]['reported_graph_qualified_review'] for r in rows)
    views={'FROZEN_MAIN46':rows,'REPORTED_GRAPH29_RETRAINED':[r for r in rows if r['reported_graph_qualified']=='True']};assert len(views[spec['views'][1]])==29
    oldparts=read(inputs[6]);oldgroup={x:r['fold_id'] for r in oldparts for x in r['test_ids'].split('|')}
    names=load(ADEQ/'analysis_spec_v0_23.json')['descriptors'];desc={r['observation_id']:np.array([float(r[k]) for k in names]) for r in read(inputs[3])};kd=load(inputs[4]);K=np.array(kd['matrix']);index={x:i for i,x in enumerate(kd['ids'])};rawall=load(inputs[5]);raw={x:rawall[x] for x in allids}
    membership=[];outer=[];inner=[];viewgroups={};viewinfo=[]
    for view,rr in views.items():
        ids=sorted(r['observation_id'] for r in rr);g=groups_for(rr);viewgroups[view]=g
        if view=='FROZEN_MAIN46':assert all((g[a]==g[b])==(oldgroup[a]==oldgroup[b]) for a in ids for b in ids)
        for x in ids:membership.append({'view':view,'observation_id':x,'group_id':g[x],'study_key':by[x]['study_key'],'constitution_smiles':by[x]['constitution_smiles'],'reported_graph_qualified':by[x]['reported_graph_qualified']})
        for group in sorted(set(g.values())):
            fit=[x for x in ids if g[x]!=group];test=[x for x in ids if g[x]==group];outer.append({'view':view,'fold_id':group,'fit_ids':'|'.join(fit),'test_ids':'|'.join(test),'fit_rows':len(fit),'test_rows':len(test),'fit_groups':len(set(g[x] for x in fit))})
            for h in sorted({g[x] for x in fit}):
                ft=[x for x in fit if g[x]!=h];val=[x for x in fit if g[x]==h]
                for field in ['study_key','constitution_smiles']:assert {by[x][field] for x in ft}.isdisjoint({by[x][field] for x in val})
                inner.append({'view':view,'outer_fold':group,'inner_group':h,'fit_ids':'|'.join(ft),'validation_ids':'|'.join(val),'outer_test_ids':'|'.join(test)})
        viewinfo.append({'view':view,'rows':len(ids),'papers':len({by[x]['study_key'] for x in ids}),'constitutions':len({by[x]['constitution_smiles'] for x in ids}),'outer_components':len(set(g.values()))})
    write(RUN/'derived/group_membership_75.csv',membership);write(RUN/'derived/outer_partitions_17.csv',outer);write(RUN/'derived/inner_partitions_128.csv',inner);write(RUN/'derived/view_summary_2.csv',viewinfo)
    dump(RUN/'reproduction/design_lock_before_fits.json',{'spec_sha256':sha(inputs[0]),'outer_partitions_sha256':sha(RUN/'derived/outer_partitions_17.csv'),'inner_partitions_sha256':sha(RUN/'derived/inner_partitions_128.csv'),'eligible_old_input_only':True,'prior_main46_outcomes_already_seen':True,'current_new_fits_started':False})
    y={x:float(by[x]['observed_assay_conditional_delta']) for x in allids};rfconfigs=[{'min_samples_leaf':l,'max_features':f} for l in spec['rf']['min_samples_leaf'] for f in spec['rf']['max_features']]
    rf_cache={};cachelog=[];fresh_rf=0;reused_rf=0;logical_fits=0;linear_kernel_fits=0
    def forest(fit,config,phase):
        nonlocal fresh_rf,reused_rf
        key=hashlib.sha256(('|'.join(fit)+'|'+cfgkey(config)).encode()).hexdigest()
        if key in rf_cache:reused_rf+=1;cached=True
        else:
            rf=RandomForestRegressor(n_estimators=256,criterion='squared_error',bootstrap=True,max_depth=None,random_state=20261031,n_jobs=1,**config);rf.fit(np.array([desc[x] for x in fit]),np.array([y[x] for x in fit]));v=rf.predict(np.array([desc[x] for x in allids]));rf_cache[key]=(rf,v);fresh_rf+=1;cached=False
        cachelog.append({'phase':phase,'cache_key':key,'fit_ids':'|'.join(fit),'config_key':cfgkey(config),'cache_hit':cached,'fit_reads_only_listed_training_responses':True})
        return key,rf_cache[key]
    scores=[];curves=[];selections=[];selectors=[];pred=[];states=[]
    for o in outer:
        view=o['view'];fold=o['fold_id'];train=o['fit_ids'].split('|');test=o['test_ids'].split('|');phase=view+'::'+fold;local=[]
        for ip in [r for r in inner if r['view']==view and r['outer_fold']==fold]:
            fit=ip['fit_ids'].split('|');val=ip['validation_ids'].split('|');yy=np.array([y[x] for x in fit]);vv=np.array([y[x] for x in val])
            def add(f,conf,p):
                r={'view':view,'outer_fold':fold,'inner_group':ip['inner_group'],'family':f,'config_key':cfgkey(conf),'inner_MAE':float(np.abs(p-vv).mean()),'fit_rows':len(fit),'validation_rows':len(val)};scores.append(r);local.append(r)
            for f,p in baseline(fit,val,y,K,index).items():add(f,{},p)
            mats,_,_=matrices(fit,val,raw,desc)
            for f,(a,b) in mats.items():
                for lam,p in path_predict(a,yy,b,spec['lambda_grid']).items():add(f,{'lambda':lam},p);logical_fits+=1;linear_kernel_fits+=1
            a,b,_=kernel_blocks(fit,val,K,index);ym=yy.mean()
            for lam in spec['lambda_grid']:add('TANIMOTO_KRR',{'lambda':lam},ym+b@np.linalg.solve(a+len(fit)*lam*np.eye(len(fit)),yy-ym));logical_fits+=1;linear_kernel_fits+=1
            for conf in rfconfigs:
                _,(rf,pp)=forest(fit,conf,phase+'::INNER::'+ip['inner_group']);add('DESCRIPTOR_RF',conf,np.array([pp[allids.index(x)] for x in val]));logical_fits+=1
        chosen={}
        for f in trainable:
            options=[]
            for ck in sorted({r['config_key'] for r in local if r['family']==f}):
                rr=[r for r in local if r['family']==f and r['config_key']==ck];score=statistics.mean(r['inner_MAE'] for r in rr);options.append({'family':f,'config':json.loads(ck),'score':score});curves.append({'view':view,'outer_fold':fold,'family':f,'config_key':ck,'equal_inner_component_MAE':score,'inner_groups':len(rr)})
            chosen[f]=choose(options,f);selections.append({'view':view,'outer_fold':fold,'family':f,'config_key':cfgkey(chosen[f]['config']),'inner_selection_MAE':chosen[f]['score']})
        sels={}
        for sel,order in [('NESTED_CHEM_FAMILY',spec['chem_selector_tie_order']),('NESTED_ALL_FAMILY',spec['all_selector_tie_order'])]:
            best=min(chosen[f]['score'] for f in order);f=next(f for f in order if abs(chosen[f]['score']-best)<=1e-12);sels[sel]=f;selectors.append({'view':view,'outer_fold':fold,'selector':sel,'selected_family':f,'config_key':cfgkey(chosen[f]['config']),'inner_selection_MAE':chosen[f]['score'],'uses_outer_response':False})
        yy=np.array([y[x] for x in train]);ym=float(yy.mean());pp=baseline(train,test,y,K,index);mats,enc,de=matrices(train,test,raw,desc)
        for f in BASES:states.append({'view':view,'outer_fold':fold,'family':f,'state_kind':'BASELINE','fit_ids':train,'test_ids':test,'config':{}})
        for f,(a,b) in mats.items():
            lam=chosen[f]['config']['lambda'];st=ridge_fit(a,yy,lam);pp[f]=b@np.array(st['coefficients'])+st['intercept'];eig=np.maximum(np.linalg.eigvalsh((a-a.mean(0))@(a-a.mean(0)).T),0)
            states.append({'view':view,'outer_fold':fold,'family':f,'state_kind':'LINEAR','fit_ids':train,'test_ids':test,'config':chosen[f]['config'],'model':st,'assay_encoder':enc,'descriptor_encoder':de,'effective_df_including_intercept':1+float(np.sum(eig/(eig+len(train)*lam)))})
            logical_fits+=1;linear_kernel_fits+=1
        a,b,ke=kernel_blocks(train,test,K,index);lam=chosen['TANIMOTO_KRR']['config']['lambda'];dual=np.linalg.solve(a+len(train)*lam*np.eye(len(train)),yy-ym);pp['TANIMOTO_KRR']=ym+b@dual;logical_fits+=1;linear_kernel_fits+=1
        states.append({'view':view,'outer_fold':fold,'family':'TANIMOTO_KRR','state_kind':'KERNEL','fit_ids':train,'test_ids':test,'config':chosen['TANIMOTO_KRR']['config'],'dual':dual.tolist(),'target_mean':ym,'kernel_encoder':ke})
        key,(rf,rp)=forest(train,chosen['DESCRIPTOR_RF']['config'],phase+'::OUTER');pp['DESCRIPTOR_RF']=np.array([rp[allids.index(x)] for x in test]);logical_fits+=1;fp=RUN/'models'/(view+'_'+fold+'_forest.joblib');joblib.dump(rf,fp)
        states.append({'view':view,'outer_fold':fold,'family':'DESCRIPTOR_RF','state_kind':'FOREST','fit_ids':train,'test_ids':test,'config':chosen['DESCRIPTOR_RF']['config'],'model_path':fp.relative_to(ROOT).as_posix(),'model_sha256':sha(fp),'cache_key':key})
        for sel,f in sels.items():pp[sel]=pp[f].copy()
        for f,p in pp.items():
            selected=sels.get(f,f)
            for x,v in zip(test,p):pred.append({'view':view,'outer_fold':fold,'family':f,'observation_id':x,'study_key':by[x]['study_key'],'compound':by[x]['compound'],'observed_delta':y[x],'predicted_delta':float(v),'absolute_error':abs(float(v)-y[x]),'train_mean_prediction':ym,'selected_family':selected,'config_key':cfgkey(chosen[selected]['config']),'unseen_categorical_fields':'|'.join(k for k in ['host_system','probe_oat1','probe_oat3'] if raw[x][k] is not None and raw[x][k] not in enc['categorical_levels'][k])})
        print(view,fold,'complete','chem',sels['NESTED_CHEM_FAMILY'],'all',sels['NESTED_ALL_FAMILY'],'RFfresh/cache',fresh_rf,reused_rf,flush=True)
    write(RUN/'derived/inner_scores_11264.csv',scores);write(RUN/'derived/objective_curves_1496.csv',curves);write(RUN/'derived/selected_hyperparameters_187.csv',selections);write(RUN/'derived/selector_choices_34.csv',selectors)
    write(RUN/'derived/outer_predictions_975.csv',pred);dump(RUN/'models/outer_states_187.json',states);write(RUN/'reproduction/RF_fit_cache_manifest.csv',cachelog)
    metrics=[];folds=[];papers=[];increments=[];pairdetails=[]
    for view in spec['views']:
        v=[r for r in pred if r['view']==view];pm={(r['family'],r['observation_id']):r for r in v};qids=sorted({r['observation_id'] for r in v});g=viewgroups[view]
        for f in families:
            rr=[r for r in v if r['family']==f];fr=[];pr=[]
            for group in sorted(set(g.values())):
                z=[r for r in rr if r['outer_fold']==group];mae=statistics.mean(r['absolute_error'] for r in z);bm=statistics.mean(pm[('TRAIN_MEAN',r['observation_id'])]['absolute_error'] for r in z);row={'view':view,'family':f,'outer_fold':group,'rows':len(z),'MAE':mae,'matched_mean_MAE':bm,'difference':mae-bm};folds.append(row);fr.append(row)
            for study in sorted({r['study_key'] for r in rr}):
                z=[r for r in rr if r['study_key']==study];vmae=statistics.mean(r['absolute_error'] for r in z);pr.append(vmae);papers.append({'view':view,'family':f,'study_key':study,'rows':len(z),'MAE':vmae})
            metrics.append({'view':view,'family':f,'queries':len(rr),'components':len(fr),'papers':len(pr),'MAE':statistics.mean(r['absolute_error'] for r in rr),'equal_component_MAE':statistics.mean(r['MAE'] for r in fr),'equal_paper_MAE':statistics.mean(pr),'matched_mean_MAE':statistics.mean(pm[('TRAIN_MEAN',r['observation_id'])]['absolute_error'] for r in rr),'components_better_than_mean':sum(r['difference']<-1e-12 for r in fr),'components_worse_than_mean':sum(r['difference']>1e-12 for r in fr)})
        for a,b in spec['paired_comparisons']:
            diffs={x:pm[(a,x)]['absolute_error']-pm[(b,x)]['absolute_error'] for x in qids};gd={h:statistics.mean(diffs[x] for x in qids if g[x]==h) for h in set(g.values())};sd={st:statistics.mean(diffs[x] for x in qids if by[x]['study_key']==st) for st in {by[x]['study_key'] for x in qids}}
            increments.append({'view':view,'candidate':a,'comparator':b,'queries':len(qids),'paired_MAE_difference':statistics.mean(diffs.values()),'equal_component_difference':statistics.mean(gd.values()),'equal_paper_difference':statistics.mean(sd.values()),'components_better':sum(v<-1e-12 for v in gd.values()),'components_worse':sum(v>1e-12 for v in gd.values()),'best_component_difference':min(gd.values()),'worst_component_difference':max(gd.values())})
            pairdetails.extend({'view':view,'candidate':a,'comparator':b,'observation_id':x,'outer_fold':g[x],'study_key':by[x]['study_key'],'paired_error_difference':diffs[x]} for x in qids)
    write(RUN/'derived/overall_metrics_26.csv',metrics);write(RUN/'derived/fold_metrics_221.csv',folds);write(RUN/'derived/paper_metrics_260.csv',papers);write(RUN/'derived/paired_assay_comparisons_10.csv',increments);write(RUN/'derived/paired_query_differences_375.csv',pairdetails)
    graphids=sorted(x for x in allids if by[x]['reported_graph_qualified']=='True');common=[]
    for view in spec['views']:
        for f in families:
            z=[r for r in pred if r['view']==view and r['family']==f and r['observation_id'] in graphids];base={r['observation_id']:r for r in pred if r['view']==view and r['family']=='TRAIN_MEAN'}
            common.append({'training_view':view,'scoring_subset':'COMMON_REPORTED_GRAPH29','family':f,'queries':len(z),'MAE':statistics.mean(r['absolute_error'] for r in z),'matched_mean_MAE':statistics.mean(base[r['observation_id']]['absolute_error'] for r in z),'interpretation':'TRAINING_COMPOSITION_AND_GROUPING_CHANGE_NOT_CAUSAL_IDENTITY_EFFECT'})
    write(RUN/'derived/common_graph29_scoring_comparison_26.csv',common)
    summary={'version':'0.26','date':'2026-10-31','status':'COMPUTED_PENDING_INDEPENDENT_REPLAY','views':viewinfo,'models':13,'outer_predictions':len(pred),'outer_partitions':len(outer),'inner_partitions':len(inner),'inner_score_rows':len(scores),'logical_model_configuration_evaluations':logical_fits,'linear_kernel_fit_evaluations':linear_kernel_fits,'fresh_RF_fits':fresh_rf,'cached_RF_reuses':reused_rf,'actual_configuration_fits':linear_kernel_fits+fresh_rf,'new_measurements':0,'independent_external_validation':False,'calibration_performed':False,'frozen_validation_candidate_sha256':sha(inputs[9]),'materials_status':'ACQUISITION_CLOSED_NO_EXTERNAL_VALIDATION','python':platform.python_version(),'numpy':np.__version__,'sklearn':sklearn.__version__,'metrics':metrics,'paired_comparisons':increments}
    dump(RUN/'summary.json',summary);print(json.dumps({'views':viewinfo,'predictions':len(pred),'logical_fits':logical_fits,'fresh_fits':linear_kernel_fits+fresh_rf,'metrics':metrics,'paired':increments}),flush=True)
if __name__=='__main__':build()

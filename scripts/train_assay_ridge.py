"""Nested grouped ridge comparisons of structure, assay facts, joint and missingness."""
from pathlib import Path
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import csv,hashlib,json,math,statistics,platform,sys
import numpy as np
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'PUBLIC_RELEASE.json').is_file())
sys.path.insert(0,str(ROOT/'scripts'))
from build_assay_features import RUN,PREV,CATS,NUMS,FIELDS,build as metadata_build,read,write,truth,sha
DESIGN=ROOT/'runs/validation_design_v1_20261030'
SPEC=RUN/'analysis_spec_v0_3.json'
MAIN='MAIN_46_MIXED_IDENTITY_TIERS'
PRIMARY='STUDY_AND_CONSTITUTION'
MODELS=['STRUCTURE_RIDGE','ASSAY_RIDGE','STRUCTURE_ASSAY_RIDGE','MISSINGNESS_RIDGE']

def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf8')
def transformed(v,key):return math.log10(v) if key in ['probe_oat1_uM','probe_oat3_uM','uptake_time_min'] and v is not None else v
def fit_encoder(ids,raw):
    enc={'fit_observation_ids':sorted(ids),'categorical_levels':{},'numeric':{}}
    for key in CATS:enc['categorical_levels'][key]=sorted({raw[x][key] for x in ids if raw[x][key] is not None})
    for key in NUMS:
        vals=[transformed(raw[x][key],key) for x in ids if raw[x][key] is not None]
        med=statistics.median(vals) if vals else 0.
        filled=np.array([transformed(raw[x][key],key) if raw[x][key] is not None else med for x in ids],float)
        sd=float(filled.std());enc['numeric'][key]={'median':med,'center':float(filled.mean()),'scale':sd if sd>1e-12 else 1.,'all_missing':not bool(vals)}
    enc['feature_names']=[key+'_z' for key in NUMS]+[key+'='+v for key in CATS for v in enc['categorical_levels'][key]]+[key+'_missing' for key in FIELDS]
    return enc

def encode(ids,raw,enc):
    out=[]
    for oid in ids:
        row=raw[oid];v=[]
        for key in NUMS:
            n=enc['numeric'][key];value=transformed(row[key],key) if row[key] is not None else n['median']
            v.append((value-n['center'])/n['scale'])
        for key in CATS:v.extend(float(row[key]==level) for level in enc['categorical_levels'][key])
        v.extend(float(row[key] is None) for key in FIELDS)
        out.append(v)
    return np.array(out,float)/math.sqrt(len(FIELDS))

def blocks(trainids,testids,raw,struct,enc=None):
    enc=enc if enc is not None else fit_encoder(trainids,raw)
    sa=np.array([struct[x] for x in trainids]);sb=np.array([struct[x] for x in testids])
    aa,ab=encode(trainids,raw,enc),encode(testids,raw,enc)
    ma=np.array([[float(raw[x][k] is None) for k in FIELDS] for x in trainids])/math.sqrt(8)
    mb=np.array([[float(raw[x][k] is None) for k in FIELDS] for x in testids])/math.sqrt(8)
    return {'STRUCTURE_RIDGE':(sa,sb),'ASSAY_RIDGE':(aa,ab),
        'STRUCTURE_ASSAY_RIDGE':(np.hstack([sa,aa])/math.sqrt(2),np.hstack([sb,ab])/math.sqrt(2)),
        'MISSINGNESS_RIDGE':(ma,mb)},enc

def path_predict(x,y,z,lambdas):
    xm=x.mean(axis=0);ym=float(y.mean());xc=x-xm;yc=y-ym
    gram=xc@xc.T;cross=(z-xm)@xc.T
    return {lam:ym+cross@np.linalg.solve(gram+len(y)*lam*np.eye(len(y)),yc) for lam in lambdas}

def ridge_fit(x,y,lam):
    xm=x.mean(axis=0);ym=float(y.mean());xc=x-xm;yc=y-ym
    dual=np.linalg.solve(xc@xc.T+len(y)*lam*np.eye(len(y)),yc)
    coef=xc.T@dual
    return {'lambda':lam,'alpha_sum_squared':len(y)*lam,'coefficients':coef.tolist(),
        'intercept':float(ym-xm@coef),'target_mean':ym,'fit_rows':len(y),'features':x.shape[1]}

def select_and_fit(trainids,testids,raw,struct,yby,groupby,lambdas,phase):
    groups=sorted({groupby[x] for x in trainids});assert len(groups)>=2
    scores=[];manifests=[]
    for group in groups:
        val=[x for x in trainids if groupby[x]==group];fit=[x for x in trainids if groupby[x]!=group]
        mats,_=blocks(fit,val,raw,struct);yy=np.array([yby[x] for x in fit]);vv=np.array([yby[x] for x in val])
        manifests.append({'phase':phase,'inner_group':group,'fit_ids':'|'.join(fit),'validation_ids':'|'.join(val),'outer_test_ids':'|'.join(testids)})
        for model,(a,b) in mats.items():
            paths=path_predict(a,yy,b,lambdas)
            for lam,pred in paths.items():scores.append({'phase':phase,'model':model,'lambda':lam,
                'inner_group':group,'inner_fit_rows':len(fit),'inner_validation_rows':len(val),'inner_MAE':float(np.abs(pred-vv).mean())})
    agg=[];selected={}
    for model in MODELS:
        means={lam:statistics.mean(r['inner_MAE'] for r in scores if r['model']==model and r['lambda']==lam) for lam in lambdas}
        best=min(means.values());chosen=max(lam for lam in lambdas if abs(means[lam]-best)<=1e-12)
        selected[model]=chosen
        for lam,value in means.items():agg.append({'phase':phase,'model':model,'lambda':lam,'inner_groups':len(groups),'equal_inner_group_MAE':value,'selected':lam==chosen})
    mats,enc=blocks(trainids,testids,raw,struct);yy=np.array([yby[x] for x in trainids]);states={};preds={}
    for model,(a,b) in mats.items():
        state=ridge_fit(a,yy,selected[model]);states[model]=state
        preds[model]=b@np.array(state['coefficients'])+state['intercept']
    return preds,{'encoder':enc,'models':states,'fit_ids':trainids,'test_ids':testids},scores,agg,manifests

def evaluate(pred):
    metrics=[];foldmetrics=[];source=[]
    for scenario,scheme,model in sorted({(r['scenario'],r['scheme'],r['model']) for r in pred}):
        rs=[r for r in pred if (r['scenario'],r['scheme'],r['model'])==(scenario,scheme,model)]
        folds=sorted({r['fold_id'] for r in rs});studies=sorted({r['study_key'] for r in rs})
        fm=[];sm=[]
        for f in folds:
            q=[r for r in rs if r['fold_id']==f];value=statistics.mean(r['absolute_error'] for r in q);fm.append(value)
            foldmetrics.append({'scenario':scenario,'scheme':scheme,'model':model,'fold_id':f,'test_rows':len(q),'fit_rows':q[0]['fit_rows'],
                'MAE':value,'lambda':q[0]['lambda'],'unseen_category_query_rows':sum(bool(r['unseen_categorical_fields']) for r in q)})
        for s in studies:
            q=[r for r in rs if r['study_key']==s];value=statistics.mean(r['absolute_error'] for r in q);sm.append(value)
            source.append({'scenario':scenario,'scheme':scheme,'model':model,'study_key':s,'rows':len(q),'MAE':value})
        metrics.append({'scenario':scenario,'scheme':scheme,'model':model,'test_rows':len(rs),'heldout_groups':len(folds),
            'micro_MAE':statistics.mean(r['absolute_error'] for r in rs),'equal_group_macro_MAE':statistics.mean(fm),
            'equal_study_macro_MAE':statistics.mean(sm),'hyperparameter_tuning':'INNER_GROUP_CV_ONLY','calibrated_coverage_measured':False})
    return metrics,foldmetrics,source

def build():
    for d in ['derived','models','reproduction']:(RUN/d).mkdir(parents=True,exist_ok=True)
    raw,rawtable,coverage=metadata_build()
    spec=json.loads(SPEC.read_text(encoding='utf8'));lambdas=spec['lambda_grid']
    data=read(PREV/'derived/exploratory_outcomes_49.csv');by={r['observation_id']:r for r in data}
    yby={k:float(r['observed_assay_conditional_delta']) for k,r in by.items()}
    bank=json.loads((PREV/'models/feature_bank_46_candidates.json').read_text(encoding='utf8'))['candidates']
    struct={}
    for oid,r in by.items():
        bits=bank[r['candidate_structure_id']]['on_bits'];x=np.zeros(2048);x[bits]=1/math.sqrt(len(bits));struct[oid]=x
    membership=read(DESIGN/'derived/group_membership.csv')
    folds=[r for r in read(DESIGN/'derived/outer_leave_one_group_options.csv') if r['scheme'] in spec['schemes']]
    allpred=[];states=[];allscores=[];allagg=[];allmanifest=[]
    for f in folds:
        scenario,scheme=f['scenario'],f['scheme'];phase=scenario+'::'+scheme+'::'+f['fold_id']
        train=sorted(f['reference_observation_ids'].split('|'));test=sorted(f['test_observation_ids'].split('|'))
        groupby={r['observation_id']:r['group_id'] for r in membership if r['scenario']==scenario and r['scheme']==scheme}
        preds,state,scores,agg,man=select_and_fit(train,test,raw,struct,yby,groupby,lambdas,phase)
        states.append({'scenario':scenario,'scheme':scheme,'fold_id':f['fold_id'],**state})
        allscores.extend(scores);allagg.extend(agg);allmanifest.extend(man)
        for model,values in preds.items():
            for oid,value in zip(test,values):
                q=by[oid];unseen=[k for k in CATS if raw[oid][k] is not None and raw[oid][k] not in state['encoder']['categorical_levels'][k]]
                allpred.append({'scenario':scenario,'scheme':scheme,'fold_id':f['fold_id'],'model':model,'observation_id':oid,
                    'study_key':q['study_key'],'compound':q['compound'],'observed_delta':yby[oid],'predicted_delta':float(value),
                    'absolute_error':abs(float(value)-yby[oid]),'lambda':state['models'][model]['lambda'],
                    'fit_rows':len(train),'test_rows':len(test),'known_assay_fields':sum(v is not None for v in raw[oid].values()),
                    'unseen_categorical_fields':'|'.join(unseen),'reported_graph_qualified':q['reported_graph_qualified'],
                    'public_name_only_added_in_sensitivity':q['public_name_only_added_in_sensitivity'],
                    'fit_ids':'|'.join(train),'calibrated_interval_lower':'','calibrated_interval_upper':''})
    # Full-fit models are separate and never used for outer evaluation.
    ids=sorted(r['observation_id'] for r in data if truth(r['in_main']))
    groupby={r['observation_id']:r['group_id'] for r in membership if r['scenario']==MAIN and r['scheme']==PRIMARY}
    _,full,sc,ag,ma=select_and_fit(ids,ids,raw,struct,yby,groupby,lambdas,'FULLFIT_MAIN46')
    full['test_ids']=[]
    for r in ma:r['outer_test_ids']=''
    allscores.extend(sc);allagg.extend(ag);allmanifest.extend(ma)
    metrics,foldmetrics,source=evaluate(allpred)
    old=read(PREV/'derived/overall_metrics_27.csv')
    joined=[{k:r.get(k,'') for k in metrics[0]} for r in old]+metrics
    oldpred=read(PREV/'derived/outer_predictions_1116.csv')
    oldby={(r['scenario'],r['scheme'],r['model'],r['observation_id']):r for r in oldpred}
    paired=[]
    for r in allpred:
        base=oldby[(r['scenario'],r['scheme'],'TRAIN_MEAN',r['observation_id'])]
        paired.append({'scenario':r['scenario'],'scheme':r['scheme'],'model':r['model'],'observation_id':r['observation_id'],
            'fold_id':r['fold_id'],'absolute_prediction_difference_from_train_mean':abs(r['predicted_delta']-float(base['predicted_delta'])),
            'absolute_error_difference_from_train_mean':r['absolute_error']-float(base['absolute_error'])})
    write(RUN/'derived/outer_predictions_1488.csv',allpred)
    write(RUN/'derived/overall_metrics_36.csv',metrics);write(RUN/'derived/comparison_metrics_63.csv',joined)
    write(RUN/'derived/fold_metrics_312.csv',foldmetrics);write(RUN/'derived/per_study_metrics.csv',source)
    write(RUN/'derived/inner_fold_scores.csv',allscores);write(RUN/'derived/lambda_selection.csv',allagg)
    write(RUN/'derived/inner_partition_manifest.csv',allmanifest)
    write(RUN/'derived/paired_comparison_to_train_mean.csv',paired)
    dump(RUN/'models/outer_states_78.json',{'states':states,'numpy':np.__version__,'spec_sha256':sha(SPEC)})
    dump(RUN/'models/fullfit_MAIN46.json',{'fit_scope':'Exploratory all46 main rows; not used in outer metrics','state':full,'spec_sha256':sha(SPEC)})
    summary={'date':'2026-10-30','models':MODELS,'outer_folds':len(states),'outer_model_fits':len(foldmetrics),
        'prediction_rows':len(allpred),'metric_comparisons':len(metrics),'inner_partition_records':len(allmanifest),
        'inner_grid_score_records':len(allscores),'lambda_candidates':lambdas,'main_metadata_coverage':coverage,
        'main_primary_metrics':[r for r in metrics if r['scenario']==MAIN and r['scheme']==PRIMARY],
        'main_primary_mean_absolute_prediction_change_from_train_mean':{m:statistics.mean(r['absolute_prediction_difference_from_train_mean'] for r in paired if r['scenario']==MAIN and r['scheme']==PRIMARY and r['model']==m) for m in MODELS},
        'python':platform.python_version(),'numpy':np.__version__,'solver':'NumPy dual ridge; n_fit*lambda penalty; unpenalized intercept',
        'nested_group_tuning_executed':True,'assay_metadata_models_fitted':True,'calibration_performed':False,
        'external_validation':False,'binary_labels_created':False,'new_measurements':0,'new_literature_requests':0,
        'prior_raw_values_and_identity_statuses_modified':False,'metadata_source_correction':'PF07328948 empty effective context supplemented only in new feature layer from preserved SI pagesS22/S42'}
    dump(RUN/'summary.json',summary)
    return raw,struct,yby,by,folds,allpred,states,metrics,summary
if __name__=='__main__':
    *_,s=build();print(json.dumps(s,ensure_ascii=False,indent=2))

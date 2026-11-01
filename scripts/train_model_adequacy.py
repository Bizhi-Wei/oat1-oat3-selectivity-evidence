"""Locked exploratory test of representation, regularization and function-class adequacy."""
from pathlib import Path
import os
os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import csv,json,hashlib,math,statistics,sys,platform
import numpy as np
from rdkit import Chem,rdBase
from rdkit.Chem import Descriptors,Lipinski,rdMolDescriptors
from sklearn.ensemble import RandomForestRegressor
import sklearn,scipy,joblib
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'PUBLIC_RELEASE.json').is_file())
sys.path.insert(0,str(ROOT/'scripts'))
from train_assay_ridge import blocks,path_predict,ridge_fit
RUN=ROOT/'runs/model_adequacy_v1_20261031';BASE=ROOT/'runs/exploratory_baselines_v1_20261030';PANEL=ROOT/'runs/panel_sensitivity_v1_20261031'
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(p,rows):
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
def load(p):return json.loads(p.read_text('utf-8-sig'))
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
BASES=['TRAIN_MEAN','TRAIN_MEDIAN','MORGAN_1NN_TIE_MEAN']
LINEAR=['JOINT_RIDGE_ORIGINAL_GRID_REPLAY','JOINT_RIDGE_EXPANDED','STRUCTURE_RIDGE_EXPANDED','DESCRIPTOR_RIDGE']
ALL=BASES+LINEAR+['TANIMOTO_KRR','DESCRIPTOR_RF','NESTED_FAMILY_SELECTION']
ORDER=['TRAIN_MEAN','TRAIN_MEDIAN','MORGAN_1NN_TIE_MEAN','STRUCTURE_RIDGE_EXPANDED','DESCRIPTOR_RIDGE','JOINT_RIDGE_EXPANDED','TANIMOTO_KRR','DESCRIPTOR_RF']
def cfgkey(d):return json.dumps(d,sort_keys=True,separators=(',',':'))
def desc_encode(train,query,desc):
    x=np.array([desc[k] for k in train],float);z=np.array([desc[k] for k in query],float)
    med=np.nanmedian(x,axis=0);x=np.where(np.isfinite(x),x,med);z=np.where(np.isfinite(z),z,med)
    mean=x.mean(0);sd=x.std(0);sd=np.where(sd>1e-12,sd,1.)
    return (x-mean)/sd/math.sqrt(10),(z-mean)/sd/math.sqrt(10),{'fit_ids':train,'median':med.tolist(),'center':mean.tolist(),'scale':sd.tolist()}
def kernel_blocks(train,query,K,index):
    a=K[np.ix_([index[x] for x in train],[index[x] for x in train])];b=K[np.ix_([index[x] for x in query],[index[x] for x in train])]
    col=a.mean(0);grand=float(a.mean());ac=a-a.mean(1)[:,None]-col[None,:]+grand;bc=b-b.mean(1)[:,None]-col[None,:]+grand
    return ac,bc,{'fit_column_mean':col.tolist(),'fit_grand_mean':grand}
def matrix_blocks(train,query,raw,struct,desc):
    m,enc=blocks(train,query,raw,struct);da,db,de=desc_encode(train,query,desc)
    return {'JOINT_RIDGE_ORIGINAL_GRID_REPLAY':m['STRUCTURE_ASSAY_RIDGE'],'JOINT_RIDGE_EXPANDED':m['STRUCTURE_ASSAY_RIDGE'],
            'STRUCTURE_RIDGE_EXPANDED':m['STRUCTURE_RIDGE'],'DESCRIPTOR_RIDGE':(da,db)},enc,de
def baseline(train,query,y,K,index):
    a=np.array([y[k] for k in train]);out={'TRAIN_MEAN':np.full(len(query),a.mean()),'TRAIN_MEDIAN':np.full(len(query),np.median(a))}
    out['MORGAN_1NN_TIE_MEAN']=np.array([np.mean([y[t] for t in train if abs(K[index[q],index[t]]-max(K[index[q],index[x]] for x in train))<=1e-12]) for q in query])
    return out
def choose(rows,family):
    best=min(r['score'] for r in rows);ties=[r for r in rows if abs(r['score']-best)<=1e-12]
    if family in LINEAR+['TANIMOTO_KRR']:return max(ties,key=lambda r:r['config']['lambda'])
    if family=='DESCRIPTOR_RF':return max(ties,key=lambda r:(r['config']['min_samples_leaf'],r['config']['max_features']=='sqrt'))
    return ties[0]
def build():
    for d in ['derived','models','reproduction','figures','metadata']:(RUN/d).mkdir(parents=True,exist_ok=True)
    spec=load(RUN/'analysis_spec_v0_23.json');assert spec['models']==ALL
    inputs=[RUN/'analysis_spec_v0_23.json',BASE/'derived/exploratory_outcomes_49.csv',BASE/'models/feature_bank_46_candidates.json',
        ROOT/'runs/assay_ridge_v1_20261030/raw_assay_metadata_49.json',PANEL/'derived/outer_partition_manifest_54.csv',PANEL/'derived/outer_predictions_2072.csv',
        ROOT/'runs/journal_selection_v1_20261031/ARTIFACT_SHA256.csv',ROOT/'scripts/train_assay_ridge.py',ROOT/'scripts/train_model_adequacy.py']
    # Discover the actual legacy prediction filename without changing any analysis input.
    inputs=[p if p.exists() else PANEL/'derived/outer_predictions_2072.csv' for p in inputs]
    if not inputs[5].exists():
        candidates=list((PANEL/'derived').glob('*prediction*.csv'));print('Legacy candidates',','.join(p.name for p in candidates),flush=True)
        inputs[5]=next(p for p in candidates if len(read(p))==2072)
    write(RUN/'INPUT_SHA256.csv',[{'path_from_project_root':p.relative_to(ROOT).as_posix(),'sha256':sha(p),'bytes':p.stat().st_size} for p in inputs])
    rows=[r for r in read(inputs[1]) if r['in_main']=='True'];by={r['observation_id']:r for r in rows};ids=sorted(by);index={x:i for i,x in enumerate(ids)}
    outer=[r for r in read(inputs[4]) if r['view']=='FROZEN46_CONTROL' and r['feature']=='NONCHIRAL' and r['fullfit']=='False'];assert len(outer)==9
    group={oid:r['fold_id'] for r in outer for oid in r['test_ids'].split('|')};assert set(group)==set(ids)
    bank=load(inputs[2])['candidates'];bits={x:set(bank[by[x]['candidate_structure_id']]['on_bits']) for x in ids};struct={x:np.array([int(j in bits[x]) for j in range(2048)],float)/math.sqrt(len(bits[x])) for x in ids}
    K=np.array([[len(bits[a]&bits[b])/len(bits[a]|bits[b]) for b in ids] for a in ids]);assert np.linalg.eigvalsh(K).min()>-1e-10
    desc={};descriptor_rows=[]
    for oid in ids:
        m=Chem.MolFromSmiles(by[oid]['canonical_isomeric_smiles']);assert m is not None
        ds=[Descriptors.MolWt(m),Descriptors.MolLogP(m),rdMolDescriptors.CalcTPSA(m),Lipinski.NumHDonors(m),Lipinski.NumHAcceptors(m),Lipinski.NumRotatableBonds(m),Lipinski.RingCount(m),rdMolDescriptors.CalcFractionCSP3(m),Chem.GetFormalCharge(m),m.GetNumHeavyAtoms()]
        assert all(math.isfinite(x) for x in ds);desc[oid]=ds;descriptor_rows.append({'observation_id':oid,'candidate_structure_id':by[oid]['candidate_structure_id'],**dict(zip(spec['descriptors'],ds))})
    write(RUN/'derived/descriptor_bank_46.csv',descriptor_rows);dump(RUN/'derived/tanimoto_matrix_46.json',{'ids':ids,'matrix':K.tolist()})
    rawall=load(inputs[3]);raw={k:rawall[k] for k in ids};partitions=[]
    for o in outer:
        train=sorted(o['fit_ids'].split('|'));test=sorted(o['test_ids'].split('|'))
        assert set(train).isdisjoint(test) and set(train)|set(test)==set(ids)
        for g in sorted({group[x] for x in train}):
            fit=[x for x in train if group[x]!=g];val=[x for x in train if group[x]==g]
            for field in ['study_key','constitution_smiles']:assert {by[x][field] for x in fit}.isdisjoint({by[x][field] for x in val})
            partitions.append({'outer_fold':o['fold_id'],'inner_group':g,'fit_ids':'|'.join(fit),'validation_ids':'|'.join(val),'outer_test_ids':'|'.join(test)})
    write(RUN/'derived/inner_partition_manifest.csv',partitions);write(RUN/'derived/outer_partition_manifest.csv',outer)
    dump(RUN/'reproduction/design_lock_before_fits.json',{'spec_sha256':sha(inputs[0]),'descriptor_bank_sha256':sha(RUN/'derived/descriptor_bank_46.csv'),'partition_manifest_sha256':sha(RUN/'derived/inner_partition_manifest.csv'),'prior_targets_already_seen':True,'no_new_fits_yet':True})
    y={x:float(by[x]['observed_assay_conditional_delta']) for x in ids}
    allscores=[];allagg=[];states=[];predictions=[];selectedrows=[];familychoices=[];fitcounter=0
    rfconfigs=[{'min_samples_leaf':l,'max_features':f} for l in spec['rf']['min_samples_leaf'] for f in spec['rf']['max_features']]
    for o in outer:
        fold=o['fold_id'];train=sorted(o['fit_ids'].split('|'));test=sorted(o['test_ids'].split('|'));inners=[r for r in partitions if r['outer_fold']==fold];local=[]
        for ip in inners:
            fit=ip['fit_ids'].split('|');val=ip['validation_ids'].split('|');yy=np.array([y[x] for x in fit]);vv=np.array([y[x] for x in val]);mats,_,_=matrix_blocks(fit,val,raw,struct,desc)
            def add(family,conf,pred):
                r={'outer_fold':fold,'inner_group':ip['inner_group'],'family':family,'config_key':cfgkey(conf),'fit_rows':len(fit),'validation_rows':len(val),'inner_MAE':float(np.abs(pred-vv).mean())};local.append(r);allscores.append(r)
            for family,p in baseline(fit,val,y,K,index).items():add(family,{},p)
            for family,(a,b) in mats.items():
                grid=spec['original_grid'] if family=='JOINT_RIDGE_ORIGINAL_GRID_REPLAY' else spec['expanded_lambda_grid']
                for lam,p in path_predict(a,yy,b,grid).items():add(family,{'lambda':lam},p);fitcounter+=1
            a,b,_=kernel_blocks(fit,val,K,index);ym=float(yy.mean())
            for lam in spec['expanded_lambda_grid']:
                pred=ym+b@np.linalg.solve(a+len(fit)*lam*np.eye(len(fit)),yy-ym);add('TANIMOTO_KRR',{'lambda':lam},pred);fitcounter+=1
            dx=np.array([desc[x] for x in fit]);dz=np.array([desc[x] for x in val])
            for conf in rfconfigs:
                rf=RandomForestRegressor(n_estimators=spec['rf']['n_estimators'],criterion='squared_error',bootstrap=True,max_depth=None,random_state=spec['rf']['random_state'],n_jobs=1,**conf)
                rf.fit(dx,yy);add('DESCRIPTOR_RF',conf,rf.predict(dz));fitcounter+=1
        chosen={}
        for family in ALL[:-1]:
            configs=sorted({r['config_key'] for r in local if r['family']==family});curve=[]
            for ck in configs:
                rr=[r for r in local if r['family']==family and r['config_key']==ck];assert len(rr)==len(inners)
                v=statistics.mean(r['inner_MAE'] for r in rr);curve.append({'family':family,'config':json.loads(ck),'score':v})
                allagg.append({'outer_fold':fold,'family':family,'config_key':ck,'equal_inner_group_MAE':v,'inner_groups':len(rr)})
            selected=choose(curve,family);chosen[family]=selected
            selectedrows.append({'outer_fold':fold,'family':family,'config_key':cfgkey(selected['config']),'equal_inner_group_MAE':selected['score']})
        best=min(chosen[f]['score'] for f in ORDER);family=next(f for f in ORDER if abs(chosen[f]['score']-best)<=1e-12)
        familychoices.append({'outer_fold':fold,'selected_family':family,'config_key':cfgkey(chosen[family]['config']),'equal_inner_group_MAE':chosen[family]['score'],'selection_uses_outer_outcome':False})
        mats,enc,de=matrix_blocks(train,test,raw,struct,desc);yy=np.array([y[x] for x in train]);vv=np.array([y[x] for x in test]);pp=baseline(train,test,y,K,index);statelist=[]
        for f in BASES:statelist.append({'outer_fold':fold,'family':f,'fit_ids':train,'test_ids':test,'config':{},'state_kind':'BASELINE','target_mean':float(yy.mean())})
        for f,(a,b) in mats.items():
            lam=chosen[f]['config']['lambda'];st=ridge_fit(a,yy,lam);pp[f]=b@np.array(st['coefficients'])+st['intercept'];xc=a-a.mean(0);eig=np.linalg.eigvalsh(xc@xc.T)
            statelist.append({'outer_fold':fold,'family':f,'fit_ids':train,'test_ids':test,'config':chosen[f]['config'],'state_kind':'LINEAR','model':st,'assay_encoder':enc if f.startswith('JOINT') else {},'descriptor_encoder':de if f=='DESCRIPTOR_RIDGE' else {},'effective_df_including_intercept':1+float(np.sum(np.maximum(eig,0)/(np.maximum(eig,0)+len(train)*lam)))})
            fitcounter+=1
        a,b,ke=kernel_blocks(train,test,K,index);lam=chosen['TANIMOTO_KRR']['config']['lambda'];ym=float(yy.mean());dual=np.linalg.solve(a+len(train)*lam*np.eye(len(train)),yy-ym);pp['TANIMOTO_KRR']=ym+b@dual;eig=np.linalg.eigvalsh(a)
        statelist.append({'outer_fold':fold,'family':'TANIMOTO_KRR','fit_ids':train,'test_ids':test,'config':chosen['TANIMOTO_KRR']['config'],'state_kind':'KERNEL','dual':dual.tolist(),'target_mean':ym,'kernel_encoder':ke,'effective_df_including_intercept':1+float(np.sum(np.maximum(eig,0)/(np.maximum(eig,0)+len(train)*lam)))})
        fitcounter+=1
        rf=RandomForestRegressor(n_estimators=spec['rf']['n_estimators'],criterion='squared_error',bootstrap=True,max_depth=None,random_state=spec['rf']['random_state'],n_jobs=1,**chosen['DESCRIPTOR_RF']['config']);rf.fit(np.array([desc[x] for x in train]),yy);pp['DESCRIPTOR_RF']=rf.predict(np.array([desc[x] for x in test]));fitcounter+=1
        modelpath=RUN/'models'/(fold+'_descriptor_forest.joblib');joblib.dump(rf,modelpath)
        statelist.append({'outer_fold':fold,'family':'DESCRIPTOR_RF','fit_ids':train,'test_ids':test,'config':chosen['DESCRIPTOR_RF']['config'],'state_kind':'FOREST','model_path':modelpath.relative_to(ROOT).as_posix(),'model_sha256':sha(modelpath),'mean_tree_depth':float(np.mean([t.tree_.max_depth for t in rf.estimators_])),'mean_tree_leaves':float(np.mean([t.tree_.n_leaves for t in rf.estimators_]))})
        pp['NESTED_FAMILY_SELECTION']=pp[family].copy()
        for f,p in pp.items():
            for oid,v,obs in zip(test,p,vv):predictions.append({'outer_fold':fold,'family':f,'observation_id':oid,'study_key':by[oid]['study_key'],'observed_delta':float(obs),'predicted_delta':float(v),'absolute_error':float(abs(v-obs)),'train_mean_prediction':float(yy.mean()),'selected_family_for_pipeline':family if f=='NESTED_FAMILY_SELECTION' else f,'config_key':cfgkey(chosen[family if f=='NESTED_FAMILY_SELECTION' else f]['config'])})
        states.extend(statelist);print('Completed',fold,'inner_groups',len(inners),'selector',family,'fits',fitcounter,flush=True)
    metrics=[];foldmetrics=[];sourcemetrics=[]
    for f in ALL:
        rr=[r for r in predictions if r['family']==f];base={r['observation_id']:r for r in predictions if r['family']=='TRAIN_MEAN'};fm=[];sm=[]
        for g in sorted(group.values()):
            if g in {r['outer_fold'] for r in fm}:continue
            r=[x for x in rr if x['outer_fold']==g];err=statistics.mean(x['absolute_error'] for x in r);be=statistics.mean(base[x['observation_id']]['absolute_error'] for x in r)
            item={'family':f,'outer_fold':g,'rows':len(r),'MAE':err,'mean_baseline_MAE_same_queries':be,'paired_MAE_difference':err-be};fm.append(item);foldmetrics.append(item)
        for s in sorted({r['study_key'] for r in rr}):
            r=[x for x in rr if x['study_key']==s];v=statistics.mean(x['absolute_error'] for x in r);sm.append(v);sourcemetrics.append({'family':f,'study_key':s,'rows':len(r),'MAE':v})
        mae=statistics.mean(r['absolute_error'] for r in rr);bm=statistics.mean(base[r['observation_id']]['absolute_error'] for r in rr)
        metrics.append({'family':f,'queries':len(rr),'outer_groups':len(fm),'micro_MAE':mae,'group_macro_MAE':statistics.mean(r['MAE'] for r in fm),'paper_macro_MAE':statistics.mean(sm),'paired_micro_difference_from_mean':mae-bm,'mean_absolute_prediction_distance_from_mean':statistics.mean(abs(r['predicted_delta']-r['train_mean_prediction']) for r in rr),'groups_better_than_mean':sum(r['paired_MAE_difference']<-1e-12 for r in fm),'groups_worse_than_mean':sum(r['paired_MAE_difference']>1e-12 for r in fm),'max_adverse_group_difference':max(r['paired_MAE_difference'] for r in fm)})
    write(RUN/'derived/inner_scores.csv',allscores);write(RUN/'derived/inner_objective_curves.csv',allagg);write(RUN/'derived/selected_hyperparameters.csv',selectedrows);write(RUN/'derived/nested_family_choices_9.csv',familychoices)
    write(RUN/'derived/outer_predictions_460.csv',predictions);write(RUN/'derived/metrics_10.csv',metrics);write(RUN/'derived/fold_metrics_90.csv',foldmetrics);write(RUN/'derived/paper_metrics_120.csv',sourcemetrics)
    dump(RUN/'models/outer_states_81.json',states)
    dump(RUN/'summary.json',{'version':'0.23','date':'2026-10-31','status':'COMPUTED_PENDING_INDEPENDENT_NUMERIC_VERIFICATION','rows':46,'groups':9,'models':10,'predictions':len(predictions),'inner_partitions':len(partitions),'inner_score_rows':len(allscores),'model_configuration_fits':fitcounter,'new_measurements':0,'model_selection_is_exploratory':True,'independent_final_test':False,'calibration_performed':False,'metrics':metrics,'python':platform.python_version(),'numpy':np.__version__,'rdkit':rdBase.rdkitVersion,'scikit_learn':sklearn.__version__,'scipy':scipy.__version__,'joblib':joblib.__version__})
    print(json.dumps({'computed_models':10,'predictions':len(predictions),'fits':fitcounter,'metrics':metrics},ensure_ascii=False),flush=True)
if __name__=='__main__':build()

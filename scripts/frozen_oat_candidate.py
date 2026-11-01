"""Outcome-free inference for a frozen validation candidate; not an eligibility certificate."""
from pathlib import Path
import os
os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import csv,json,math,hashlib,argparse
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors,Lipinski,rdMolDescriptors,rdFingerprintGenerator
import joblib
CATS=['host_system','probe_oat1','probe_oat3'];NUMS=['probe_oat1_uM','probe_oat3_uM','uptake_time_min','uptake_temperature_C','pH'];FIELDS=CATS+NUMS
FEATURE_COLUMNS=['record_id','canonical_isomeric_smiles']+FIELDS
FPARGS={'radius':2,'fpSize':2048,'includeChirality':False,'countSimulation':False,'useBondTypes':True,'includeRingMembership':True,'onlyNonzeroInvariants':False,'includeRedundantEnvironments':False}
def chemical(smiles):
    mol=Chem.MolFromSmiles(smiles)
    if mol is None:raise ValueError('Invalid graph; do not substitute a name-derived structure silently.')
    bits=list(rdFingerprintGenerator.GetMorganGenerator(**FPARGS).GetFingerprint(mol).GetOnBits());assert bits
    desc=[Descriptors.MolWt(mol),Descriptors.MolLogP(mol),rdMolDescriptors.CalcTPSA(mol),Lipinski.NumHDonors(mol),Lipinski.NumHAcceptors(mol),Lipinski.NumRotatableBonds(mol),Lipinski.RingCount(mol),rdMolDescriptors.CalcFractionCSP3(mol),Chem.GetFormalCharge(mol),mol.GetNumHeavyAtoms()]
    parent=Chem.Mol(mol);Chem.RemoveStereochemistry(parent)
    return {'on_bits':bits,'descriptors':desc,'constitution_smiles':Chem.MolToSmiles(parent,isomericSmiles=False,canonical=True)}
def assay(rows,enc):
    out=[]
    for row in rows:
        x={}
        for k in CATS:x[k]=None if row.get(k) in ['',None] else str(row[k])
        for k in NUMS:
            z=row.get(k);x[k]=None if z in ['',None] else float(z)
            if x[k] is not None and (not math.isfinite(x[k]) or (k in NUMS[:3] and x[k]<=0)):raise ValueError('Invalid explicit assay number: '+k)
        vals=[]
        for k in NUMS:
            e=enc['numeric'][k];z=x[k]
            if z is None:z=e['median']
            elif k in NUMS[:3]:z=math.log10(z)
            vals.append((z-e['center'])/e['scale'])
        for k in CATS:vals.extend(float(x[k]==v) for v in enc['categorical_levels'][k])
        vals.extend(float(x[k] is None) for k in FIELDS);out.append(vals)
    return np.array(out)/math.sqrt(8)
def predict(bundle_path,rows):
    bundle_path=Path(bundle_path);b=json.loads(bundle_path.read_text('utf8'));assert b['role']=='FROZEN_EXTERNAL_VALIDATION_CANDIDATE_NOT_VALIDATED_DEPLOYMENT_MODEL'
    for r in rows:
        unexpected=set(r)-set(FEATURE_COLUMNS)
        if unexpected:raise ValueError('Inference accepts features only; remove outcome/audit fields: '+','.join(sorted(unexpected)))
        if not r.get('record_id') or not r.get('canonical_isomeric_smiles'):raise ValueError('Explicit record ID and reviewed graph are required.')
    cs=[chemical(r['canonical_isomeric_smiles']) for r in rows];family=b['selected_family'];st=b['fitted_state'];n=len(rows)
    if family in ['TRAIN_MEAN','TRAIN_MEDIAN']:p=np.full(n,st['constant'])
    elif family=='MORGAN_1NN_TIE_MEAN':
        p=[]
        for c in cs:
            q=set(c['on_bits']);sims=[len(q&set(x))/len(q|set(x)) for x in b['training_on_bits']];best=max(sims)
            p.append(np.mean([v for v,z in zip(b['training_responses'],sims) if abs(z-best)<=1e-12]))
    elif family in ['STRUCTURE_RIDGE_EXPANDED','JOINT_RIDGE_EXPANDED','DESCRIPTOR_RIDGE']:
        if family=='DESCRIPTOR_RIDGE':
            e=st['descriptor_encoder'];z=np.array([c['descriptors'] for c in cs]);a=(z-np.array(e['center']))/np.array(e['scale'])/math.sqrt(10)
        else:
            a=np.zeros((n,2048))
            for i,c in enumerate(cs):a[i,c['on_bits']]=1/math.sqrt(len(c['on_bits']))
            if family=='JOINT_RIDGE_EXPANDED':a=np.hstack([a,assay(rows,st['assay_encoder'])])/math.sqrt(2)
        p=a@np.array(st['model']['coefficients'])+st['model']['intercept']
    elif family=='TANIMOTO_KRR':
        a=[]
        for c in cs:
            q=set(c['on_bits']);a.append([len(q&set(x))/len(q|set(x)) for x in b['training_on_bits']])
        a=np.array(a);a=a-a.mean(1)[:,None]-np.array(st['kernel_encoder']['fit_column_mean'])+st['kernel_encoder']['fit_grand_mean'];p=st['target_mean']+a@np.array(st['dual'])
    else:
        path=bundle_path.parent/st['model_file'];actual=hashlib.sha256(path.read_bytes()).hexdigest()
        if actual!=st['model_sha256']:raise ValueError('Frozen forest hash mismatch.')
        rf=joblib.load(path);p=rf.predict(np.array([c['descriptors'] for c in cs]))
    return [{'record_id':r['record_id'],'predicted_assay_conditional_delta':float(v),'frozen_train_mean_prediction':b['training_mean'],'selected_family':family,
             'same_training_constitution':c['constitution_smiles'] in b['training_constitutions'],'nearest_training_tanimoto':max(len(set(c['on_bits'])&set(x))/len(set(c['on_bits'])|set(x)) for x in b['training_on_bits']),
             'scope':'POINT_PREDICTION_ONLY_DATA_QUALIFICATION_REQUIRED_NO_INTERVAL'} for r,c,v in zip(rows,cs,p)]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--bundle',required=True);ap.add_argument('--features',required=True);ap.add_argument('--output',required=True);a=ap.parse_args()
    with Path(a.features).open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
    if not rows:raise ValueError('Feature file has no eligible records; no external prediction was made.')
    out=predict(a.bundle,rows)
    with Path(a.output).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,list(out[0]));w.writeheader();w.writerows(out)
    print(json.dumps({'predictions':len(out),'external_validation_performed':False,'note':'Inference alone is not source/endpoint qualification or performance evaluation.'}))
if __name__=='__main__':main()

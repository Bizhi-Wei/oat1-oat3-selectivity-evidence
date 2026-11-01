"""Project explicit assay metadata without outcomes, study-ID features, or guessed fields."""
from pathlib import Path
import csv,json,math,hashlib
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'PUBLIC_RELEASE.json').is_file())
RUN=ROOT/'runs/assay_ridge_v1_20261030'
PREV=ROOT/'runs/exploratory_baselines_v1_20261030'
CATS=['host_system','probe_oat1','probe_oat3']
NUMS=['probe_oat1_uM','probe_oat3_uM','uptake_time_min','uptake_temperature_C','pH']
FIELDS=CATS+NUMS
PF_PDF=ROOT/'runs/expansion_v1_20261021/sources/PMID39560668_SI_27854910_file_50635191.pdf'
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(p,rows):
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
def truth(v):return str(v).lower()=='true'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def normalize_probe(s):
    if not s:return None
    s=str(s).upper()
    if 'CARBOXYFLUORESCEIN' in s:return '6_CF'
    if 'TENOFOVIR' in s:return 'TENOFOVIR'
    if 'ESTRONE' in s:return 'ESTRONE_SULFATE'
    if 'PAH' in s:return 'PAH'
    if s=='FLUORESCEIN':return 'FLUORESCEIN'
    raise ValueError('Unreviewed probe '+s)

def build():
    (RUN/'derived').mkdir(parents=True,exist_ok=True)
    raw={};audit=[];table=[]
    for r in read(PREV/'derived/exploratory_outcomes_49.csv'):
        c=json.loads(r['effective_context_json']);a={k:None for k in FIELDS};basis={k:'UNKNOWN_OR_NOT_CONFIRMED' for k in FIELDS}
        host=c.get('host','')
        if 'S2' in host:a['host_system']='S2_MOUSE'
        elif 'HEK293' in host:a['host_system']='HEK293'
        elif 'proximal tubule epithelial' in host.lower():a['host_system']='CIPTEC'
        if a['host_system']:basis['host_system']='effective_context_json.host; normalized cell system, not culture temperature'
        for t in ['oat1','oat3']:
            key='probe_'+t
            a[key]=normalize_probe(c.get(key,c.get('probe')))
            if a[key]:basis[key]='effective_context_json.'+(key if key in c else 'probe')+'; isotope label removed, chemical probe identity retained'
            for alias in [key+'_uM','probe_concentration_'+t+'_uM','probe_concentration_uM']:
                if alias in c:
                    # Explicit null is not replaced by a screening concentration or a different endpoint.
                    if c[alias] is not None:a[key+'_uM']=float(c[alias]);basis[key+'_uM']='effective_context_json.'+alias
                    break
            if a[key+'_uM'] is None:
                target=next((x for x in c.get('targets',[]) if x.get('target')==t.upper() and x.get('endpoint')=='IC50'),None)
                if target and target.get('probe_concentration_uM') is not None:
                    a[key+'_uM']=float(target['probe_concentration_uM']);basis[key+'_uM']='effective_context_json.targets['+t.upper()+'].probe_concentration_uM; IC50 scoped'
        for key in ['uptake_time_min','uptake_temperature_C','pH']:
            aliases=[key,'temperature_C'] if key=='uptake_temperature_C' else [key]
            for alias in aliases:
                if alias in c:
                    if c[alias] is not None:a[key]=float(c[alias]);basis[key]='effective_context_json.'+alias
                    break
        notes=[]
        if r['observation_id']=='PF07328948_S7_PAIR':
            assert sha(PF_PDF)==r['source_sha256']
            vals={'host_system':'HEK293','probe_oat1':'PAH','probe_oat3':'ESTRONE_SULFATE','probe_oat1_uM':.5,'probe_oat3_uM':.1}
            a.update(vals)
            for key in vals:basis[key]='Existing supplied SI TableS7 pS22; probe concentrations also methods pS42, visually rechecked2026-10-30'
            notes.append('Retain at-Km footnote versus below-Km methods conflict; explicit concentrations agree, no S/Km feature. Uptake1-4min range not resolved per target, no midpoint or temperature/pH imputation as evidence.')
        if r['study_key']=='pmid:27467266':notes.append('Screening concentration1uM not transferred to IC50; both probe concentrations remain missing.')
        if r['study_key'] in ['pmid:30892891','pmid:30108866']:notes.append('Actual uptake temperature unconfirmed; buffer preheat/culture37C not transferred.')
        for key in NUMS:
            if a[key] is not None:assert math.isfinite(a[key])
        for key in ['probe_oat1_uM','probe_oat3_uM','uptake_time_min']:
            if a[key] is not None:assert a[key]>0
        raw[r['observation_id']]=a
        table.append({'observation_id':r['observation_id'],'study_key':r['study_key'],'compound':r['compound'],
            'in_main':r['in_main'],'reported_graph_qualified':r['reported_graph_qualified'],**a,
            'known_fields':sum(v is not None for v in a.values()),'source_file':r['source_file'],'source_sha256':r['source_sha256'],'notes':' '.join(notes)})
        for key in FIELDS:audit.append({'observation_id':r['observation_id'],'study_key':r['study_key'],
            'field':key,'value':a[key],'known':a[key] is not None,'basis':basis[key],
            'source_file':r['source_file'],'source_sha256':r['source_sha256'],'metadata_evidence_upgraded_to_training_qualification':False})
    coverage=[]
    for key in FIELDS:
        main=[r for r in table if truth(r['in_main'])]
        coverage.append({'field':key,'main_rows':len(main),'known_main_rows':sum(r[key] is not None for r in main),
            'known_main_studies':len({r['study_key'] for r in main if r[key] is not None}),
            'missing_main_rows':sum(r[key] is None for r in main),
            'distinct_observed_values':'|'.join(sorted({str(r[key]) for r in main if r[key] is not None}))})
    write(RUN/'derived/raw_assay_metadata_49.csv',table)
    write(RUN/'derived/metadata_field_evidence_392.csv',audit)
    write(RUN/'derived/metadata_coverage.csv',coverage)
    (RUN/'raw_assay_metadata_49.json').write_text(json.dumps(raw,ensure_ascii=False,indent=2),encoding='utf8')
    return raw,table,coverage
if __name__=='__main__':
    raw,table,coverage=build();print(json.dumps(coverage,ensure_ascii=False,indent=2))

"""Static scientific figure of all frozen comparisons and the CA4 numeric decomposition."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'PUBLIC_RELEASE.json').is_file())
RUN=ROOT/'runs/prediction_robustness_v1_20261031'
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def build():
    (RUN/'figures').mkdir(exist_ok=True)
    rows=read(RUN/'derived/stratum_metrics_80.csv');order=list(dict.fromkeys(r['family'] for r in rows));idx={(r['family'],r['stratum']):r for r in rows};cs=read(RUN/'derived/CA4_feature_contributions_10.csv')
    labels=['Train mean','Train median','Morgan 1-NN','Joint ridge: original grid','Joint ridge: expanded grid','Structure ridge: expanded','Descriptor ridge','Tanimoto kernel ridge','Descriptor random forest','Nested family selection']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none'})
    fig,(a,b)=plt.subplots(1,2,figsize=(15,7),gridspec_kw={'width_ratios':[1.18,1]});fig.subplots_adjust(left=.20,right=.96,top=.85,bottom=.26,wspace=.6)
    yy=np.arange(10)
    for j,(stratum,label,color,marker) in enumerate([('ALL46','All 46','#333333','o'),('REPORTED_GRAPH29','Reported graph 29','#2b75a1','s'),('OTHER_IDENTITY17','Other identity 17','#bc642b','^')]):
        vals=[float(idx[(f,stratum)]['paired_MAE_difference']) for f in order];a.scatter(vals,yy+(j-1)*.2,s=34,c=color,marker=marker,label=label,zorder=3)
    a.set_yticks(yy,labels);a.set_ylim(9.6,-.6);a.axvline(0,lw=1,color='.65');a.grid(axis='x',alpha=.15);a.set_xlabel('MAE difference from the mean on the same queries');a.set_title('a  Frozen-prediction identity strata',loc='left',fontweight='bold',pad=14);a.legend(loc='upper center',bbox_to_anchor=(.5,-.16),frameon=False,ncol=1,fontsize=9)
    vals=np.array([float(r['additive_contribution_vs_training_mean']) for r in cs]);names=[r['feature'] for r in cs];y=np.arange(10);b.barh(y,vals,color=np.where(vals<0,'#337da5','#c77642'),height=.66);b.set_yticks(y,names);b.set_ylim(9.6,-.6);b.axvline(0,lw=1,color='.65');b.grid(axis='x',alpha=.15);b.set_xlabel('Additive contribution to prediction minus train mean');b.set_title('b  CA4: saved descriptor ridge, λ = 10⁻⁸',loc='left',fontweight='bold',pad=14)
    for i,v in enumerate(vals):b.text(v+(.015 if v>=0 else -.015),i,f'{v:+.3f}',ha='left' if v>=0 else 'right',va='center',fontsize=8)
    b.set_xlim(-.68,.63);b.text(.5,-.17,'All 10 raw descriptors are within training marginal ranges.\nMean −0.344 → prediction −1.581; observed −0.403.',transform=b.transAxes,ha='center',va='top',fontsize=9)
    for ax in [a,b]:ax.spines[['right','top']].set_visible(False)
    fig.suptitle('Predictive gains vary by evidence stratum; marginal ranges do not guarantee reliability',x=.56,y=.97,fontsize=14)
    fig.text(.20,.035,'Evaluation subsets use the existing mixed-identity training fits. No subset retraining or intervals are shown.\nFeature contributions are correlated, parameterization-specific arithmetic; they do not establish biological mechanisms.',fontsize=9,color='.3')
    for ext in ['png','svg']:fig.savefig(RUN/'figures'/('Figure_prediction_robustness.'+ext),dpi=300,facecolor='white')
    plt.close(fig)
    (RUN/'figures/robustness_figure_scope.json').write_text(json.dumps({'all_families':order,'strata':['ALL46','REPORTED_GRAPH29','OTHER_IDENTITY17'],'new_fits':0,'identity_certification':False,'feature_decomposition_is_causal':False,'intervals_or_significance_claim':False},indent=2)+'\n',encoding='utf8')
    print('Saved PNG/SVG')
if __name__=='__main__':build()

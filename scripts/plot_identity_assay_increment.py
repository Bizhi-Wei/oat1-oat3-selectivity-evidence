"""Python-only figures for the two prespecified retraining and assay comparisons."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'PUBLIC_RELEASE.json').is_file())
RUN=ROOT/'runs/identity_assay_increment_v1_20261031'
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

def main():
    out=RUN/'figures';out.mkdir(exist_ok=True)
    contract={'conclusion':'Chemical-family gains do not persist after graph29 retraining; assay increments depend on the view, weighting and controls.',
      'archetype':'quantitative grid','backend':'Python/matplotlib','dimensions_inches':{'Figure6':[7.2,5.4],'FigureS3':[7.2,3.5]},
      'panels':'Separate46/29 views; all13 procedures in Figure6, all5 fixed comparisons in FigureS3.',
      'n_definition':'46 records/12 papers/9 components and29 records/8 papers/8 components; records are not independent laboratory replicates.',
      'split':'Outer joint source-and-constitution components, fresh graph29 components; inner component tuning only.',
      'metric':'Query-weighted and equal-component MAE; FigureS3 is candidate minus comparator.',
      'intervals':'None; dots are estimands, not uncertainty intervals. No significance tests.',
      'source_data':['derived/overall_metrics_26.csv','derived/paired_assay_comparisons_10.csv'],
      'exports':['editable SVG','PDF','400dpi PNG'],'journal':'Exploratory manuscript draft; no journal-specific final compliance claim.'}
    (out/'identity_assay_figure_contract.json').write_text(json.dumps(contract,indent=2)+'\n',encoding='utf8')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.2,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'legend.frameon':False})
    views=['FROZEN_MAIN46','REPORTED_GRAPH29_RETRAINED']
    titles=['a  Mixed identity: 46 records\n    9 strict components','b  Graph-associated: 29 records\n    8 strict components']
    data=read(RUN/'derived/overall_metrics_26.csv');families=list(dict.fromkeys(r['family'] for r in data));ix={(r['view'],r['family']):r for r in data}
    labels=['Training mean','Training median','Morgan 1-NN','Descriptor ridge (D)','D / sqrt(2): scale control','D + missingness','D + assay','Assay ridge','Missingness ridge','Tanimoto kernel ridge','Descriptor random forest','Inner chemical selection','Inner full-menu selection']
    colors=['#276b92','#b86d32'];markers=['o','D'];metrics=['MAE','equal_component_MAE']
    fig,axs=plt.subplots(1,2,figsize=(7.2,5.4),sharey=True)
    fig.subplots_adjust(left=.29,right=.98,top=.83,bottom=.20,wspace=.24)
    yy=np.arange(len(families))
    for ax,v,title in zip(axs,views,titles):
        for j,(key,c,marker,lab) in enumerate(zip(metrics,colors,markers,['Query-weighted MAE','Equal-component MAE'])):
            ax.scatter([float(ix[(v,f)][key]) for f in families],yy+(j-.5)*.20,s=17,color=c,marker=marker,label=lab,zorder=3)
            ax.axvline(float(ix[(v,'TRAIN_MEAN')][key]),color=c,ls=':',lw=.9,alpha=.8)
        ax.set_xlim(.36,.84);ax.set_xticks([.4,.6,.8]);ax.set_ylim(len(families)-.4,-.6)
        ax.set_yticks(yy,labels);ax.set_title(title,loc='left',fontsize=7.1,fontweight='bold',pad=14)
        ax.set_xlabel('MAE in delta pIC50 units');ax.grid(axis='x',alpha=.15)
    axs[0].legend(loc='upper left',bbox_to_anchor=(-.02,-.12),ncol=2,fontsize=7,handletextpad=.3,columnspacing=1)
    fig.text(.29,.975,'Retraining and assay increments under strict transfer',fontsize=10,weight='bold',va='top')
    fig.text(.29,.035,'Each view uses its own training means and nested component splits.\nDotted lines mark mean-reference MAEs. View errors are not a causal identity comparison.',fontsize=6.8,color='.25')
    for ext in ['png','svg','pdf']:fig.savefig(out/('Figure6.'+ext),dpi=400,facecolor='white')
    plt.close(fig)
    pairs=read(RUN/'derived/paired_assay_comparisons_10.csv')
    pair_order=[('DESCRIPTOR_ASSAY_RIDGE','DESCRIPTOR_RIDGE'),('DESCRIPTOR_ASSAY_RIDGE','DESCRIPTOR_HALF_BLOCK_RIDGE'),('DESCRIPTOR_ASSAY_RIDGE','DESCRIPTOR_MISSING_RIDGE'),('DESCRIPTOR_MISSING_RIDGE','DESCRIPTOR_RIDGE'),('NESTED_ALL_FAMILY','NESTED_CHEM_FAMILY')]
    plabel=['D + assay minus D','D + assay minus scale control','D + assay minus (D + missingness)','D + missingness minus D','Full-menu minus chemical selection']
    pix={(r['view'],r['candidate'],r['comparator']):r for r in pairs}
    fig,axs=plt.subplots(1,2,figsize=(7.2,3.5),sharey=True)
    fig.subplots_adjust(left=.36,right=.98,top=.82,bottom=.31,wspace=.24)
    for ax,v,title in zip(axs,views,['a  Mixed identity: 46','b  Graph-associated: 29']):
        for j,(key,c,marker,lab) in enumerate(zip(['paired_MAE_difference','equal_component_difference'],colors,markers,['Query-weighted','Equal-component'])):
            ax.scatter([float(pix[(v,*p)][key]) for p in pair_order],np.arange(5)+(j-.5)*.20,s=20,c=c,marker=marker,label=lab,zorder=3)
        ax.set_yticks(np.arange(5),plabel);ax.set_ylim(4.55,-.55);ax.axvline(0,c='.5',lw=.8)
        ax.set_xlim(-.065,.065);ax.set_xticks([-.05,0,.05]);ax.grid(axis='x',alpha=.15)
        ax.set_title(title,loc='left',weight='bold',fontsize=8,pad=12);ax.set_xlabel('Paired MAE difference')
    axs[0].legend(loc='upper left',bbox_to_anchor=(-.04,-.24),ncol=2,fontsize=7,handletextpad=.25,columnspacing=.65)
    fig.text(.36,.98,'Assay increments depend on weighting and controls',fontsize=10,weight='bold',va='top')
    fig.text(.36,.035,'Negative values favor the candidate; positive values favor the comparator.\nIdentical queries within each comparison; no intervals or independent testing.',fontsize=6.8,color='.25')
    for ext in ['png','svg','pdf']:fig.savefig(out/('FigureS3.'+ext),dpi=400,facecolor='white')
    plt.close(fig)
    print('Figure6 and FigureS3: PNG/SVG/PDF saved')
if __name__=='__main__':main()

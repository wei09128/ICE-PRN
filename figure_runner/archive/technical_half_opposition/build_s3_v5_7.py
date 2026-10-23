#!/usr/bin/env python3
"""Compare each chemical's reference-opposition scores across MCF7 and PC3."""
from pathlib import Path
from datetime import datetime
import argparse,json,zipfile
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
BREAST='MCF7_vs_MCF10A'
PROSTATE=['PC3_vs_RWPE1','prostate_tumor_vs_adjacent_15pairs']
CATEGORIES=['MCF7_only','PC3_only','both','neither','uncertain']
COLORS={'MCF7_only':'#0072B2','PC3_only':'#D55E00','both':'#009E73','neither':'#666666','uncertain':'#BBBBBB'}
NAMES={'MCF7_only':'MCF7 opposition only','PC3_only':'PC3 opposition only','both':'Opposition in both','neither':'Opposition in neither','uncertain':'Half disagreement / zero / missing'}

def quadrant(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float)
    result=np.full(x.shape,'uncertain',dtype='<U12')
    valid=np.isfinite(x)&np.isfinite(y)&(x!=0)&(y!=0)
    result[valid&(x>0)&(y<0)]='MCF7_only'
    result[valid&(x<0)&(y>0)]='PC3_only'
    result[valid&(x>0)&(y>0)]='both'
    result[valid&(x<0)&(y<0)]='neither'
    return result

def preference(x0,x1,y0,y1):
    d0=np.asarray(x0,float)-np.asarray(y0,float)
    d1=np.asarray(x1,float)-np.asarray(y1,float)
    ans=np.full(d0.shape,'uncertain',dtype='<U12')
    valid=np.isfinite(d0)&np.isfinite(d1)
    ans[valid&(d0>0)&(d1>0)]='MCF7_higher'
    ans[valid&(d0<0)&(d1<0)]='PC3_higher'
    return d0,d1,ans

def compare(scores,out,expected=4185):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    required={'reference','cell','identity_group','reversal_spearman_half0','reversal_spearman_half1','reversal_spearman_halfmean'}
    if not required.issubset(scores):raise ValueError('Missing required score columns')
    if scores.duplicated(['reference','identity_group']).any():raise ValueError('Duplicate identity/reference rows')
    breast=scores[scores.reference.eq(BREAST)].set_index('identity_group')
    if len(breast)!=expected:raise ValueError('Expected identity cohort not found')
    if not breast.cell.eq('MCF7').all():raise ValueError('Breast score cell mismatch')
    details=[];summaries=[];transitions=[]
    for pref in PROSTATE:
        prostate=scores[scores.reference.eq(pref)].set_index('identity_group')
        if len(prostate)!=expected or set(prostate.index)!=set(breast.index):raise ValueError('Paired identity cohort differs')
        if not prostate.cell.eq('PC3').all():raise ValueError('Prostate score cell mismatch')
        prostate=prostate.loc[breast.index]
        x0=breast.reversal_spearman_half0.to_numpy();x1=breast.reversal_spearman_half1.to_numpy()
        y0=prostate.reversal_spearman_half0.to_numpy();y1=prostate.reversal_spearman_half1.to_numpy()
        q0=quadrant(x0,y0);q1=quadrant(x1,y1)
        stable=np.where((q0==q1)&(q0!='uncertain'),q0,'uncertain')
        d0,d1,p=preference(x0,x1,y0,y1)
        for i,name in enumerate(breast.index):
            row=dict(identity_group=name,prostate_reference=pref,breast_reference=BREAST,
                MCF7_score_half0=x0[i],MCF7_score_half1=x1[i],PC3_score_half0=y0[i],PC3_score_half1=y1[i],
                MCF7_score_halfmean=breast.reversal_spearman_halfmean.iloc[i],
                PC3_score_halfmean=prostate.reversal_spearman_halfmean.iloc[i],
                discovery_half0=q0[i],confirmation_half1=q1[i],stable_pattern=stable[i],
                MCF7_minus_PC3_half0=d0[i],MCF7_minus_PC3_half1=d1[i],stable_score_preference=p[i])
            if 'compound_names' in breast:row['compound_names']=breast.compound_names.iloc[i]
            details.append(row)
        for a in CATEGORIES:
            for b in CATEGORIES:
                transitions.append(dict(prostate_reference=pref,half0_pattern=a,half1_pattern=b,n=int(np.sum((q0==a)&(q1==b)))))
        for c in CATEGORIES:
            selected=q0==c
            summaries.append(dict(prostate_reference=pref,discovery_pattern=c,discovery_n=int(selected.sum()),
                same_pattern_in_half1=int(np.sum(selected&(q1==c))) if c!='uncertain' else None,
                confirmation_fraction=float(np.mean(q1[selected]==c)) if c!='uncertain' and selected.any() else None,
                stable_n=int(np.sum(stable==c))))
    detail=pd.DataFrame(details).sort_values(['prostate_reference','identity_group'])
    detail.to_csv(out/'paired_context_scores.tsv',sep='\t',index=False)
    pd.DataFrame(summaries).to_csv(out/'context_group_counts.tsv',sep='\t',index=False)
    pd.DataFrame(transitions).to_csv(out/'half_transition_counts.tsv',sep='\t',index=False)
    plot(detail,out)
    return {'status':'COMPLETE_EXPLORATORY_CONTEXT_COMPARISON','identities_per_comparison':expected,
        'paired_reference_comparisons':2,'labels_use_strict_zero_threshold_no_tuning':True,
        'stable_counts':{ref:{c:int(np.sum(d.stable_pattern==c)) for c in CATEGORIES} for ref,d in detail.groupby('prostate_reference')},
        'interpretation':'Reference-specific transcriptional opposition, not treatment efficacy or a pure drug-by-cancer interaction',
        'gene_universes':'Original reference-specific masks retained; scores are not calibrated effect sizes across contexts',
        'limitations':['Both halves reuse experiments; agreement is technical rather than biological validation',
        'The PC3 reference is unreplicated and disagrees with patient tissue; comparisons stay separate',
        'MCF7_higher/PC3_higher is relative score preference; either score can be negative',
        'Pattern labels use the same data shown; no separate group-enrichment significance is claimed']}

def plot(detail,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    fig,axes=plt.subplots(2,2,figsize=(11,10),sharex=True,sharey=True)
    values=detail[['MCF7_score_half0','MCF7_score_half1','PC3_score_half0','PC3_score_half1']].to_numpy()
    bound=max(.25,float(np.ceil(np.nanmax(np.abs(values))/.1)*.1))
    for ri,ref in enumerate(PROSTATE):
        q=detail[detail.prostate_reference.eq(ref)]
        for half,ax in enumerate(axes[ri]):
            for c in ['uncertain','neither','both','PC3_only','MCF7_only']:
                d=q[q.stable_pattern.eq(c)]
                ax.scatter(d[f'MCF7_score_half{half}'],d[f'PC3_score_half{half}'],s=9,color=COLORS[c],alpha=.5,edgecolors='none',rasterized=True)
            ax.axhline(0,color='.35',lw=.8);ax.axvline(0,color='.35',lw=.8)
            ax.plot([-bound,bound],[-bound,bound],ls=':',lw=.8,color='.6')
            ax.set_xlim(-bound,bound);ax.set_ylim(-bound,bound);ax.set_aspect('equal',adjustable='box')
            ax.set_title(('PC3/RWPE1 reference' if ri==0 else 'Prostate patient-tissue proxy')+f' | half {half}',fontsize=10)
            ax.set_xlabel('MCF7 reference-opposition score')
            ax.set_ylabel('PC3 reference-opposition score')
    handles=[Line2D([],[],marker='o',ls='',color=COLORS[c],label=NAMES[c]) for c in CATEGORIES]
    fig.legend(handles=handles,loc='lower center',ncol=3,frameon=False,fontsize=9,bbox_to_anchor=(.5,.025))
    fig.suptitle('Each point is one chemical identity; colors require agreement across technical halves',fontsize=12)
    fig.text(.5,.005,'Positive = opposition to the specified reference. No phenotype reversal or therapeutic selectivity established.',ha='center',fontsize=9,color='.4')
    fig.tight_layout(rect=(0,.10,1,.95))
    for ext in ['png','pdf','svg']:fig.savefig(out/('Figure_S3_chemical_context_opposition.'+ext),dpi=200,bbox_inches='tight')
    plt.close(fig)

def main():
    source=ROOT/'inputs/cancer_reversal'
    out=ROOT/'figures'
    out.mkdir(exist_ok=True)
    source_summary=json.loads((source/'SUMMARY.json').read_text())
    if source_summary['status']!='COMPLETE_EXPLORATORY_OPPOSITION_ONLY':
        raise ValueError('Unqualified source score export')
    scores=pd.read_csv(source/'identity_reversal_scores.tsv',sep='\t')
    result=compare(scores,out)
    expected=json.loads((source/'expected_context_summary.json').read_text())
    if result['stable_counts']!=expected['stable_counts']:
        raise ValueError('Frozen chemical grouping does not reproduce')
    result['new_models_fitted']=False
    result['new_statistical_tests']=False
    (out/'Figure_S4_build_record.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('S4 complete: 4,185 identities in each reference comparison; frozen categories reproduced')
if __name__=='__main__':main()

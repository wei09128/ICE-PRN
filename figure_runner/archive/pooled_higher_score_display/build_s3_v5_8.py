#!/usr/bin/env python3
"""Pooled observed-response opposition: distributions and paired score preference."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
P=Path(__file__).resolve().parent
BREAST='MCF7_vs_MCF10A'
PROSTATE=['PC3_vs_RWPE1','prostate_tumor_vs_adjacent_15pairs']
BLUE='#1b6fa8';ORANGE='#c55f2c'

def main():
 source=P/'inputs/cancer_reversal';out=P/'figures';out.mkdir(exist_ok=True)
 s=pd.read_csv(source/'identity_reversal_scores.tsv',sep='\t')
 assert not s.duplicated(['reference','identity_group']).any()
 x=s[s.reference.eq(BREAST)].set_index('identity_group').reversal_spearman_halfmean
 assert len(x)==4185 and np.isfinite(x).all()
 detail=[];summary=[];series=[('MCF7 / MCF10A',x.to_numpy(),BLUE,'-')]
 for ref,label,ls in zip(PROSTATE,['PC3 / RWPE1','PC3 / prostate tissue'],['-','--']):
  q=s[s.reference.eq(ref)].set_index('identity_group');assert set(q.index)==set(x.index)
  y=q.loc[x.index,'reversal_spearman_halfmean'].to_numpy();assert np.isfinite(y).all()
  delta=x.to_numpy()-y
  detail.extend(dict(identity_group=key,prostate_reference=ref,MCF7_pooled_score=float(v),PC3_pooled_score=float(w),MCF7_minus_PC3=float(v-w)) for key,v,w in zip(x.index,x,y))
  summary.append(dict(prostate_reference=ref,n=len(x),median_MCF7=float(np.median(x)),median_PC3=float(np.median(y)),median_paired_difference=float(np.median(delta)),mean_paired_difference=float(np.mean(delta)),MCF7_higher=int(np.sum(delta>0)),PC3_higher=int(np.sum(delta<0)),tied=int(np.sum(delta==0)),MCF7_positive=int(np.sum(x>0)),PC3_positive=int(np.sum(y>0))))
  series.append((label,y,ORANGE,ls))
 pd.DataFrame(detail).to_csv(out/'pooled_paired_scores.tsv',sep='\t',index=False)
 pd.DataFrame(summary).to_csv(out/'pooled_context_summary.tsv',sep='\t',index=False)
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
 fig,(ax,bx)=plt.subplots(1,2,figsize=(12,4.6),gridspec_kw={'width_ratios':[1.15,1]})
 fig.subplots_adjust(left=.08,right=.98,bottom=.18,top=.86,wspace=.48)
 for label,values,color,ls in series:
  v=np.sort(values);ax.step(v,np.arange(1,len(v)+1)/len(v),where='post',color=color,ls=ls,lw=2,label=f'{label} (median {np.median(v):+.3f})')
 ax.axvline(0,color='.5',lw=.8,ls=':');ax.axhline(.5,color='.8',lw=.8,ls=':')
 bound=max(.25,np.ceil(max(np.max(np.abs(v)) for _,v,_,_ in series)*10)/10)
 ax.set_xlim(-bound,bound);ax.set_ylim(0,1);ax.yaxis.set_major_formatter(PercentFormatter(1))
 ax.set_xlabel('Pooled reference-opposition score');ax.set_ylabel('Cumulative share of chemical identities')
 ax.legend(frameon=False,fontsize=8,loc='upper left');ax.set_title('A  Overall score distributions',loc='left',weight='bold',pad=15)
 for i,z in enumerate(summary):
  m=z['MCF7_higher']/z['n'];v=z['PC3_higher']/z['n'];t=z['tied']/z['n']
  bx.barh(i,m,color=BLUE,height=.45,label='MCF7 score higher' if i==0 else None)
  bx.barh(i,v,left=m,color=ORANGE,height=.45,label='PC3 score higher' if i==0 else None)
  if t:bx.barh(i,t,left=m+v,color='.6',height=.45)
  bx.text(m/2,i,f'{m:.1%}',ha='center',va='center',color='white',weight='bold')
  bx.text(m+v/2,i,f'{v:.1%}',ha='center',va='center',color='white',weight='bold')
 bx.axvline(.5,color='.3',ls=':',lw=1);bx.set_xlim(0,1);bx.set_ylim(1.7,-.8)
 bx.set_yticks([0,1],['Cell-line\nreference','Tissue\nreference']);bx.xaxis.set_major_formatter(PercentFormatter(1))
 bx.set_xlabel('Share of paired chemical identities');bx.set_ylabel('Prostate reference')
 bx.set_title('B  Which cell line has the higher score?',loc='left',weight='bold',pad=15)
 bx.legend(frameon=False,fontsize=9,loc='upper center',bbox_to_anchor=(.5,.97),ncol=2)
 for ext in ['png','pdf','svg']:fig.savefig(out/('Figure_S3_chemical_context_opposition.'+ext),dpi=240,bbox_inches='tight')
 plt.close(fig)
 status=dict(status='COMPLETE_EXPLORATORY_POOLED_DESCRIPTIVE',chemical_identities=4185,averaging='Opposite technical-half response vectors averaged equally before Spearman scoring; not a mean of correlations',new_models_fitted=False,new_significance_tests=False,summary=summary,interpretation='Reference-specific score distributions and paired preferences; not cancer-selective efficacy or a pure cancer interaction')
 (out/'Figure_S3_pooled_build_record.json').write_text(json.dumps(status,indent=2)+'\n')
 print(pd.DataFrame(summary).to_string(index=False))
if __name__=='__main__':main()

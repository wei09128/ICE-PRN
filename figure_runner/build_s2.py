from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=Path(__file__).resolve().parent
out=P/'figures';out.mkdir(exist_ok=True)
scores=pd.read_csv(P/'inputs/identity_program/group_scores.tsv',sep='\t')
contrast=pd.read_csv(P/'inputs/identity_program/program_vs_complement_contrasts.tsv',sep='\t')
programs=['WELCH_GATA1_TARGETS','HALLMARK_HEME_METABOLISM','HALLMARK_APICAL_JUNCTION','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION']
labels=['GATA1-responsive (10)','Heme metabolism (99)','Apical junction (74)','EMT (39)']
directions=['K562_to_RPE1','RPE1_to_K562'];colors=['#1b6fa8','#c55f2c']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
fig,ax=plt.subplots(figsize=(7,3.1));axs=[ax]
fig.subplots_adjust(left=.34,right=.98,bottom=.22,top=.86,wspace=.3)
for di,direction in enumerate(directions):
    values=scores[scores.direction.eq(direction)&scores.method.eq('diagonal_plus_cross_gene')].set_index('group').loc[programs]
    x=contrast[contrast.direction.eq(direction)].set_index('program').loc[programs]
    y=np.arange(4)+(di-.5)*.22
    axs[0].errorbar(values.skill*100,y,xerr=np.array([values.skill-values.lower95,values.upper95-values.skill])*100,fmt='o',color=colors[di],capsize=3,label=direction.replace('_to_',' → '))
for ax in axs:
    ax.set_ylim(3.5,-.5)
    ax.set_yticks(np.arange(4))
    ax.axvline(0,lw=.9,color='.5',ls='--')
axs[0].set_yticklabels(labels)
axs[0].set_xlim(0,23)
axs[0].set_xlabel('Learned improvement over destination mean (%)')
axs[0].set_title('Skill within programs',loc='left',fontweight='bold',fontsize=11,pad=15)
axs[0].legend(loc='lower right',frameon=False,fontsize=9)
for ext in ['png','pdf','svg']:fig.savefig(out/f'Figure_S2_identity_programs.{ext}',dpi=220,facecolor='white')
plt.close(fig)

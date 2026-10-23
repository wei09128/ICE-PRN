"""Figure 3 (v5.11): biological associations of learned transfer.
A-B: external p53-hit coverage across outcomes and its sensitivity to matching.
C: identity-related response programs. D: continuous cancer-pathway prediction skill.
Rank panels moved to Figure 4.
Saved estimates only; no refitting."""
from pathlib import Path
import numpy as np, pandas as pd
from matplotlib.ticker import PercentFormatter
from _style import plt, DIR_COLORS, DIR_NAMES, GROUP_COLORS, GROUP_TICKS, MUTED, panel, save

P = Path(__file__).resolve().parent
t = pd.read_csv(P / 'inputs/rank/full_rank_curve.tsv', sep='\t')
cov = pd.read_csv(P / 'inputs/coverage_intervals.tsv', sep='\t')
gap = pd.read_csv(P / 'inputs/recovery_gap_intervals.tsv', sep='\t')
unm = __import__('json').loads((P / 'uncertainty_summary.json').read_text())['unmatched_B2']
dirs = ['K562_to_RPE1', 'RPE1_to_K562']; ranks = [0, 1, 2, 4, 8, 16, 32]

fig = plt.figure(figsize=(16, 12))
outer = fig.add_gridspec(2, 2, left=.08, right=.98, top=.93, bottom=.10, wspace=.40, hspace=.45, height_ratios=[1,1.05])
def heading(slot,letter,title):
    box=slot.get_position(fig)
    fig.text(box.x0-.023,box.y1+.032,letter,weight='bold',fontsize=16,va='bottom')
    fig.text(box.x0,box.y1+.035,title,weight='bold',fontsize=11,va='bottom')
heading(outer[0,0],'A','External p53 hits across prediction outcomes')
heading(outer[0,1],'B','p53 rescue gap and matching sensitivity')

ax = fig.add_subplot(outer[0,0])
groups = ['rescued', 'persistent_failure', 'both_succeed', 'other_mixed_or_tied']
counts = [(314, 1751), (71, 625), (19, 247), (98, 1125)]
q = cov[(cov.scope == 'pooled') & (cov.denominator == 'all_cases')].set_index('group').loc[groups]
for j, (g, (n, den)) in enumerate(zip(groups, counts)):
    r = q.loc[g]; v = n / den; assert np.isclose(v, r.coverage)
    ax.bar(j, v, color=GROUP_COLORS[j], width=.62)
    ax.errorbar(j, v, yerr=[[v - r.lower95], [r.upper95 - v]], fmt='none', color='.25', capsize=3, lw=.9)
    ax.text(j, r.upper95 + .008, f'{v:.1%}\n{n}/{den:,}', ha='center', fontsize=7.8)
ax.set_xticks(range(4), ['Rescued', 'Persistent\nfailure', 'Both\nsucceed', 'Other'], fontsize=9); ax.set_ylim(0, .27)
ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0)); ax.set_ylabel('External p53 hits / target–direction cases')


ax = fig.add_subplot(outer[0,1])
rows = [('Unmatched', 'pooled', unm['estimate'] * 100, unm['lower'] * 100, unm['upper'] * 100)]
for var, lab in [('original_3cov_full_panel', 'Magnitude + cells'),
                 ('original_3cov_complete_baselines', 'Baseline subset'),
                 ('expanded_5cov_complete_baselines', '+ baseline expr.')]:
    r = gap[(gap.variant == var) & (gap.direction == 'pooled')].iloc[0]
    rows.append((lab, 'pooled', r.gap_pp, r.lower95_pp, r.upper95_pp))
for d in dirs:
    r = gap[(gap.variant == 'expanded_5cov_complete_baselines') & (gap.direction == d)].iloc[0]
    rows.append((f'+ baseline expr.\n{DIR_NAMES[d]}', d, r.gap_pp, r.lower95_pp, r.upper95_pp))
for j, (lab, d, v, lo, hi) in enumerate(rows):
    c = '#333333' if d == 'pooled' else DIR_COLORS[d]
    ax.errorbar(v, j, xerr=[[v - lo], [hi - v]], fmt='o' if d == 'pooled' else 'D', color=c, ms=6, capsize=3,
                mfc=c if d == 'pooled' else 'white')
    ax.annotate(f'{v:.1f}', (v, j), xytext=(0, 7), textcoords='offset points', ha='center', fontsize=7.8)
ax.axhline(3.5, color='.85', lw=.8); ax.axvline(0, color='.4', ls='--', lw=.8)
ax.set_yticks(range(len(rows)), [r[0] for r in rows], fontsize=8.2); ax.invert_yaxis(); ax.set_xlim(-10, 25)
ax.set_xlabel('p53-hit − below-threshold rescue rate\n(percentage points)', fontsize=8)



identity=pd.read_csv(P/'inputs/identity_program/program_vs_complement_contrasts.tsv',sep='\t')
drug=pd.read_csv(P/'inputs/drug_pathway_skill/primary_endpoints.tsv',sep='\t')
assert __import__('json').loads((P/'inputs/drug_pathway_skill/SUMMARY.json').read_text())['status']=='COMPLETE_EXPLORATORY_REVIEW_REQUIRED'
assert len(drug)==4 and drug[['estimate','lower98_75','upper98_75']].notna().all().all()
COLORS=[DIR_COLORS[d] for d in dirs];DIRS=dirs;NAMES=[DIR_NAMES[d] for d in dirs]

def forest(ax,values,lower,upper,ys,color,marker='o'):
 for v,lo,hi,y in zip(values,lower,upper,ys):
  # Draw endpoints directly: bootstrap intervals need not bracket the point estimate.
  ax.hlines(y,lo,hi,color=color,lw=1.7)
  ax.plot([lo,hi],[y,y],linestyle='none',marker='|',color=color,ms=6)
  ax.plot(v,y,marker,color=color,ms=6,zorder=3)


# C: identity readouts. Both directions and all pre-existing contrasts are retained.
heading(outer[1,0],'C','Identity-related response programs')
ax=fig.add_subplot(outer[1,0])
plist=['WELCH_GATA1_TARGETS','HALLMARK_HEME_METABOLISM','HALLMARK_APICAL_JUNCTION','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION']
plabels=['GATA1-responsive (10)','Heme metabolism (99)','Apical junction (74)','EMT (39)']
for di,d in enumerate(DIRS):
 q=identity[identity.direction==d].set_index('program').loc[plist]
 yy=np.arange(4)+(di-.5)*.24
 forest(ax,q.estimate*100,q.lower99_375*100,q.upper99_375*100,yy,COLORS[di],'o' if di==0 else 's')
 ax.plot([],[],'o' if di==0 else 's',color=COLORS[di],label=NAMES[di])
 for pi,v in enumerate(q.estimate):
  ax.annotate(f'{v*100:+.1f}',(v*100,yy[pi]),xytext=(0,-12 if di else 9),textcoords='offset points',fontsize=7,ha='center',color=COLORS[di])
ax.set_yticks(range(4),plabels,fontsize=9);ax.set_ylim(3.75,-.7)
ax.axvline(0,color='.5',ls='--',lw=.8);ax.set_xlim(-4,11)
ax.set_xlabel('Learned skill: program − other measured genes (points)')
ax.legend(frameon=False,fontsize=8,loc='lower right')

# D: broad cancer and estrogen contrasts share the same relative-skill question.
heading(outer[1,1],'D','Chemical response programs and matched controls')
e=pd.read_csv(P/'inputs/readout_followups/estrogen_primary_endpoints.tsv',sep='\t')
dm=pd.read_csv(P/'inputs/readout_followups/drug_matched_endpoints.tsv',sep='\t')
ax=fig.add_subplot(outer[1,1]);labels=[]
items=[('MCF7_to_PC3','prostate_cancer','Prostate cancer (23)\nMCF7 → PC3'),('PC3_to_MCF7','breast_cancer','Breast cancer (33)\nPC3 → MCF7'),('PC3_to_MCF7','HALLMARK_ESTROGEN_RESPONSE_EARLY','Estrogen early (30)\nPC3 → MCF7'),('PC3_to_MCF7','HALLMARK_ESTROGEN_RESPONSE_LATE','Estrogen late (34)\nPC3 → MCF7')]
for j,(d,gs,lab) in enumerate(items):
 labels.append(lab);c=DIR_COLORS[d]
 if j<2:
  u=drug[(drug.direction==d)&drug.endpoint.eq('pathway_minus_complement_skill')].iloc[0];m=dm[dm.direction==d].iloc[0]
  rows=[(u.estimate,u.lower98_75,u.upper98_75,'complement','UNMATCHED'),(m.matched_skill_difference,m.lower97_5,m.upper97_5,'matched',m.matching_status)]
 else:
  q=e[e.gene_set==gs];rows=[(r.estimate,r.lower98_75,r.upper98_75,r.reference,r.matching_status) for r in q.itertuples()]
 for v,lo,hi,ref,status in rows:
  yy=j+(-.15 if ref=='complement' else .15);mk='o' if ref=='complement' else 'D'
  ax.hlines(yy,100*lo,100*hi,color=c,lw=1.6)
  ax.plot(100*v,yy,marker=mk,color=c,mfc='white' if status=='INADEQUATE_MATCHING_SENSITIVITY_ONLY' else c,ms=6)
ax.set_yticks(range(4),labels,fontsize=8);ax.set_ylim(3.7,-.65);ax.set_xlim(-2,4.3);ax.axvline(0,color='.5',ls='--',lw=.8)
ax.set_xlabel('Learned skill: program − reference readouts (points)',fontsize=8)
ax.plot([],[],'o',color='.3',label='Other readouts');ax.plot([],[],'D',color='.3',label='Matched readouts');ax.plot([],[],'D',color='.3',mfc='white',label='Inadequate balance')
ax.legend(frameon=False,fontsize=7.5,loc='upper center',bbox_to_anchor=(.5,-.17),ncol=3)

(P/'figures').mkdir(exist_ok=True)
for ext in ['png','pdf','svg']:
    fig.savefig(P/'figures'/f'Figure_3_biological_associations.{ext}',dpi=240,facecolor='white',bbox_inches='tight',pad_inches=.2)
plt.close(fig)
print('Figure 3 complete: p53, identity and drug panels')

"""Pooled and direction-specific p53 coverage from frozen counts and intervals."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
P=Path(__file__).resolve().parent
O=P/'figures';O.mkdir(exist_ok=True)
ci=pd.read_csv(P/'coverage_intervals.tsv',sep='\t')
groups=['rescued','persistent_failure','both_succeed','other_mixed_or_tied']
scopes=['pooled','K562_to_RPE1','RPE1_to_K562']
titles=['A   Pooled','B   K562 → RPE1','C   RPE1 → K562']
counts=[[(314,1751),(71,625),(19,247),(98,1125)],[(153,762),(33,293),(18,229),(47,590)],[(161,989),(38,332),(1,18),(51,535)]]
colors=['#008f8a','#b97a22','#5b72b5','#9a9a9a']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
fig,axs=plt.subplots(1,3,figsize=(11.5,3.9),sharey=True)
fig.subplots_adjust(left=.07,right=.98,top=.88,bottom=.22,wspace=.22)
rows=[]
for ax,scope,title,nn in zip(axs,scopes,titles,counts):
 q=ci[ci.scope.eq(scope)&ci.denominator.eq('all_cases')].set_index('group')
 for j,(group,(n,den),color) in enumerate(zip(groups,nn,colors)):
  r=q.loc[group];v=n/den;assert np.isclose(v,r.coverage)
  ax.bar(j,v,color=color,width=.65)
  ax.errorbar(j,v,yerr=[[v-r.lower95],[r.upper95-v]],fmt='none',color='.25',capsize=3)
  ax.text(j,r.upper95+.009,f'{v:.1%}\n{n}/{den:,}',ha='center',fontsize=8)
  rows.append({'scope':scope,'group':group,'hits':n,'cases':den,'coverage':v,'lower95':r.lower95,'upper95':r.upper95})
 ax.set_xticks(range(4),['Rescued','Persistent\nfailure','Both\nsucceeded','Other'],fontsize=8)
 ax.set_ylim(0,.35);ax.yaxis.set_major_formatter(PercentFormatter(1));ax.set_title(title,loc='left',fontweight='bold')
axs[0].set_ylabel('External p53-hit annotation coverage')
for ext in ['pdf','svg','png']:fig.savefig(O/f'Figure_S1_p53_direction_details.{ext}',dpi=240,facecolor='white')
plt.close(fig)
pd.DataFrame(rows).to_csv(O/'Figure_S1_values.tsv',sep='\t',index=False)
print('S1 complete: p53 coverage by direction, replacing historical context assembly.')

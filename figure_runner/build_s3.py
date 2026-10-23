"""Figure 4 (v5.11): capacity of the learned correction and where added capacity goes.
A-B: held-out skill across preset correction ranks (formerly Figure 3A).
C: skill within fixed Fischer (2017) p53 target-gene readouts versus all other readouts, by rank and direction.
D: p53-minus-other gain over affine by rank (12 fixed contrasts; adjusted 99.583% intervals).
Saved estimates only; no refitting. Exploratory readout endpoint (contract fixed before rank-specific p53 losses)."""
from pathlib import Path
import json, numpy as np, pandas as pd
from _style import plt, DIR_COLORS, DIR_NAMES, MUTED

P = Path(__file__).resolve().parent
t = pd.read_csv(P / 'inputs/rank/full_rank_curve.tsv', sep='\t')
summ = pd.read_csv(P / 'inputs/rank/rank_capacity_summary.tsv', sep='\t')
cur = pd.read_csv(P / 'inputs/rank_p53/rank_readout_curves.tsv', sep='\t')
con = pd.read_csv(P / 'inputs/rank_p53/primary_contrasts.tsv', sep='\t')
meta = json.loads((P / 'inputs/rank_p53/SUMMARY.json').read_text())
assert meta['primary_family_size'] == 12 and len(con) == 12
dirs = ['K562_to_RPE1', 'RPE1_to_K562']; ranks = [0, 1, 2, 4, 8, 16, 32]
OTHER = '#8a8a8a'

sc=pd.read_csv(P/'inputs/ribosome_p53/secondary_rank_contrasts.tsv',sep='\t')
fig=plt.figure(figsize=(14,5))
outer=fig.add_gridspec(1,2,left=.07,right=.98,top=.83,bottom=.16,wspace=.4,width_ratios=[1.4,1])
x=np.arange(7)
def heading(slot,letter,title):
 box=slot.get_position(fig);fig.text(box.x0-.025,box.y1+.06,letter,weight='bold',fontsize=15);fig.text(box.x0,box.y1+.06,title,weight='bold',fontsize=10)
# C: p53 readouts vs other readouts
n53 = meta['readout_counts']['p53']; nother = meta['readout_counts']['complement']
heading(outer[0], 'A', 'Skill within p53 target readouts versus other readouts')
cg = outer[0].subgridspec(1, 2, wspace=.12)
for k, d in enumerate(dirs):
    ax = fig.add_subplot(cg[k])
    for grp, col, mk, lab in [('p53', DIR_COLORS[d], 'o', f'p53 targets ({n53})'), ('complement', OTHER, 's', f'Other readouts ({nother:,})')]:
        q = cur[(cur.direction == d) & (cur.group == grp) & cur.method.str.match(r'rank_\d+$')].set_index('rank').loc[ranks]
        ax.plot(x, q.skill * 100, mk + '-', color=col, ms=5.5, lw=1.8, label=lab, mfc=col if grp == 'p53' else 'white')
        ax.fill_between(x, q.lower95 * 100, q.upper95 * 100, color=col, alpha=.12)
        ax.annotate(f'{q.skill.iloc[-1]*100:.1f}%', (6, q.skill.iloc[-1] * 100), xytext=(6, -3), textcoords='offset points', fontsize=8, color=col)
    ax.set_xticks(x, ['0\naffine', '1', '2', '4', '8', '16', '32']); ax.set_xlim(-.3, 7.3); ax.set_ylim(0, 21.5)
    ax.set_title(DIR_NAMES[d], fontsize=10, color=DIR_COLORS[d], pad=8)
    ax.set_xlabel('Correction rank')
    if k == 0: ax.set_ylabel('Skill against destination mean (%)'); ax.legend(frameon=False, fontsize=8.2, loc='upper left')
    else: ax.set_yticklabels([])

ax = fig.add_subplot(outer[1]); rk = [1, 2, 4, 8, 16, 32]; x = np.arange(6)
for j, d in enumerate(dirs):
    q = sc[sc.direction == d].set_index('rank').loc[rk]; xs = x + (j - .5) * .2
    ax.vlines(xs, q.lower99_583333 * 100, q.upper99_583333 * 100, color=DIR_COLORS[d], lw=1.8)
    ax.plot(xs, q.estimate * 100, 'o' if j == 0 else 's', color=DIR_COLORS[d], ms=6, label=DIR_NAMES[d])
ax.axhline(0, color='.5', ls='--', lw=.8); ax.set_xticks(x, [str(r) for r in rk]); ax.set_xlabel('Common-basis correction rank')
ax.set_ylabel('Preferential p53 gain difference (points)'); ax.legend(frameon=False, fontsize=8, loc='upper left')
ax.set_title('B   Ribosomal target contrast across ranks', loc='left', fontweight='bold', fontsize=10, pad=10)

(P / 'figures').mkdir(exist_ok=True)
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(P / 'figures' / f'Figure_S3_capacity_details.{ext}', dpi=240, facecolor='white', bbox_inches='tight', pad_inches=.2)
plt.close(fig)
print('Figure 4 complete: rank capacity and p53 readout panels')

"""Figure S2 (v5): identity-TF and drug-mechanism annotations (formerly v4.16 Fig. 2B-C)."""
from pathlib import Path
import numpy as np, pandas as pd
from matplotlib.ticker import PercentFormatter
from _style import plt, GROUP_COLORS, GROUP_TICKS, MUTED, panel, save

P = Path(__file__).resolve().parent
groups = ['rescued', 'persistent_failure', 'both_succeeded', 'other_mixed_or_tied']
ic = pd.read_csv(P / 'inputs/identity_targets/broad_atlas_outcome_coverage.tsv', sep='\t')
ic = ic[ic.scope == 'pooled'].set_index('group').loc[groups]
q = pd.read_csv(P / 'inputs/drug/moa_all_64_contrasts.tsv', sep='\t'); classes = sorted(q.moa.unique())
def short(m):
    m = m.replace('bacterial cell wall synthesis inhibitor', 'bacterial cell-wall inhibitor').replace('receptor', 'R').replace('channel', '-channel').replace(' -', '-')
    return (m[0].upper() + m[1:]).replace('dna', 'DNA')

fig = plt.figure(figsize=(11.5, 6.6))
a = fig.add_axes([.07, .25, .22, .6]); c = fig.add_axes([.52, .14, .21, .71]); d = fig.add_axes([.765, .14, .21, .71], sharey=c)
for j, g in enumerate(groups):
    r = ic.loc[g]; n, den = int(r.core_TF_cases), int(r.all_scored_TF_cases); v = n / den
    a.bar(j, v, color=GROUP_COLORS[j], width=.62)
    a.errorbar(j, v, yerr=[[v - r.lower95], [r.upper95 - v]], fmt='none', color='.25', capsize=3, lw=.9)
    a.text(j, r.upper95 + .02, f'{v:.0%}\n{n}/{den}', ha='center', fontsize=7.8)
a.set_xticks(range(4), GROUP_TICKS, fontsize=8); a.set_ylim(0, 1.0); a.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
a.set_ylabel('Core-TF candidates / atlas-scored TF cases'); panel(a, 'A', 'Broad identity-TF annotation')
for ax, direction, title, n in [(c, 'MCF7_to_PC3', 'MCF7 → PC3', 14), (d, 'PC3_to_MCF7', 'PC3 → MCF7', 12)]:
    for i, contrast in enumerate(['rescued_minus_both_success', 'persistent_failure_minus_both_success']):
        r = q[(q.direction == direction) & (q.contrast == contrast)].set_index('moa').loc[classes]
        v = r.coverage_difference.to_numpy() * 100; lo = r.conservative_family_lower.to_numpy() * 100; hi = r.conservative_family_upper.to_numpy() * 100
        ax.errorbar(v, np.arange(16) + (i - .5) * .25, xerr=[v - lo, hi - v], fmt='os'[i], ms=3.5, lw=.8, capsize=1.5,
                    color=GROUP_COLORS[i], label=['Rescued − both succeeded', 'Persistent − both succeeded'][i])
    ax.axvline(0, color='.4', lw=.7); ax.set_xlim(-55, 16); ax.set_xticks([-50, -25, 0])
    ax.set_title(f'{title}\n(both-success annotated n = {n})', fontsize=8.5); ax.set_xlabel('Coverage difference (points)', fontsize=8)
c.set_yticks(np.arange(16), [short(m) for m in classes], fontsize=7.5); c.invert_yaxis(); d.tick_params(axis='y', labelleft=False)
fig.text(.52, .93, 'B   Drug-mechanism coverage', fontweight='bold', fontsize=10)
h, l = c.get_legend_handles_labels(); fig.legend(h, l, loc='lower center', bbox_to_anchor=(.74, .0), ncol=2, frameon=False, fontsize=8)
save(fig, 'Figure_S2_identity_and_drug_annotations')
print('S2 complete')

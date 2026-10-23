"""Figure 1 (v5): copying fails and learned correction recovers part of the response,
in genetic (CRISPRi) and chemical (LINCS) transfer. Saved estimates only; no refitting."""
from pathlib import Path
import numpy as np, pandas as pd
from _style import plt, DIR_COLORS, DIR_NAMES, MUTED, panel, save

P = Path(__file__).resolve().parent
rank = pd.read_csv(P / 'inputs/rank/full_rank_curve.tsv', sep='\t')
bench = pd.read_csv(P / 'benchmark_CI.tsv', sep='\t')
lincs = pd.read_csv(P / 'inputs/drug/transfer_skill.tsv', sep='\t')
norm = pd.read_csv(P / 'inputs/normalized_cross_gene_skill.tsv', sep='\t')

rows = []  # system, direction, method, skill, lo, hi
for d in ['K562_to_RPE1', 'RPE1_to_K562']:
    z = rank[(rank.direction == d) & (rank.method == 'zero')].iloc[0]
    rows.append(('CRISPRi', d, 'zero', z.skill, z.lower95, z.upper95))
    for m, key in [('source_copy', 'copy'), ('diagonal_affine', 'affine'), ('diagonal_plus_cross_gene', 'learned')]:
        r = bench[(bench.direction == d) & (bench.method == m)].iloc[0]
        rows.append(('CRISPRi', d, key, r.skill, r.lower, r.upper))
    n = norm[norm.direction == d].iloc[0]
    rows.append(('CRISPRi, library-size normalized', d, 'learned', n.skill, n.lower95, n.upper95))
for d in ['MCF7_to_PC3', 'PC3_to_MCF7']:
    for m, key in [('zero', 'zero'), ('source_copy', 'copy'), ('diagonal_affine', 'affine'), ('diagonal_plus_rank16', 'learned')]:
        r = lincs[(lincs.direction == d) & (lincs.method == m)].iloc[0]
        rows.append(('LINCS chemical', d, key, r.skill_vs_destination_mean, r.lower95, r.upper95))
T = pd.DataFrame(rows, columns=['system', 'direction', 'predictor', 'skill', 'lower95', 'upper95'])

fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.3), gridspec_kw={'width_ratios': [1, 1.15]})
fig.subplots_adjust(left=.115, right=.985, top=.9, bottom=.17, wspace=.42)

# A: naive predictors. Both systems, both directions; copy is worse than predicting no change everywhere.
ax = axs[0]
groups = [('CRISPRi', 'K562_to_RPE1'), ('CRISPRi', 'RPE1_to_K562'), ('LINCS chemical', 'MCF7_to_PC3'), ('LINCS chemical', 'PC3_to_MCF7')]
for gi, (s, d) in enumerate(groups):
    for k, (pred, mk) in enumerate([('zero', 's'), ('copy', 'o')]):
        r = T[(T.system == s) & (T.direction == d) & (T.predictor == pred)].iloc[0]
        y = gi + (k - .5) * .32
        ax.errorbar(r.skill * 100, y, xerr=[[100 * (r.skill - r.lower95)], [100 * (r.upper95 - r.skill)]],
                    fmt=mk, ms=7, color=DIR_COLORS[d], mfc='white' if pred == 'zero' else DIR_COLORS[d],
                    capsize=2.5, lw=1.2)
        ax.annotate(f'{r.skill*100:.1f}%', (r.skill * 100, y), xytext=(0, 7), textcoords='offset points',
                    ha='center', fontsize=7.5, color='#333333')
ax.set_yticks(range(4), [f'{DIR_NAMES[d]}\n{"genetic" if s=="CRISPRi" else "chemical"}' for s, d in groups], fontsize=8.5)
ax.invert_yaxis(); ax.axvline(0, color='.55', lw=.8); ax.set_xlim(-215, 15)
ax.axhline(1.5, color='.85', lw=.8)
ax.set_xlabel('Skill against destination mean (%)')
panel(ax, 'A', 'Copying and zero prediction')
ax.plot([], [], 's', mfc='white', mec='.3', ms=7, label='Predict no change (zero)')
ax.plot([], [], 'o', color='.3', ms=7, label='Copy source response')
ax.legend(frameon=False, fontsize=8, loc='lower left')

# B: learned predictors.
ax = axs[1]
cats = [('CRISPRi', 'affine'), ('CRISPRi', 'learned'), ('CRISPRi, library-size normalized', 'learned'),
        ('LINCS chemical', 'affine'), ('LINCS chemical', 'learned')]
labels = ['Genetic: per-gene affine', 'Genetic: + cross-gene', 'Genetic: + cross-gene,\nnormalized counts',
          'Chemical: per-gene affine', 'Chemical: + cross-gene']
for ci, (s, pred) in enumerate(cats):
    sub = T[(T.system == s) & (T.predictor == pred)]
    for j, r in enumerate(sub.itertuples()):
        y = ci + (j - .5) * .3
        ax.errorbar(r.skill * 100, y, xerr=[[100 * (r.skill - r.lower95)], [100 * (r.upper95 - r.skill)]],
                    fmt='o', ms=7, color=DIR_COLORS[r.direction], capsize=2.5, lw=1.2)
        ax.annotate(f'{r.skill*100:+.1f}%', (r.upper95 * 100, y), xytext=(5, -3), textcoords='offset points',
                    fontsize=7.5, color='#333333')
ax.set_yticks(range(5), labels, fontsize=8.5); ax.invert_yaxis()
ax.axhline(2.5, color='.85', lw=.8); ax.axvline(0, color='.55', lw=.8); ax.set_xlim(-0.5, 19)
ax.set_xlabel('Skill against destination mean (%)')
panel(ax, 'B', 'Affine and cross-gene prediction')
ax.plot([], [], 'o', color=DIR_COLORS['K562_to_RPE1'], label='K562 → RPE1  /  MCF7 → PC3')
ax.plot([], [], 'o', color=DIR_COLORS['RPE1_to_K562'], label='RPE1 → K562  /  PC3 → MCF7')
ax.legend(frameon=False, fontsize=8, loc='lower right')

save(fig, 'Figure_1_two_system_transfer')
sd = P / 'source_data'; sd.mkdir(exist_ok=True)
T.to_csv(sd / 'Figure_1_values.tsv', sep='\t', index=False)
print('Figure 1 complete')

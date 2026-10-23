"""Figure 2 (v5): how comparator loss and response magnitude shape transfer scores.
Panels re-use the saved v4.16 Figure 1 B-F estimates; no refitting."""
from pathlib import Path
import numpy as np, pandas as pd
from _style import plt, DIR_COLORS, DIR_NAMES, MUTED, panel, save

P = Path(__file__).resolve().parent
b = pd.read_csv(P / 'benchmark_CI.tsv', sep='\t'); risk = pd.read_csv(P / 'risk_CI.tsv', sep='\t')
dirs = ['K562_to_RPE1', 'RPE1_to_K562']; rmap = {'K562->RPE1': 'K562_to_RPE1', 'RPE1->K562': 'RPE1_to_K562'}
PROG = {'HALLMARK_E2F_TARGETS': ('E2F', 'o', '-'), 'HALLMARK_G2M_CHECKPOINT': ('G2M', 's', '--')}

fig = plt.figure(figsize=(11.8, 7.6))
gs = fig.add_gridspec(2, 2, left=.095, right=.98, top=.92, bottom=.12, wspace=.48, hspace=.68)

ax = fig.add_subplot(gs[0, 0])
for i, d in enumerate(dirs):
    r = b[(b.direction == d) & (b.method == 'source_copy')].iloc[0]
    ax.bar(i - .17, r.baseline_loss, width=.32, color=DIR_COLORS[d])
    ax.bar(i + .17, r.loss, width=.32, facecolor='white', edgecolor=DIR_COLORS[d], hatch='//', lw=1.2)
    ax.text(i - .17, r.baseline_loss + .0006, f'{r.baseline_loss:.4f}', ha='center', fontsize=7.5)
    ax.text(i + .17, r.loss + .0006, f'{r.loss:.4f}', ha='center', fontsize=7.5)
sk = [b[(b.direction == d) & (b.method == 'source_copy')].iloc[0].skill * 100 for d in dirs]
ax.set_xticks([0, 1], [f'{DIR_NAMES[d]}\ncopy skill {v:.1f}%' for d, v in zip(dirs, sk)])
ax.set_ylim(0, .037); ax.set_ylabel('Whole-response MSE')
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color='.45', label='Destination mean'), Patch(facecolor='white', edgecolor='.45', hatch='//', label='Copy')], frameon=False, fontsize=8, loc='upper center', ncol=2)
panel(ax, 'A', 'Copy error and comparator loss')

ax = fig.add_subplot(gs[0, 1])
x = np.linspace(-.6, 2.6, 400); y = 1 - (x - 1) ** 2
ax.plot(x, y, color='#333333', lw=2); ax.fill_between(x, y, 0, where=y > 0, color='#1b6fa8', alpha=.15)
ax.axhline(0, color='.4', lw=.8)
for v in (0, 2): ax.axvline(v, color='.6', ls=':', lw=1)
ax.text(1, .25, 'copy beats\n“no change”', ha='center', fontsize=8.5)
ax.text(.12, -1.05, '← wrong sign', ha='left', fontsize=8.5); ax.text(1.88, -1.05, 'overshoot loses →', ha='right', fontsize=8.5)
ax.set_xlabel('Source / destination response  s / y'); ax.set_ylabel('Scalar copy skill against zero')
panel(ax, 'B', 'Scalar response and copy skill')

container = fig.add_subplot(gs[1, 0]); container.axis('off')
panel(container, 'C', 'Risk-ranked coverage')
sub = gs[1, 0].subgridspec(2, 2, height_ratios=[.14, 1], hspace=0, wspace=.65)
for i, d in enumerate(dirs):
    ax = fig.add_subplot(sub[1, i])
    q = risk[risk.direction.map(rmap) == d]
    for prog, (lab, mk, ls) in PROG.items():
        s = q[q.program == prog]
        ax.plot(s.coverage * 100, s.skill * 100, marker=mk, ls=ls, color=DIR_COLORS[d], label=lab, ms=6)
        ax.fill_between(s.coverage * 100, s.lower * 100, s.upper * 100, color=DIR_COLORS[d], alpha=.10)
    ax.axhline(0, color='.4', lw=.8, ls='--'); ax.set_xticks([20, 40, 60, 80, 100])
    ax.set_xlabel('Targets retained (%)', fontsize=8); ax.set_ylabel('Copy skill vs zero (%)', fontsize=8)
    ax.tick_params(labelsize=8)
    ax.set_xticks([20, 60, 100])
    ax.legend(frameon=False, fontsize=8, loc='upper left' if i == 0 else 'lower left')
    ax.set_title(DIR_NAMES[d], fontsize=8.5, pad=6)

ax = fig.add_subplot(gs[1, 1])
q = risk[risk.coverage == .2].copy(); q['d'] = q.direction.map(rmap)
order = [(d, p) for d in dirs for p in PROG]
for j, (d, p) in enumerate(order):
    r = q[(q.d == d) & (q.program == p)].iloc[0]
    ax.errorbar(r.skill * 100, j, xerr=[[100 * (r.skill - r.lower)], [100 * (r.upper - r.skill)]], fmt=PROG[p][1],
                color=DIR_COLORS[d], capsize=3, ms=6)
    ax.annotate(f'{r.skill*100:+.2f}%', (r.skill * 100, j), xytext=(0, 7), textcoords='offset points', ha='center', fontsize=7.5)
ax.set_yticks(range(4), [f'{DIR_NAMES[d]}\n{PROG[p][0]}' for d, p in order], fontsize=8.5); ax.invert_yaxis()
ax.axvline(0, color='.4', ls='--', lw=.8); ax.set_xlim(-30, 5); ax.set_xlabel('Copy skill against zero (%)')
panel(ax, 'D', 'Lowest-risk 20%')

save(fig, 'Figure_2_scoring')
print('Figure 2 complete')

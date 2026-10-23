"""Figure S5 (v5.11): ribosomal-protein knockdowns and p53 target readouts.
A: p53-readout and other-readout skill with affine and the original rank-16 cross-gene model, for KEGG ribosome targets
   (n=72) and other annotated targets (n=1,796), by direction (descriptive).
B: primary contrast, ribosome minus other targets in preferential p53 gain over affine, original rank-16 model (97.5%).
C: secondary common-basis rank grid (99.583%). Saved estimates; no refitting."""
from pathlib import Path
import json, numpy as np, pandas as pd
from _style import plt, DIR_COLORS, DIR_NAMES

P = Path(__file__).resolve().parent
g = pd.read_csv(P / 'inputs/ribosome_p53/all_group_skills.tsv', sep='\t')
pc = pd.read_csv(P / 'inputs/ribosome_p53/primary_contrasts.tsv', sep='\t')
sc = pd.read_csv(P / 'inputs/ribosome_p53/secondary_rank_contrasts.tsv', sep='\t')
meta = json.loads((P / 'inputs/ribosome_p53/SUMMARY.json').read_text())
assert all(meta['parity'].values()) and len(pc) == 2 and len(sc) == 12
dirs = ['K562_to_RPE1', 'RPE1_to_K562']; OTHER = '#8a8a8a'

fig = plt.figure(figsize=(15, 4.8))
gs = fig.add_gridspec(1, 3, width_ratios=[1.5, .7, 1.1], wspace=.42, left=.06, right=.98, top=.85, bottom=.17)

ax = fig.add_subplot(gs[0]); labels = []; y = 0
for d in dirs:
    for grp, glab in [('ribosome', f"Ribosome targets ({meta['group_sizes']['ribosome']})"), ('nonmember', f"Other targets ({meta['group_sizes']['nonmember']:,})")]:
        for ro, col, off in [('p53', DIR_COLORS[d], -.13), ('complement', OTHER, .13)]:
            q = g[(g.target_group == grp) & (g.direction == d) & (g.readouts == ro)].set_index('method').skill * 100
            a, b = q['rank_0'], q['original_rank16']
            ax.annotate('', xy=(b, y + off), xytext=(a, y + off), arrowprops=dict(arrowstyle='->', color=col, lw=1.6))
            ax.plot(a, y + off, 'o', mfc='white', mec=col, ms=5.5); ax.plot(b, y + off, 'o', color=col, ms=6)
            if grp == 'ribosome' and ro == 'p53':
                ax.annotate(f'{a:.1f}% → {b:.1f}%', (b, y + off), xytext=(6, 0), textcoords='offset points', va='center', fontsize=8, color=col)
        labels.append(f'{glab}\n{DIR_NAMES[d]}'); y += 1
ax.set_yticks(range(4), labels, fontsize=8.5); ax.set_ylim(3.6, -.6); ax.axvline(0, color='.5', ls='--', lw=.8); ax.set_xlim(-2, 60)
ax.set_xlabel('Skill against destination mean (%)')
ax.plot([], [], 'o', color=DIR_COLORS[dirs[0]], label='p53 target readouts'); ax.plot([], [], 'o', color=OTHER, label='Other readouts')
ax.legend(frameon=False, fontsize=8, loc='lower right')
ax.set_title('A   Readout skill by target group', loc='left', fontweight='bold', fontsize=10, pad=10)

ax = fig.add_subplot(gs[1])
for j, d in enumerate(dirs):
    r = pc[pc.direction == d].iloc[0]
    ax.vlines(j, r.lower97_5 * 100, r.upper97_5 * 100, color=DIR_COLORS[d], lw=2)
    ax.plot(j, r.estimate * 100, 'o' if j == 0 else 's', color=DIR_COLORS[d], ms=7)
    ax.annotate(f'{r.estimate*100:+.2f}', (j, r.estimate * 100), xytext=(9, 0), textcoords='offset points', va='center', fontsize=8.5, color=DIR_COLORS[d])
ax.axhline(0, color='.5', ls='--', lw=.8); ax.set_xticks([0, 1], [DIR_NAMES[d] for d in dirs], fontsize=8.5); ax.set_xlim(-.5, 1.7)
ax.set_ylabel('Ribosome − other targets:\npreferential p53 gain (points)')
ax.set_title('B   Original rank-16 model', loc='left', fontweight='bold', fontsize=10, pad=10)

ax = fig.add_subplot(gs[2]); rk = [1, 2, 4, 8, 16, 32]; x = np.arange(6)
for j, d in enumerate(dirs):
    q = sc[sc.direction == d].set_index('rank').loc[rk]; xs = x + (j - .5) * .2
    ax.vlines(xs, q.lower99_583333 * 100, q.upper99_583333 * 100, color=DIR_COLORS[d], lw=1.8)
    ax.plot(xs, q.estimate * 100, 'o' if j == 0 else 's', color=DIR_COLORS[d], ms=6, label=DIR_NAMES[d])
ax.axhline(0, color='.5', ls='--', lw=.8); ax.set_xticks(x, [str(r) for r in rk]); ax.set_xlabel('Common-basis correction rank')
ax.set_ylabel('Preferential p53 gain difference (points)'); ax.legend(frameon=False, fontsize=8, loc='upper left')
ax.set_title('C   Rank grid', loc='left', fontweight='bold', fontsize=10, pad=10)

(P / 'figures').mkdir(exist_ok=True)
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(P / 'figures' / f'Figure_S5_ribosome_p53_readouts.{ext}', dpi=240, facecolor='white', bbox_inches='tight', pad_inches=.2)
plt.close(fig); print('Figure S5 complete')

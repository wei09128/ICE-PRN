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

fig = plt.figure(figsize=(16, 11.5))
outer = fig.add_gridspec(2, 2, left=.07, right=.98, top=.93, bottom=.09, wspace=.28, hspace=.45)
def heading(slot, letter, title):
    box = slot.get_position(fig)
    fig.text(box.x0 - .023, box.y1 + .032, letter, weight='bold', fontsize=16, va='bottom')
    fig.text(box.x0, box.y1 + .035, title, weight='bold', fontsize=11, va='bottom')

# A: skill by rank, all readouts
heading(outer[0, 0], 'A', 'Skill by correction rank')
ax = fig.add_subplot(outer[0, 0]); x = np.arange(7)
for d in dirs:
    q = t[(t.direction == d) & t.method.str.startswith('rank_')].set_index('rank').loc[ranks]
    ax.plot(x, q.skill * 100, 'o-', color=DIR_COLORS[d], label=DIR_NAMES[d], ms=6, lw=2)
    ax.fill_between(x, q.lower95 * 100, q.upper95 * 100, color=DIR_COLORS[d], alpha=.13)
    a = t[(t.direction == d) & (t.method == 'original_rank16')].iloc[0]
    ax.plot(5.18, a.skill * 100, 'D', mfc='white', mec=DIR_COLORS[d], ms=7, zorder=5)
    ax.annotate(f'{q.skill.iloc[-1]*100:.1f}%', (6, q.skill.iloc[-1] * 100), xytext=(8, -3), textcoords='offset points', fontsize=8.5, color='#333')
ax.plot([], [], 'D', mfc='white', mec='.3', label='Benchmark model (rank 16)')
ax.set_xticks(x, ['0\naffine', '1', '2', '4', '8', '16', '32']); ax.set_xlim(-.3, 6.8); ax.set_ylim(0, 16)
ax.set_xlabel('Preset cross-gene correction rank'); ax.set_ylabel('Skill against destination mean (%)')
ax.legend(frameon=False, fontsize=8.5, loc='upper left')

# D: p53 minus other gain over affine
heading(outer[0, 1], 'B', 'Extra gain over affine: p53 targets minus other readouts')
ax = fig.add_subplot(outer[0, 1]); rk = [1, 2, 4, 8, 16, 32]; xx = np.arange(6)
for j, d in enumerate(dirs):
    q = con[con.direction == d].set_index('rank').loc[rk]; xs = xx + (j - .5) * .22
    ax.vlines(xs, q.lower99_583333 * 100, q.upper99_583333 * 100, color=DIR_COLORS[d], lw=1.8)
    ax.plot(xs, q.estimate * 100, 'o' if j == 0 else 's', color=DIR_COLORS[d], ms=6.5, label=DIR_NAMES[d], zorder=3)
    v = q.estimate.iloc[-1] * 100
    ax.annotate(f'{v:+.2f}', (xs[-1], v), xytext=(10 if j else -10, 0), textcoords='offset points', ha='left' if j else 'right', va='center', fontsize=8.5, color=DIR_COLORS[d])
ax.axhline(0, color='.5', ls='--', lw=.8)
ax.set_xticks(xx, [str(r) for r in rk]); ax.set_xlabel('Correction rank'); ax.set_ylabel('Gain difference (percentage points)')
ax.legend(frameon=False, fontsize=8.5, loc='upper left')


g=pd.read_csv(P/'inputs/ribosome_p53/all_group_skills.tsv',sep='\t')
pc=pd.read_csv(P/'inputs/ribosome_p53/primary_contrasts.tsv',sep='\t')
meta_r=json.loads((P/'inputs/ribosome_p53/SUMMARY.json').read_text());meta=meta_r
ax = fig.add_subplot(outer[1,0]); labels = []; y = 0
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
ax.set_title('C   Readout skill by target group', loc='left', fontweight='bold', fontsize=10, pad=10)

ax = fig.add_subplot(outer[1,1])
for j, d in enumerate(dirs):
    r = pc[pc.direction == d].iloc[0]
    ax.vlines(j, r.lower97_5 * 100, r.upper97_5 * 100, color=DIR_COLORS[d], lw=2)
    ax.plot(j, r.estimate * 100, 'o' if j == 0 else 's', color=DIR_COLORS[d], ms=7)
    ax.annotate(f'{r.estimate*100:+.2f}', (j, r.estimate * 100), xytext=(9, 0), textcoords='offset points', va='center', fontsize=8.5, color=DIR_COLORS[d])
ax.axhline(0, color='.5', ls='--', lw=.8); ax.set_xticks([0, 1], [DIR_NAMES[d] for d in dirs], fontsize=8.5); ax.set_xlim(-.5, 1.7)
ax.set_ylabel('Ribosome − other targets:\npreferential p53 gain (points)')
ax.set_title('D   Original rank-16 model', loc='left', fontweight='bold', fontsize=10, pad=10)

(P / 'figures').mkdir(exist_ok=True)
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(P / 'figures' / f'Figure_4_capacity_and_p53_readouts.{ext}', dpi=240, facecolor='white', bbox_inches='tight', pad_inches=.2)
plt.close(fig)
print('Figure 4 complete: rank capacity and p53 readout panels')

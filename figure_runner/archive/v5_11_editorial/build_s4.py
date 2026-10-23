"""Figure S4 (v5.11): exploratory destination-program readout follow-ups in chemical transfer.
A: Hallmark estrogen-response readouts toward ER-positive MCF7, against all other readouts and against
   expression/variability-matched readouts (family of 4; 98.75% intervals).
B: Cancer-pathway readouts against all other readouts (original, 98.75%) and against matched readouts (97.5%).
Open markers: matching judged inadequate (sensitivity only). Saved estimates; no refitting."""
from pathlib import Path
import numpy as np, pandas as pd
from _style import plt, DIR_COLORS, DIR_NAMES

P = Path(__file__).resolve().parent
e = pd.read_csv(P / 'inputs/readout_followups/estrogen_primary_endpoints.tsv', sep='\t')
dm = pd.read_csv(P / 'inputs/readout_followups/drug_matched_endpoints.tsv', sep='\t')
dp = pd.read_csv(P / 'inputs/drug_pathway_skill/primary_endpoints.tsv', sep='\t')
assert len(e) == 4 and len(dm) == 2

fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw=dict(wspace=.75))
def row(ax, y, v, lo, hi, col, filled, label=None):
    ax.hlines(y, lo, hi, color=col, lw=1.8); ax.plot([lo, hi], [y, y], '|', color=col, ms=7)
    ax.plot(v, y, 'o', color=col, mfc=col if filled else 'white', ms=7, zorder=3, label=label)
    ax.annotate(f'{v:+.2f}', (v, y), xytext=(0, 8), textcoords='offset points', ha='center', fontsize=8, color=col)

ax = axs[0]; c = DIR_COLORS['PC3_to_MCF7']; labs = []; y = 0
for gs, short in [('HALLMARK_ESTROGEN_RESPONSE_EARLY', 'Estrogen early'), ('HALLMARK_ESTROGEN_RESPONSE_LATE', 'Estrogen late')]:
    for ref in ['complement', 'matched']:
        r = e[(e.gene_set == gs) & (e.reference == ref)].iloc[0]
        row(ax, y, r.estimate * 100, r.lower98_75 * 100, r.upper98_75 * 100, c, r.matching_status != 'INADEQUATE_MATCHING_SENSITIVITY_ONLY')
        labs.append(f'{short} ({int(r.n_readouts)})\nvs {"all other" if ref == "complement" else "matched"} readouts'); y += 1
ax.set_yticks(range(4), labs, fontsize=8.5); ax.set_ylim(3.6, -.6); ax.axvline(0, color='.5', ls='--', lw=.8); ax.set_xlim(-.5, 4)
ax.set_xlabel('Skill difference (percentage points)'); ax.set_title('A   Estrogen-response readouts, PC3 → MCF7', loc='left', fontweight='bold', fontsize=10, pad=10)

ax = axs[1]; labs = []; y = 0
for d, dis, short in [('MCF7_to_PC3', 'prostate_cancer', 'Prostate cancer (23)'), ('PC3_to_MCF7', 'breast_cancer', 'Breast cancer (33)')]:
    c = DIR_COLORS[d]
    r = dp[(dp.direction == d) & (dp.endpoint == 'pathway_minus_complement_skill')].iloc[0]
    row(ax, y, r.estimate * 100, r.lower98_75 * 100, r.upper98_75 * 100, c, True); labs.append(f'{short}, {DIR_NAMES[d]}\nvs all other readouts'); y += 1
    m = dm[dm.direction == d].iloc[0]
    row(ax, y, m.matched_skill_difference * 100, m.lower97_5 * 100, m.upper97_5 * 100, c, m.matching_status != 'INADEQUATE_MATCHING_SENSITIVITY_ONLY')
    labs.append(f'{short}, {DIR_NAMES[d]}\nvs matched readouts'); y += 1
ax.set_yticks(range(4), labs, fontsize=8.5); ax.set_ylim(3.6, -.6); ax.axvline(0, color='.5', ls='--', lw=.8); ax.set_xlim(-2, 4)
ax.set_xlabel('Skill difference (percentage points)'); ax.set_title('B   Cancer-pathway readouts with matched references', loc='left', fontweight='bold', fontsize=10, pad=10)
ax.plot([], [], 'o', color='.3', label='Matching adequate or not applicable'); ax.plot([], [], 'o', color='.3', mfc='white', label='Matching inadequate (sensitivity only)')
fig.legend(frameon=False, fontsize=8, loc='upper left', bbox_to_anchor=(.30, -.02), ncol=2)
(P / 'figures').mkdir(exist_ok=True)
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(P / 'figures' / f'Figure_S4_destination_program_readouts.{ext}', dpi=240, facecolor='white', bbox_inches='tight', pad_inches=.2)
plt.close(fig); print('Figure S4 complete')

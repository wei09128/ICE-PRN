"""Shared style for ICE-PRN v5 figures. Palettes checked with the dataviz validator
(lightness, chroma, colour-vision separation and contrast all pass on white)."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DIR_COLORS = {'K562_to_RPE1': '#1b6fa8', 'RPE1_to_K562': '#c55f2c',
              'MCF7_to_PC3': '#1b6fa8', 'PC3_to_MCF7': '#c55f2c'}
DIR_NAMES = {'K562_to_RPE1': 'K562 → RPE1', 'RPE1_to_K562': 'RPE1 → K562',
             'MCF7_to_PC3': 'MCF7 → PC3', 'PC3_to_MCF7': 'PC3 → MCF7'}
GROUP_COLORS = ['#008f8a', '#b97a22', '#5b72b5', '#9a9a9a']   # rescued, persistent, both, other (neutral)
GROUP_TICKS = ['Rescued', 'Persistent\nfailure', 'Both\nsucceeded', 'Other']
INK = '#222222'; MUTED = '#666666'; GRID = '#bbbbbb'

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.edgecolor': '#444444', 'axes.labelcolor': INK,
                     'xtick.color': '#444444', 'ytick.color': '#444444', 'svg.fonttype': 'none',
                     'pdf.fonttype': 42, 'axes.titlesize': 10})


def panel(ax, letter, title):
    ax.set_title(f'{letter}   {title}', loc='left', fontweight='bold', fontsize=10, pad=10, color=INK)


def save(fig, stem):
    out = Path(__file__).resolve().parent / 'figures'
    out.mkdir(exist_ok=True)
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(out / f'{stem}.{ext}', dpi=300, facecolor='white')
    plt.close(fig)

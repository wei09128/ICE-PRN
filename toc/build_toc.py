"""Rebuild the ICE-PRN TOC graphic (PNG, editable SVG, vector PDF).

Requires matplotlib. Artwork size: 3.25 x 1.75 inches; PNG: 600 dpi.
Numbers summarize original held-out MSE skill against the destination mean.
"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'DejaVu Sans', 'svg.fonttype':'none',
                     'pdf.fonttype':42, 'font.size':7})
fig, ax = plt.subplots(figsize=(3.25,1.75))
fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
ax.set(xlim=(0,3.25),ylim=(0,1.75)); ax.axis('off')
ink='#203044'; blue='#2166A5'; teal='#137F78'; muted='#526172'
ax.text(.13,1.59,'What transfers between cell lines?',fontsize=9,
        weight='bold',color=ink,va='center')
def box(x,y,w,h,fill):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.018,rounding_size=0.05',
                               facecolor=fill,edgecolor='none'))
box(.13,.56,.83,.83,'#EFF3F7')
ax.text(.22,1.25,'Source response',weight='bold',color=ink,fontsize=6)
ax.text(.22,1.02,'Genetic',weight='bold',color=blue,fontsize=7)
ax.text(.22,.87,'K562 ↔ RPE1',fontsize=6.6,color=ink)
ax.text(.22,.70,'Chemical',weight='bold',color=teal,fontsize=7)
ax.text(.22,.57,'MCF7 ↔ PC3',fontsize=6.6,color=ink)
for y in (1.0,.69):
    ax.add_patch(FancyArrowPatch((1.00,y),(1.29,y),arrowstyle='-|>',mutation_scale=9,
                                linewidth=1.1,color=muted))
box(1.32,.56,1.80,.83,'#F0F7F6')
ax.text(1.43,1.25,'Learned destination prediction',weight='bold',color=ink,fontsize=6.6)
ax.text(1.44,.99,'10–12%',weight='bold',fontsize=12,color=blue,va='center')
ax.text(2.24,.99,'genetic skill',fontsize=7,color=ink,va='center')
ax.text(1.44,.69,'3.2–3.6%',weight='bold',fontsize=12,color=teal,va='center')
ax.text(2.39,.69,'chemical skill',fontsize=6.5,color=ink,va='center')
ax.text(1.33,.415,'MSE reduction vs destination mean',fontsize=6.6,color=muted)
ax.plot([.13,3.12],[.33,.33],color='#D7DEE5',lw=.6)
ax.text(.13,.205,'p53 association weakens after baseline-expression matching',
        fontsize=6.6,color=muted)
for ext in ('png','svg','pdf'):
    fig.savefig(OUT/f'ICE_PRN_TOC_graphic.{ext}',dpi=600,facecolor='white')
plt.close(fig)

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

weight = 5


def removeSpines(ax, loc):
    for pos in loc:
        ax.spines[pos].set_visible(False)


def setAxColor(ax, col, main=False):
    ax.yaxis.label.set_color(col)
    ax.tick_params(axis='y', colors=col)
    ax.spines['top'].set_visible(False)
    for l in ['left', 'right', 'bottom']:
        ax.spines[l].set(linewidth=1.5)
    if main:
        ax.spines['left'].set_color(col)
        ax.spines['bottom'].set_color((0, 0, 0))
    else:
        ax.spines['right'].set_color(col)


mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.family'] = 'CMU serif'

sizeFactor = 2
px = 1.0 / plt.rcParams['figure.dpi']
plt.rcParams['figure.figsize'] = [700*px*sizeFactor, 500*px*sizeFactor]
mpl.rcParams['font.size'] = 18*sizeFactor
plt.rcParams['lines.linewidth'] = 6
legendFontSize = 14*sizeFactor

color = {}
color['b'] = (0/255,   0/255,   0/255)
color['c1'] = (255/255, 225/255, 137/255)
color['c2'] = (226/255, 143/255,  83/255)
color['c3'] = (72/255, 133/255, 191/255)
color['c4'] = (181/255,  88/255,  91/255)
color['c5'] = (169/255, 209/255, 142/255)
color['w'] = (255/255, 255/255, 255/255)
color['g1'] = (240/255, 240/255, 240/255)
color['g2'] = (160/255, 160/255, 160/255)

wstr = '{:02d}'.format(weight)

df = pd.read_csv('{}/_log/results.csv'.format(wstr))

fig0, ax0 = plt.subplots()
fig0.subplots_adjust(right=0.75)

ax1 = ax0.twinx()
ax2 = ax0.twinx()

ax2.spines.right.set_position(("axes", 1.2))

p1, = ax0.plot(
    df['iter'],
    df[' vfrac'],
    color=color['b'],
    linestyle="-",
    linewidth=2.5,
)

p2, = ax1.plot(
    df['iter'],
    df[' comp'],
    color=color['c2'],
    linestyle="-",
    linewidth=2.5,
)

p3, = ax2.plot(
    df['iter'],
    df[' enrr'],
    color=color['c3'],
    linestyle="-",
    linewidth=2.5,
)

ax0.set_xlim(0, 350)
ax0.set_ylim(0, 1)
ax0.set_xticks([0, 50, 100, 150, 200, 250, 300, 350])
ax0.set_yticks([0, 0.25, 0.5, 0.75, 1.0])

ax0.set_xlabel('Design iterations')
ax0.set_ylabel('Volume fraction')
ax1.set_ylabel('Mean compliance [$\mathrm{J}$]')
ax2.set_ylabel('Energy release rate [$\mathrm{J} / \mathrm{mm}^2$]')

setAxColor(ax0, p1.get_color(), True)
setAxColor(ax1, p2.get_color())
setAxColor(ax2, p3.get_color())

plt.savefig('{}.pdf'.format(wstr))

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

plot = 2

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

enrr = []
comp = []
for i in range(6+1):
    istr = '{:02d}'.format(i)
    df = pd.read_csv('{}/_log/results.csv'.format(istr))
    enrr.append(df[' enrr'].iloc[-1])
    comp.append(df[' comp'].iloc[-1])

fig0, ax0 = plt.subplots()

if plot == 0:
    p1, = ax0.plot(
        [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
        enrr,
        color=color['b'],
        linestyle="-",
        marker="o",
        linewidth=2.5,
        markersize=10
    )
    ax0.set_ylabel('Energy release rate [$\mathrm{J} / \mathrm{mm}^2$]')
    ax0.set_xlabel('Weight $\\alpha$')
elif plot == 1:
    p1, = ax0.plot(
        [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
        comp,
        color=color['b'],
        linestyle="-",
        marker="o",
        linewidth=2.5,
        markersize=10
    )
    ax0.set_ylabel('Mean compliance [$\mathrm{J}$]')
    ax0.set_xlabel('Weight $\\alpha$')
elif plot == 2:
    p1, = ax0.plot(
        enrr,
        comp,
        color=color['b'],
        linestyle="-",
        marker="o",
        linewidth=2.5,
        markersize=10
    )
    ax0.set_xlabel('Energy release rate [$\mathrm{J} / \mathrm{mm}^2$]')
    ax0.set_ylabel('Mean compliance [$\mathrm{J}$]')

# ax0.set_xlim(0, 0.6)

ax0.set_ylabel('Energy release rate [$\mathrm{J} / \mathrm{mm}^2$]')

ax0.yaxis.label.set_color(p1.get_color())
ax0.tick_params(axis='y', colors=p1.get_color())
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
for l in ['left', 'bottom']:
    ax0.spines[l].set(linewidth=1.5)
ax0.spines['left'].set_color(p1.get_color())
ax0.spines['bottom'].set_color((0, 0, 0))

if plot == 0:
    plt.savefig('enrr.pdf')
elif plot == 1:
    plt.savefig('comp.pdf')
elif plot == 2:
    plt.savefig('pareto.pdf')

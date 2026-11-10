"""Draw the hadron production/decay tree from pp_single_event_full.log."""
from pathlib import Path
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.path import Path as PlotPath

folder = Path(__file__).resolve().parent
text = (folder / 'pp_single_event_full.log').read_text()
text = text.split('Complete event history')[1]
particles = {}
for line in text.splitlines():
    fields = line.split()
    if len(fields) != 15 or not re.fullmatch(r'\d+', fields[0]):
        continue
    i = int(fields[0])
    particles[i] = {'pid': int(fields[1]), 'name': fields[2].strip('()'),
                    'status': int(fields[3]), 'm1': int(fields[4]), 'm2': int(fields[5]),
                    'd1': int(fields[6]), 'd2': int(fields[7]), 'E': float(fields[13])}

# In this event, status 83/84 marks primary string-fragmentation hadrons.
# Collapse all parton/string history into one shared hadronization node;
# only draw explicit hadron decay relations below those primary hadrons.
roots = [i for i, p in particles.items() if abs(p['status']) in (83, 84)]
children = {}
for i, p in particles.items():
    children[i] = [j for j, q in particles.items() if q['m1'] == i and q['m2'] == 0]

def leaf_count(i):
    if not children[i]:
        return 1
    return sum(leaf_count(j) for j in children[i])

labels = {22: r'$\gamma$', 111: r'$\pi^0$', 211: r'$\pi^+$', -211: r'$\pi^-$',
          2212: r'$p$', 2112: r'$n$', -2112: r'$\bar n$', -2212: r'$\bar p$',
          321: r'$K^+$', -321: r'$K^-$', 221: r'$\eta$', 331: r"$\eta'$",
          223: r'$\omega$', 113: r'$\rho^0$', 213: r'$\rho^+$', -323: r'$K^{*-}$'}
N_leaves = sum(leaf_count(i) for i in roots)
fig, ax = plt.subplots(figsize=(19, max(12, N_leaves * 0.43)))
positions = {}
next_y = 0

def place(i, depth):
    global next_y
    if not children[i]:
        y = next_y
        next_y += 1
    else:
        ys = [place(j, depth + 1) for j in children[i]]
        y = np.mean(ys)
    positions[i] = (1.9 + depth * 3.25, y)
    return y

for i in roots:
    place(i, 0)
origin = (-1.0, (N_leaves - 1) / 2)

def draw_edge(start, i, curved=False):
    p = particles[i]
    end = positions[i]
    photon = p['pid'] == 22
    color = '#bf8500' if photon else 'black'
    neutral = p['pid'] in (22, 111, 2112, -2112, 221, 331, 223, 113)
    width = 0.8 + 1.5 * np.sqrt(p['E'] / 200.0)
    if curved:
        # Bend the production branches toward separate horizontal endpoints,
        # leaving particle labels in vertically separated rows.
        path = PlotPath([start, (start[0] + 0.05, end[1]),
                         (end[0] - 1.0, end[1]), end],
                        [PlotPath.MOVETO, PlotPath.CURVE4, PlotPath.CURVE4, PlotPath.CURVE4])
        arrow = FancyArrowPatch(path=path, arrowstyle='-|>', mutation_scale=10,
                                linewidth=width, linestyle='--' if neutral else '-', color=color)
        x = end[0] - 0.85
        y = end[1]
    else:
        arrow = FancyArrowPatch(start, end, arrowstyle='-|>', mutation_scale=10,
                                linewidth=width, linestyle='--' if neutral else '-', color=color)
        x = 0.4 * start[0] + 0.6 * end[0]
        y = 0.4 * start[1] + 0.6 * end[1]
    ax.add_patch(arrow)
    label = labels.get(p['pid'], p['name'])
    ax.text(x, y - 0.16, f"{label}  {p['E']:.3g} GeV  [{i}]", fontsize=9,
            ha='center', va='bottom', color=color,
            bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.9, 'pad': 1})
    if children[i]:
        ax.plot(*end, 'ko', markersize=3)
        for j in children[i]:
            draw_edge(end, j)

for i in roots:
    draw_edge(origin, i, True)
for dy, label in [(-1.4, r'$p$: $T=1000$ GeV'), (1.4, r'$p$: stationary target')]:
    start = (-3.9, origin[1] + dy)
    ax.add_patch(FancyArrowPatch(start, origin, arrowstyle='-|>', mutation_scale=18,
                                linewidth=2.5, color='black'))
    ax.text(start[0], start[1] - 0.25, label, fontsize=11, va='bottom')
ax.plot(*origin, 'ko', markersize=6)
ax.text(-1.15, origin[1] + 2.5, 'pp interaction\n+ hadronization', rotation=90,
        ha='right', va='center', fontsize=11)
ax.set_xlim(-4.2, max(x for x, y in positions.values()) + 0.5)
ax.set_ylim(N_leaves - 0.2, -1.4)
ax.axis('off')
ax.set_title('One PYTHIA pp event: hadron production and decays\n'
             '1 TeV projectile kinetic energy, stationary proton; non-diffractive, sqrt(s) = 43.36 GeV',
             fontsize=14, pad=28)
fig.text(0.5, 0.012,
         'Labels: LAB total energy; [index] = complete-event record. Gold: photons. Dashed: neutral particles.\n'
         'Schematic connections, not spatial tracks. Parton history collapsed. Terminal particles follow enabled PYTHIA decays.\n'
         'Values read from the rounded text log; not all terminal particles are stable on astrophysical timescales.',
         ha='center', fontsize=10)
fig.subplots_adjust(left=0.03, right=0.99, top=0.94, bottom=0.07)
fig.savefig(folder / 'pp_event_channels.png', dpi=180)
fig.savefig(folder / 'pp_event_channels.svg')

# Check that each final particle is included exactly once in the drawn tree.
terminal = [i for i in positions if not children[i]]
expected = [i for i, p in particles.items() if p['status'] > 0]
assert sorted(terminal) == sorted(expected)
print('Primary hadrons:', len(roots), 'Final particles drawn:', len(terminal))
print('Saved pp_event_channels.png and pp_event_channels.svg')

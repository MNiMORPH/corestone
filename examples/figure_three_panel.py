#! /usr/bin/python3
"""
The three-panel figure: the plumbing, the equation, and what is left.

Follows the causal chain rather than showing the answer alone. The middle panel
is the point of the figure: the affinity term is the only field in the model,
so putting it on screen is the difference between seeing the equation and being
told it.

    python3 examples/figure_three_panel.py            # writes the PNG
    python3 examples/figure_three_panel.py --show     # and opens it
"""
import sys

import numpy as np
import matplotlib
if "--show" not in sys.argv:
    matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from scipy.ndimage import uniform_filter

from corestone import (FractureNetwork, Weathering, orthogonal_grid,
                       uniform_grid_shape)

# ---- what to draw -------------------------------------------------------------
DX = 0.05                                  # cell size [m]
SPACING = 1.5                              # joint spacing [m]
# Cell counts that let the joint set tile the section exactly, with a joint on
# each wall and no odd block left over at the edges.
NZ, NX = uniform_grid_shape(20.0, 15.0, DX, SPACING)
KYR = 100.0                                # elapsed time
OUT = "examples/figure_three_panel.png"

net = FractureNetwork(NZ, NX, DX).seed(sets=orthogonal_grid(SPACING),
                                       rng=np.random.default_rng(12345))
model = Weathering(net).run(years=KYR * 1e3)

X = model.dissolved_fraction
q = model.q
affinity = model.affinity
L_eq = model.equilibration_length
LX, LZ = net.lx, net.lz

INK, MUTED = "#1a1a1a", "#6b6b6b"
EXT = [0.0, LX, LZ, 0.0]                   # depth increases downward

fig = plt.figure(figsize=(16.8, 6.9))
gs = GridSpec(2, 3, height_ratios=[1.0, 0.04], figure=fig,
              left=0.045, right=0.988, top=0.735, bottom=0.235,
              wspace=0.13, hspace=0.30)
axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
cbax = [fig.add_subplot(gs[1, i]) for i in range(3)]
for a in axes[1:]:
    a.sharey(axes[0])

CAPTIONS = []


def joints(ax, color="#1a1a1a", lw=0.9, alpha=0.6):
    """
    Draw the two sets distinguishably: throughgoing solid and heavier, abutting
    dashed. Told apart by dash pattern as well as weight, never by colour alone.
    """
    for (p0, p1), name in zip(net.segments, net.segment_set):
        through = name == "J1"
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color,
                lw=lw if through else lw * 0.75,
                alpha=alpha if through else alpha * 0.85,
                ls="-" if through else (0, (3.5, 2.2)),
                solid_capstyle="round", zorder=3)


def dress(ax, title, caption):
    ax.set_xlim(0, LX)
    ax.set_ylim(LZ, 0)
    ax.set_aspect("equal")
    ax.set_xlabel("Distance [m]", color=MUTED, fontsize=9, labelpad=2)
    ax.set_title(title, color=INK, fontsize=12.5, loc="left", pad=9,
                 weight="bold")
    for sp in ax.spines.values():
        sp.set_color("#d9d9d9")
    ax.tick_params(colors=MUTED, labelsize=8.5, length=3)
    CAPTIONS.append((ax, caption))


def bar(im, cax, label, ticks=None, ticklabels=None):
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label(label, color=MUTED, fontsize=9, labelpad=4)
    cb.ax.tick_params(colors=MUTED, labelsize=8, length=2, pad=2)
    cb.outline.set_visible(False)
    if ticks is not None:
        cb.set_ticks(ticks)
        if ticklabels is not None:
            cb.set_ticklabels(ticklabels)


# ---- 1: where the water goes -------------------------------------------------
ax = axes[0]
qn = q / q[0, :].mean()
im = ax.imshow(np.log10(np.maximum(qn, 1e-3)), extent=EXT, origin="upper",
               cmap="Blues", vmin=-1.5, vmax=1.5, interpolation="nearest")
joints(ax)
dress(ax, "1   The plumbing",
      "Steady Darcy flow: the head field decides where water goes. Joints\n"
      "conduct, and water runs ALONG the horizontal ones to the next vertical.")
ax.set_ylabel("Depth [m]", color=MUTED, fontsize=9)
bar(im, cbax[0], "water flux, relative to mean infiltration",
    ticks=[-1.5, 0, 1.5], ticklabels=["0.03×", "1×", "30×"])

# ---- 2: where it can still dissolve ------------------------------------------
ax = axes[1]
im = ax.imshow(affinity, extent=EXT, origin="upper", cmap="Greens",
               vmin=0, vmax=1, interpolation="nearest")
joints(ax, color="#0b3d20")
dress(ax, "2   The equation, made visible",
      "Dark green is fresh water with capacity left. Pale is water already\n"
      "saturated: the bracket is zero and weathering has stopped.")
bar(im, cbax[1], r"affinity term  $1 - C/C_{eq}$")
ax.text(0.982, 0.045, "$L_{eq}$ = %.2f m" % L_eq,
        transform=ax.transAxes, ha="right", va="bottom", fontsize=10.5,
        color=INK, bbox=dict(boxstyle="round,pad=0.35", fc="white",
                             ec="#c9c9c9", alpha=0.94))

# ---- 3: what is left ----------------------------------------------------------
ax = axes[2]
im = ax.imshow(X, extent=EXT, origin="upper", cmap="Oranges", vmin=0, vmax=1,
               interpolation="nearest")
ax.contour(np.linspace(DX / 2, LX - DX / 2, NX),
           np.linspace(DX / 2, LZ - DX / 2, NZ),
           X, levels=[model.x_grus], colors="#7a2e00", linewidths=1.4, zorder=4)
joints(ax, color="#3a1c00", lw=0.8, alpha=0.42)
dress(ax, "3   Corestones in grus",
      "Rock further from a joint than $L_{eq}$ is never reached by\n"
      "undersaturated water. Not tougher granite – unvisited granite.")
bar(im, cbax[2], "fraction of the soluble phase dissolved")

# Label a corestone genuinely surrounded by grus, and a grus patch, keeping the
# two labels on opposite halves so they cannot collide.
enclosed = uniform_filter(model.is_grus.astype(float),
                          size=max(int(round(0.35 / DX)), 3))
half = np.zeros_like(model.is_corestone)
half[:, :NX // 2] = True
core = np.unravel_index(
    np.argmax(np.where(model.is_corestone & half, enclosed, -1.0)), X.shape)
gi = np.unravel_index(
    np.argmax(np.where(model.is_grus & ~half, enclosed, -1.0)), X.shape)
for (iz, ix), label, off in ((core, "corestone", 4.2), (gi, "grus", -3.4)):
    x0, z0 = ix * DX + DX / 2, iz * DX + DX / 2
    ax.annotate(label, xy=(x0, z0),
                xytext=(np.clip(x0, 2.2, LX - 2.2),
                        np.clip(z0 + off, 1.0, LZ - 1.0)),
                ha="center", fontsize=10.5, color=INK, zorder=6,
                bbox=dict(boxstyle="round,pad=0.34", fc="white", ec="#a8a8a8"),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.2,
                                shrinkA=3, shrinkB=4))

ax.legend(handles=[Line2D([0], [0], color="#3a1c00", lw=1.0,
                          label="joint, throughgoing"),
                   Line2D([0], [0], color="#3a1c00", lw=0.8,
                          ls=(0, (3.5, 2.2)), label="joint, abutting"),
                   Line2D([0], [0], color="#7a2e00", lw=1.4,
                          label="grus threshold, X = %.1f" % model.x_grus)],
          loc="lower left", fontsize=8.5, frameon=True, framealpha=0.93,
          edgecolor="#d9d9d9", labelcolor=INK)

# ---- header and captions, placed from the final axes geometry ----------------
fig.text(0.045, 0.962, "corestone – fracture-controlled granite weathering",
         fontsize=15.5, color=INK, ha="left", va="top", weight="bold")
fig.text(0.045, 0.915,
         "%.2f × %.2f m section at %.0f cm, %.0f kyr at %.0f K, joints every "
         "%.2f m. Every parameter is a placeholder – see "
         "design/02-teaching-scope.md."
         % (LX, LZ, DX * 100, KYR, model.T, SPACING),
         fontsize=9.5, color=MUTED, ha="left", va="top")

p2 = axes[1].get_position()
fig.text(p2.x0 + p2.width / 2, 0.845,
         r"$R \;=\; k(T)\,A\,\left(1 - C/C_{eq}\right)$",
         fontsize=17, color="#0b3d20", ha="center", va="center")

for ax, caption in CAPTIONS:
    pos = ax.get_position()
    fig.text(pos.x0, 0.118, caption, fontsize=9.5, color=MUTED,
             ha="left", va="top", linespacing=1.55)

fig.savefig(OUT, dpi=170, facecolor="white")
print("grus %.1f %%   corestone %.1f %%   L_eq %.2f m"
      % (model.is_grus.mean() * 100, model.is_corestone.mean() * 100, L_eq))
print("wrote %s" % OUT)
if "--show" in sys.argv:
    plt.show()

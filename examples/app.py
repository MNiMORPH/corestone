#! /usr/bin/python3
"""
The interactive demo: four sliders and a Run button.

Serve it locally with

    PYTHONPATH=src panel serve examples/app.py --show

or compile it to a self-contained browser page with artesian. numpy, scipy and
matplotlib are all bundled by Pyodide -- checked with ``artesian check``, not
assumed -- so the model itself runs in the browser rather than on a server.

Why the angle and spacing sliders SNAP. A conjugate pair rotated by theta has
along-x periods ``S/cos(theta)`` and ``S/sin(theta)``. The section is periodic
left-to-right -- it has no side walls, because a no-flow wall manufactures a
domain-scale circulation -- and both periods have to divide the width for the
joints to close on themselves. That needs ``tan(theta) = b/a`` for integers,
with a and b dividing the cell count, and it quantises the spacing to
``lx / (k * hypot(a, b))``. Off those values the joints do not line up across
the seam. Measured: on the snapped set the distance-to-joint at the seam is
within 13 % of the interior, so there is effectively no defect.
"""
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib.figure import Figure
import panel as pn

from corestone import (FractureNetwork, Weathering, orthogonal_grid,
                       tiling_angles, tiling_spacing, YEAR)

pn.extension()

# ---- the section ------------------------------------------------------------
# Resolution is chosen for how many cells cross a BLOCK, which is what makes a
# corestone look round rather than stepped: spacing/dx. At dx = 0.05 m a 1.5 m
# joint spacing gives 30 cells per block and a 2 m spacing gives 40. A smaller
# section at finer resolution beats a larger one coarsely resolved -- two or
# three blocks is plenty to show the rounding.
#
#   section   dx     cells   cells/block   200 kyr
#   6 x 6 m  0.05    14520        18         12.8 s
#   4 x 4 m  0.05     6480        30          3.4 s
#   3 x 3 m  0.05     3660        30          1.6 s   <- here
#   3 x 3 m  0.02    22650        67         30.2 s
#
# The width also decides how many orientations tile: the lattice angles need
# a and b to divide the cell count, so 60 (with divisors 1..6) gives thirteen
# angles where 80 gives seven. Narrower is both faster AND better sampled.
DX = 0.05                       # cell size [m]
NX, NZ = 60, 61                 # 3.0 x 3.05 m
LX, LZ = NX * DX, NZ * DX

#: Rotations at which the pair tiles the periodic width. Low indices only --
#: high ones force an unusably fine spacing (0.1 m at 0.95 degrees).
ANGLES = tiling_angles(NX, max_index=6)
INK, MUTED = "#1a1a1a", "#6b6b6b"


def spacings_for(angle_deg):
    """The tiling spacings available at this angle, coarse to fine."""
    a, b = next((a, b) for ang, a, b in ANGLES if abs(ang - angle_deg) < 1e-6)
    out = []
    for k in range(1, 9):
        S = tiling_spacing(LX, a, b, LX / (k * np.hypot(a, b)))
        if 0.3 <= S <= 3.0 and not any(abs(S - v) < 1e-9 for v in out):
            out.append(S)
    return sorted(out, reverse=True) or [LX / np.hypot(a, b)]


# ---- widgets ----------------------------------------------------------------
angle = pn.widgets.DiscreteSlider(
    name="Joint orientation [degrees from vertical/horizontal]",
    options={"%.2f°" % a: a for a, _, _ in ANGLES}, value=0.0)
spacing = pn.widgets.DiscreteSlider(
    name="Joint spacing [m]  (isotropic)",
    options={"%.3f m" % s: s for s in spacings_for(0.0)},
    value=spacings_for(0.0)[0])
velocity = pn.widgets.FloatSlider(
    name="Infiltration rate [m/yr]", start=0.05, end=1.0, step=0.05, value=0.30)
elapsed = pn.widgets.IntSlider(
    name="Time elapsed [kyr]", start=0, end=200, step=5, value=50)
run = pn.widgets.Button(name="▶  Run", button_type="primary", width=140)
status = pn.pane.Markdown("", styles={"color": MUTED, "font-size": "90%"})
plot = pn.pane.Matplotlib(dpi=110, tight=False, width=980, height=440)


@pn.depends(angle.param.value, watch=True)
def _resnap_spacing(a):
    """The tiling spacings depend on the angle, so the options move with it."""
    opts = spacings_for(a)
    nearest = min(opts, key=lambda s: abs(s - spacing.value))
    spacing.options = {"%.3f m" % s: s for s in opts}
    spacing.value = nearest


def simulate():
    net = FractureNetwork(NZ, NX, DX, periodic_x=True).seed(
        sets=orthogonal_grid(spacing.value, rotation=angle.value),
        rng=np.random.default_rng(12345))
    m = Weathering(net)
    m.set_infiltration(velocity.value / YEAR)
    if elapsed.value > 0:
        m.run(years=elapsed.value * 1e3)
    else:
        m.initialize()
        m.c = m.solve_solute(m.reaction_coefficient)
    return net, m


def draw(net, m):
    fig = Figure(figsize=(9.0, 4.0), layout="constrained")
    extent = [0.0, LX, LZ, 0.0]                     # depth increases downward
    for ax, field, cmap, title, label in (
            (fig.add_subplot(1, 2, 1), m.affinity, "Greens",
             "Where the water can still dissolve", r"$1 - C/C_{eq}$"),
            (fig.add_subplot(1, 2, 2), m.dissolved_fraction, "Oranges",
             "What is left of the rock", "fraction dissolved")):
        im = ax.imshow(field, extent=extent, origin="upper", cmap=cmap,
                       vmin=0, vmax=1, interpolation="nearest")
        for (p0, p1), name in zip(net.segments, net.segment_set):
            ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="#2a2a2a",
                    lw=0.8 if name == "J1" else 0.6,
                    ls="-" if name == "J1" else (0, (3, 2)), alpha=0.45)
        ax.set_xlim(0, LX)
        ax.set_ylim(LZ, 0)
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=10, color=INK, loc="left")
        ax.set_xlabel("Distance [m]", fontsize=8, color=MUTED)
        ax.tick_params(labelsize=7, colors=MUTED)
        for sp in ax.spines.values():
            sp.set_color("#d9d9d9")
        cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
        cb.set_label(label, fontsize=8, color=MUTED)
        cb.ax.tick_params(labelsize=7, colors=MUTED)
    fig.axes[0].set_ylabel("Depth [m]", fontsize=8, color=MUTED)
    return fig


def go(event=None):
    run.disabled = True
    status.object = "running…"
    t0 = time.perf_counter()
    net, m = simulate()
    plot.object = draw(net, m)
    status.object = (
        "**%.0f kyr** · joints every **%.3f m** at **%.2f°** · "
        "infiltration **%.2f m/yr**  \n"
        "grus %.1f %% · corestone %.1f %% · saturation length %.2f m at mean "
        "infiltration  \n*%.1f s*"
        % (elapsed.value, spacing.value, angle.value, velocity.value,
           100 * m.is_grus.mean(), 100 * m.is_corestone.mean(),
           m.saturation_length, time.perf_counter() - t0))
    run.disabled = False


run.on_click(go)
go()

app = pn.Column(
    pn.pane.Markdown(
        "## corestone — fracture-controlled granite weathering\n"
        r"$$\nabla\cdot(qc) - \nabla\cdot(D\nabla c) = r\,(1-c)$$"
        "\n\nWater runs down the joints and dissolves the rock it touches. "
        "Where it has equilibrated it stops, however soluble the rock. "
        "Blocks weather inward from every face and fastest at the corners, "
        "which shed solute to two joints rather than one — that is the "
        "rounding. **Every parameter is a placeholder.**"),
    pn.Row(pn.Column(angle, spacing, velocity, elapsed,
                     pn.Row(run), status, width=380),
           plot),
)
app.servable()

"""
Fracture-controlled granite weathering, in the browser.

The exercise this serves has one payload, and it is a misconception:

    A CORESTONE IS NOT TOUGHER ROCK.

Same granite, same minerals, same temperature as the crumbling grus around it.
It survives because the water never reached it, or reached it already carrying
all the solute it can hold. Weathering here is a race between how fast joints
deliver fresh water and how fast rock dissolves into it, and the reader runs
that race by pressing play.

What is drawn, left and right:

* **Where the water can still dissolve**, ``1 - C/C_eq``. Dark green is hungry
  water, white is water at saturation. Watch it: fresh water enters at the
  surface, runs down the joints, and saturates as it goes.
* **What is left of the rock**, the fraction dissolved. The joints go first,
  then the block faces, then the corners -- a corner sheds solute to two
  joints where a face sheds to one -- and what is left in the middle is a
  corestone.

Build and view::

    artesian build interactive_demo/corestone_panel.py -o _build \\
        -p . -r numpy -r scipy --serve

**Every parameter in the model is a placeholder.** None is measured. The
demo teaches the mechanism, and no number out of it is a result.
"""
import numpy as np
import panel as pn
from bokeh.models import ColorBar, ColumnDataSource, LinearColorMapper, Range1d
from bokeh.palettes import Greens256, Oranges256
from bokeh.plotting import figure

from artesian.live import animator, reset_button, responsive

from corestone import (FractureNetwork, Weathering, orthogonal_grid,
                       tiling_angles, tiling_spacings, YEAR)

pn.extension()

# ---- the section ------------------------------------------------------------
# 3.0 x 3.05 m at 5 cm cells. Resolution is chosen for how many cells cross a
# BLOCK, since that is what makes a corestone look round rather than stepped:
# a 1 m joint spacing is 20 cells across. A small section finely resolved beats
# a large one coarsely resolved, and two or three blocks is plenty to show the
# rounding.
DX = 0.05                       # cell size [m]
NX, NZ = 60, 61                 # 3.0 x 3.05 m
LX, LZ = NX * DX, NZ * DX

#: Rotations at which the joint pair tiles the periodic width exactly. The
#: index cap is not arbitrary: a high-index angle tiles only at a very fine
#: spacing (nine divisions of a three-metre section), so its spacing slider
#: would carry a single choice. Capping at 4 keeps the angles that offer a
#: real range of spacings -- 0, 14.0, 18.4, 26.6, 33.7, 36.9 and 45 degrees.
ANGLES = tiling_angles(NX, max_index=4)

#: Which snapped spacings to offer. Below 0.3 m a block is 6 cells across and
#: looks square however long it runs; above 3 m there is no block inside the
#: section at all.
SPACING_LOW, SPACING_HIGH = 0.3, 3.0

#: Stop here. 200 kyr dissolves the section at the default settings, so there
#: is nothing further to watch.
END_KYR = 200.0

#: Tighter than the model's own default of 0.03, for two reasons that happen
#: to agree. One frame is one step, so the budget sets how long the animation
#: lasts: 0.03 is 58 frames, under two seconds at 30 fps, which is over before
#: a reader has focused on it. 0.01 is 180 frames, about six seconds. It also
#: cuts the error by roughly three, since the error is close to linear in this.
#: A step costs 2.9 ms here, so six seconds of animation is 0.5 s of arithmetic
#: and the frame budget is nowhere near threatened.
C_DRIFT_MAX = 0.01

# Lay the app out to look right at this width; the embedding page scales the
# whole thing above it, so everything enlarges together rather than the figures
# growing while the sliders stay 18 px tall.
#
# It is also a HARD CAP on the layout, and that is not decoration. Left
# uncapped, a stretch-to-fit app inside an iframe that sizes itself to its
# content is a feedback loop with no fixed point: the app is as wide as it is
# given, the frame is as wide as the app, and nothing settles it. Desktop
# browsers pin the frame and hide the problem; iOS Safari and every iPad
# browser (they are all WebKit) size an iframe to its content, and the demo
# ran away wider than the page. Capping here gives the loop a fixed point,
# from the app's side, in one place, whatever the browser does.
DESIGN_WIDTH = 900
#: A third of the width each, since the three sit on one row. Stacked, each
#: took a label line and a track line and the three of them were 150 px of an
#: 809 px app -- and the embedding page scales that height along with the
#: width, so every pixel here is multiplied on a wide screen.
SLIDER_WIDTH = DESIGN_WIDTH // 3 - 20
#: Wider than it is tall, because a figure is not its data area: the depth
#: axis and its label take about 55 px on the left and the colour bar another
#: 60 on the right, while only the distance axis (~55 px) is below. Sized 1:1,
#: as this was, a SQUARE section comes out visibly taller than it is wide.
FIG_W, FIG_H = 460, 400


def _spacings(angle_deg):
    """The snapped spacings available at this angle, coarse to fine."""
    a, b = next((a, b) for ang, a, b in ANGLES if abs(ang - angle_deg) < 1e-6)
    return tiling_spacings(LX, a, b, SPACING_LOW, SPACING_HIGH)


# ---- widgets ----------------------------------------------------------------
# The angle and spacing sliders SNAP. A joint pair rotated by theta has
# along-x periods S/cos(theta) and S/sin(theta), and the section is periodic
# left-to-right -- it has no side walls, because a no-flow wall manufactures a
# domain-scale circulation. Both periods must divide the width for the joints
# to close on themselves, which needs tan(theta) = b/a for integers and
# quantises the spacing too. Off those values the joints fail to line up
# across the seam.
angle = pn.widgets.DiscreteSlider(
    name="Joint orientation [°]",
    options={"%.1f°" % a: a for a, _, _ in ANGLES}, value=0.0,
    sizing_mode="stretch_width", max_width=SLIDER_WIDTH)
spacing = pn.widgets.DiscreteSlider(
    name="Joint spacing [m]",
    options={"%.2f m" % s: s for s in _spacings(0.0)}, value=1.0,
    sizing_mode="stretch_width", max_width=SLIDER_WIDTH)
infiltration = pn.widgets.FloatSlider(
    name="Infiltration [m/yr]", start=0.05, end=1.00, step=0.05,
    value=0.30, format="0.00",
    sizing_mode="stretch_width", max_width=SLIDER_WIDTH)


def _build():
    """A fresh network and a fresh model at the current slider settings."""
    net = FractureNetwork(NZ, NX, DX, periodic_x=True).seed(
        sets=orthogonal_grid(spacing.value, rotation=angle.value),
        rng=np.random.default_rng(12345))
    m = Weathering(net)
    m.set_infiltration(infiltration.value / YEAR)
    m.c_drift_max = C_DRIFT_MAX
    m.initialize()
    m.c = m.solve_solute(m.reaction_coefficient)
    return net, m


# `sim`, never `state`: panel exports pn.state, and shadowing it fails silently.
sim = {}


def step():
    """One weathering step per frame, and stop at the end of the run."""
    m = sim["model"]
    target = END_KYR * 1e3 * YEAR
    if m.t >= target - 1e-9 * YEAR:
        run.value = False                      # reached the end; pause
        return
    m.update(dt_limit=target - m.t)
    _redraw()


def do_reset():
    """Rebuild from the sliders. Every slider here is structural."""
    net, m = _build()
    sim["net"], sim["model"] = net, m
    _joints()
    _redraw()


def _joints():
    """The joint traces, as segments in metres."""
    net = sim["net"]
    seg = np.array([[p0[0], p0[1], p1[0], p1[1]]
                    for p0, p1 in net.segments]) if net.segments else \
        np.zeros((0, 4))
    for src in (joints_left, joints_right):
        src.data = {"x0": seg[:, 0], "y0": seg[:, 1],
                    "x1": seg[:, 2], "y1": seg[:, 3]}


def _redraw():
    m = sim["model"]
    affinity.data = {"image": [m.affinity]}
    dissolved.data = {"image": [m.dissolved_fraction]}
    fig_left.title.text = "Where the water can still dissolve"
    fig_right.title.text = "What is left of the rock"
    readout.object = (
        "**%.0f kyr** &nbsp;·&nbsp; grus **%.0f %%** &nbsp;·&nbsp; "
        "corestone **%.0f %%**"
        % (m.t / YEAR / 1e3, 100 * m.is_grus.mean(),
           100 * m.is_corestone.mean()))


# ---- figures ----------------------------------------------------------------
# Depth increases DOWNWARD, so the y range runs from LZ at the bottom of the
# axis to 0 at the top. Row 0 of the array is the ground surface, and bokeh
# draws row 0 at the anchor and later rows at increasing y, which with a
# reversed range puts the surface at the top where it belongs -- no flip.
affinity = ColumnDataSource(data={"image": [np.zeros((NZ, NX))]})
dissolved = ColumnDataSource(data={"image": [np.zeros((NZ, NX))]})
joints_left = ColumnDataSource(data={"x0": [], "y0": [], "x1": [], "y1": []})
joints_right = ColumnDataSource(data={"x0": [], "y0": [], "x1": [], "y1": []})


def _panel(source, joints, palette, label):
    """One map of the section, with its colour bar and joint traces."""
    fig = figure(width=FIG_W, height=FIG_H,
                 x_axis_label="Distance [m]", y_axis_label="Depth [m]",
                 x_range=Range1d(0, LX), y_range=Range1d(LZ, 0),
                 tools="", toolbar_location=None)
    mapper = LinearColorMapper(palette=palette, low=0.0, high=1.0)
    fig.image(image="image", x=0, y=0, dw=LX, dh=LZ, source=source,
              color_mapper=mapper)
    fig.segment("x0", "y0", "x1", "y1", source=joints,
                color="#2a2a2a", line_width=1, alpha=0.45)
    bar = ColorBar(color_mapper=mapper, width=8, title=label,
                   label_standoff=6, padding=4)
    fig.add_layout(bar, "right")
    # Half the design width each, since the two sit side by side. Without a
    # bound, responsive() defaults to 1200 PER FIGURE, so the row is entitled
    # to 2400 -- far wider than the app is laid out for.
    responsive(fig, aspect_ratio=float(FIG_W) / FIG_H,
               max_width=DESIGN_WIDTH // 2)
    return fig


# Palettes reversed so that 0 is pale and 1 is saturated: bokeh's 256-step
# ramps run dark to light.
fig_left = _panel(affinity, joints_left, Greens256[::-1], "1 − C/Ceq")
fig_right = _panel(dissolved, joints_right, Oranges256[::-1],
                   "fraction dissolved")

readout = pn.pane.Markdown("", sizing_mode="stretch_width")
run = animator(step)
reset = reset_button(do_reset, name="Fresh rock")

# Every slider rebuilds: the joint geometry is the initial condition, and the
# infiltration rate sets a flow field that is solved once and held. None of
# them is a forcing that can be turned while the rock evolves, so changing one
# starts the clock again rather than pretending otherwise.
for w in (angle, spacing, infiltration):
    w.param.watch(lambda event: do_reset(), "value")


@pn.depends(angle.param.value, watch=True)
def _resnap_spacing(a):
    """Which spacings tile depends on the angle, so the options move with it."""
    opts = _spacings(a)
    nearest = min(opts, key=lambda s: abs(s - spacing.value))
    spacing.options = {"%.2f m" % s: s for s in opts}
    spacing.value = nearest


do_reset()

# ONE line, not four paragraphs. This app is embedded in a page that already
# explains the mechanism directly above the frame, so prose here is read twice
# and paid for once in height: the three paragraphs that used to sit here were
# 180 px of an 809 px app, and the frame is sized to its content. What has to
# stay is the placeholder warning, which belongs with the numbers rather than
# with the teaching, and enough of a title that the app still makes sense
# opened on its own.
pn.Column(
    pn.pane.Markdown(
        "**Why a corestone survives** — press **▶** and watch the blocks "
        "round inward. *Every parameter is a placeholder; this teaches the "
        "mechanism, not a rate.*",
        margin=(0, 10, 5, 10), sizing_mode="stretch_width"),
    pn.Row(run, reset, readout),
    pn.Row(angle, spacing, infiltration, sizing_mode="stretch_width",
           max_width=DESIGN_WIDTH),
    pn.Row(fig_left, fig_right, sizing_mode="stretch_width",
           max_width=DESIGN_WIDTH),
    sizing_mode="stretch_width", max_width=DESIGN_WIDTH,
).servable(title="corestone — fracture-controlled granite weathering")

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
# BLOCK, since that is what makes a corestone look round rather than stepped --
# a 1 m joint spacing is 20 cells across.
#
# NEARLY SQUARE, and it stays that way. A 2:1 section was tried, because the
# section's shape is the only thing that sets how tall the demo is on a page:
# two panels at a 900 px design width are 450 px each, less ~55 for the depth
# axis and ~65 for the colour bar, so the data area is 330 px wide and a square
# one is 330 px tall. Halving the depth halves that. But it also halves the
# rock: the blocks come out as letterbox slots rather than the roughly
# equidimensional joint-bounded cubes that a granite outcrop actually has, and
# a corestone that is twice as wide as it is tall is not the thing this demo
# exists to show. Reverted. The box is taller; the geometry is right.
LX = LZ = 3.00                  # the section, in metres. Fixed: it is the cell
                                # SIZE that varies, not the piece of rock.

#: Cell sizes offered. The section is 3 m square, so these are 60, 120 and 150
#: cells across. Finer is not simply better: what makes a corestone look round
#: rather than stepped is cells per BLOCK, and a 1 m joint spacing is already
#: 20 cells at 5 cm. What 2 cm buys is a sharper weathering rind; what it costs
#: is 6.3x the cells and about 18x the time, measured below.
#:
#:     cell    cells   per frame   200 kyr
#:     5 cm     3600      1.9 ms     0.34 s
#:     2.5 cm  14400     25.2 ms     5.21 s
#:     2 cm    22500     34.9 ms     7.72 s
#:
#: One frame is one step, so at 2 cm a frame already costs more than the 33 ms
#: animation budget here and several times that in a browser. It animates, just
#: slowly; Show result is the way to use it.
CELL_SIZES = {"5 cm": 0.05, "2.5 cm": 0.025, "2 cm": 0.02}

#: The index cap is not arbitrary: a high-index angle tiles only at a very
#: fine spacing (nine divisions of a three-metre section), so its spacing
#: slider would carry a single choice. Capping at 4 keeps the angles that offer
#: a real range of spacings.
MAX_INDEX = 4

#: Which snapped spacings to offer. Below 0.3 m a block is 6 cells across and
#: looks square however long it runs; above 3 m there is no block inside the
#: section at all.
SPACING_LOW, SPACING_HIGH = 0.3, 3.0

#: The longest run offered. 200 kyr dissolves the section at the default
#: settings, so there is nothing further to watch.
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
SLIDER_WIDTH = DESIGN_WIDTH // 5 - 14
#: Wider than it is tall, because a figure is not its data area: the depth
#: axis and its label take about 55 px on the left and the colour bar another
#: 60 on the right, while only the distance axis (~55 px) is below. Sized 1:1
#: the DATA would come out visibly taller than it is wide; 460 x 400 makes the
#: data square, which is what a 3.0 x 3.05 m section is.
FIG_W, FIG_H = 460, 400


#: The model's reference temperature, in the units the slider speaks. Both
#: temperature factors are exactly 1 here by construction.
T_REF_C = 285.0 - 273.15


def _cells(dx):
    """Cells across the section at this cell size."""
    return int(round(LX / dx))


def _angles(dx):
    """
    Rotations at which the joint pair tiles the periodic width exactly.

    A FUNCTION OF THE CELL COUNT, not of the section: the lattice indices have
    to divide it. 60 cells across has divisors 1..6 and offers seven angles;
    150 does not divide by 4, so 14.0 and 36.9 degrees are simply unavailable
    at 2 cm cells. That is why the cell size cannot be swapped as a constant --
    it changes what the other two sliders may offer.
    """
    return tiling_angles(_cells(dx), max_index=MAX_INDEX)


def _spacings(angle_deg, dx):
    """The snapped spacings available at this angle and cell size."""
    a, b = next((a, b) for ang, a, b in _angles(dx)
                if abs(ang - angle_deg) < 1e-6)
    return tiling_spacings(LX, a, b, SPACING_LOW, SPACING_HIGH, dx=dx)


# ---- widgets ----------------------------------------------------------------
# The angle and spacing sliders SNAP. A joint pair rotated by theta has
# along-x periods S/cos(theta) and S/sin(theta), and the section is periodic
# left-to-right -- it has no side walls, because a no-flow wall manufactures a
# domain-scale circulation. Both periods must divide the width for the joints
# to close on themselves, which needs tan(theta) = b/a for integers and
# quantises the spacing too. Off those values the joints fail to line up
# across the seam.
cell = pn.widgets.DiscreteSlider(
    name="Cell size", options=CELL_SIZES, value=0.05,
    sizing_mode="stretch_width", max_width=SLIDER_WIDTH)
angle = pn.widgets.DiscreteSlider(
    name="Joint orientation [°]",
    options={"%.1f°" % a: a for a, _, _ in _angles(0.05)}, value=0.0,
    sizing_mode="stretch_width", max_width=SLIDER_WIDTH)
spacing = pn.widgets.DiscreteSlider(
    name="Joint spacing [m]",
    options={"%.2f m" % s: s for s in _spacings(0.0, 0.05)}, value=1.0,
    sizing_mode="stretch_width", max_width=SLIDER_WIDTH)
infiltration = pn.widgets.FloatSlider(
    name="Infiltration [m/yr]", start=0.05, end=1.00, step=0.05,
    value=0.30, format="0.00",
    sizing_mode="stretch_width", max_width=SLIDER_WIDTH)
# Offered in degrees Celsius, because a reader thinks in a climate rather than
# in kelvin; the model is given kelvin. The default is the model's own
# reference temperature, 285 K, so the demo opens with both temperature
# factors at exactly 1 and the slider is the only thing that moves them.
temperature = pn.widgets.FloatSlider(
    name="Temperature [°C]", start=0.0, end=30.0, step=1.0,
    value=T_REF_C, format="0",
    sizing_mode="stretch_width", max_width=SLIDER_WIDTH)


def _build():
    """A fresh network and a fresh model at the current slider settings."""
    n = _cells(cell.value)
    net = FractureNetwork(n, n, cell.value, periodic_x=True).seed(
        sets=orthogonal_grid(spacing.value, rotation=angle.value),
        rng=np.random.default_rng(12345))
    m = Weathering(net)
    m.set_infiltration(infiltration.value / YEAR)
    m.set_temperature(temperature.value + 273.15)
    m.c_drift_max = C_DRIFT_MAX
    m.initialize()
    m.c = m.solve_solute(m.reaction_coefficient)
    return net, m


# `sim`, never `state`: panel exports pn.state, and shadowing it fails silently.
sim = {}


def step():
    """One weathering step per frame, stopping exactly at the chosen time."""
    m = sim["model"]
    target = stop_at.value * 1e3 * YEAR
    if m.t >= target - 1e-9 * YEAR:
        run.value = False                      # arrived; pause
        return
    m.update(dt_limit=target - m.t)
    _redraw()


def do_reset():
    """Rebuild from the sliders. Every slider here is structural."""
    net, m = _build()
    sim["net"], sim["model"] = net, m
    _joints()
    _redraw()


def show_result(event=None):
    """Jump straight to the state at the chosen time, without animating.

    Always from fresh rock, even when the model has not yet reached that time
    and could simply be advanced. Watching it evolve is one question; asking
    what the rock looks like at 50 kyr is another, and the answer to the second
    should not depend on what was pressed before it. Rebuilding makes this a
    function of the sliders alone, which is what lets two settings be compared:
    same time, same answer, every time.
    """
    run.value = False                          # stop animating, if it was
    do_reset()
    m = sim["model"]
    target = stop_at.value * 1e3 * YEAR
    while m.t < target - 1e-9 * YEAR:
        m.update(dt_limit=target - m.t)
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
affinity = ColumnDataSource(data={"image": [np.zeros((2, 2))]})
dissolved = ColumnDataSource(data={"image": [np.zeros((2, 2))]})
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
#: How far to run. A property of the RUN, not of the rock, which is why it
#: sits by the buttons rather than in the parameter row and why it does not
#: reset the model: a student sets 50 kyr, runs, changes the joint spacing --
#: which does reset -- and runs to 50 kyr again, comparing like with like.
#: Left where it is, it also continues: run to 50, raise it to 100, run again.
stop_at = pn.widgets.IntSlider(
    name="Run to [kyr]", start=10, end=int(END_KYR), step=10,
    value=int(END_KYR), sizing_mode="stretch_width", max_width=260)

run = animator(step)
reset = reset_button(do_reset, name="Fresh rock")
jump = pn.widgets.Button(name="Show result", button_type="default", width=120)
jump.on_click(show_result)

# Every slider rebuilds: the joint geometry is the initial condition, and the
# infiltration rate sets a flow field that is solved once and held. None of
# them is a forcing that can be turned while the rock evolves, so changing one
# starts the clock again rather than pretending otherwise.
for w in (angle, spacing, infiltration, temperature):
    w.param.watch(lambda event: do_reset(), "value")


def _resnap_spacing(a=None):
    """Which spacings tile depends on the angle AND the cell size."""
    opts = _spacings(angle.value if a is None else a, cell.value)
    nearest = min(opts, key=lambda s: abs(s - spacing.value))
    spacing.options = {"%.2f m" % s: s for s in opts}
    spacing.value = nearest


angle.param.watch(lambda event: _resnap_spacing(event.new), "value")


@pn.depends(cell.param.value, watch=True)
def _resnap_grid(dx):
    """
    A new cell size changes what the other two sliders may offer, so they are
    re-snapped before the model is rebuilt. Nearest surviving value each time
    rather than a default: someone comparing 26.6 degrees at two resolutions
    should not have the angle silently reset under them.
    """
    opts = _angles(dx)
    nearest = min((a for a, _, _ in opts), key=lambda a: abs(a - angle.value))
    angle.options = {"%.1f°" % a: a for a, _, _ in opts}
    angle.value = nearest
    _resnap_spacing()
    do_reset()


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
        "**Why a corestone survives** – press **▶** and watch the blocks "
        "round inward, or **Show result** to jump straight to the state at "
        "**Run to**. *Every parameter is a placeholder; this teaches the "
        "mechanism, not a rate.*",
        margin=(0, 10, 5, 10), sizing_mode="stretch_width"),
    pn.Row(run, reset, stop_at, jump, readout),
    pn.Row(angle, spacing, infiltration, temperature, cell,
           sizing_mode="stretch_width", max_width=DESIGN_WIDTH),
    pn.Row(fig_left, fig_right, sizing_mode="stretch_width",
           max_width=DESIGN_WIDTH),
    # Centred, not jammed left. The cap means the app can be narrower than the
    # frame -- whenever the embedding page has not scaled the frame to the
    # design width -- and left-aligned that reads as a broken layout with a
    # slab of empty space beside it rather than as a demo.
    sizing_mode="stretch_width", max_width=DESIGN_WIDTH, align="center",
).servable(title="corestone – fracture-controlled granite weathering")

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
import asyncio

import numpy as np
import panel as pn
from bokeh.models import (ColorBar, ColumnDataSource, FixedTicker,
                         LinearColorMapper, Range1d)
from bokeh.palettes import Blues256, Oranges256
from bokeh.plotting import figure

from artesian.live import animator, reset_button, responsive

from corestone import (FractureNetwork, Weathering, orthogonal_grid,
                       tiling_angles, tiling_spacings, YEAR)

pn.extension(loading_indicator=True)

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
#: is 6.3x the cells and about 18x the time.
#:
#:     cell     cells   median ms/step    vs 5 cm
#:     5 cm      3600      1.75  1.13        --
#:     2.5 cm   14400     17.51 10.84    10.0x  9.6x
#:     2 cm     22500     31.52 18.79    18.0x 16.6x
#:
#: Two independent median runs, disagreeing by about 10 % on the ratio. That
#: is the honest precision of this figure on a shared machine, so it is quoted
#: as "about eighteen times" and should NOT be re-chased: the scatter is
#: larger than any improvement worth making, and chasing it is what put a
#: wrong number on the exercise page once already.
#:
#: Measured as the MEDIAN over 250 steps, not as one timed run, and the step
#: count to a given time is identical at every cell size (800 steps and 86
#: flow solves to 200 kyr), so this ratio is the whole story. That method is
#: not fussiness: single timed runs of this on a loaded workstation returned
#: anywhere from 14x to 30x, and a 15x from one such run was published to the
#: exercise page before the scatter was noticed. At 5 cm the p90 step is
#: 13.81 ms against a 1.75 ms median -- contention, not arithmetic.
#:
#: Against the 33 ms frame budget, and the 3.0x this app runs slower under
#: Pyodide, only 5 cm keeps time: 2.5 cm and 2 cm come to roughly 53 and 95 ms
#: a frame in a browser. They animate, slowly; Show is the way to use them.
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

#: The longest moment Show will jump to. NOT a limit on Run, which keeps going
#: until it is paused -- and that distinction is what makes a finite cap here
#: safe, because NO cap finishes every setting. An earlier note claiming 1000
#: kyr did was simply wrong: temperature and infiltration multiply, so the time
#: to dissolve the section spans two orders of magnitude across the sliders.
#: Measured on the 3 m section at 5 cm, kyr to reach 50 / 90 / 99 % dissolved:
#:
#:     1.0 m, 0.30 m/yr, 12 C   default        1670   3713    5414
#:     1.0 m, 0.30 m/yr, 30 C   warm            455   1025    1490
#:     1.0 m, 0.30 m/yr,  0 C   cold           4563  10044   14678
#:     1.0 m, 0.05 m/yr, 12 C   dry            4933  10876   14340
#:     1.0 m, 0.05 m/yr,  0 C   both          14286 >20000  >20000
#:     3.0 m, 0.05 m/yr,  0 C   and coarse   >20000 >20000  >20000
#:
#: 200 kyr would cut the temperature comparison in half -- the default finishes
#: and 0 C does not -- which reads as the tool giving up rather than as a rate
#: difference, and comparing rates is what the slider is for. 15000 carries
#: the default well past 99 % at 5414 kyr and the cold case to 99 % at 14678,
#: which was the criterion that chose the original cap. The compound-slow
#: corners it does not reach; Run does, unbounded.
#:
#: The cap has moved 500 -> 4000 -> 15000 because the MODEL moved, twice:
#: deriving tau from the mineralogy, then correcting the matrix transport.
#: Each time it is the same criterion applied to new numbers, not a new
#: judgement.
#:
#: Re-measured after E_a and delta_H_r were sourced (oligoclase and quartz).
#: The 12 C row barely moved, because T_ref is 285 K and both factors are 1
#: there by construction; the ends spread, which is the point of sourcing them.
#:
#: It also bounds what Show costs, since the jump computes the whole run: at
#: 2 cm, 34.3 s here against 66.2 s for 1000 kyr. And it keeps the default's
#: interesting range in the first 40 % of a slider stepped in 10 kyr.
END_KYR = 15000.0

#: How much MODEL TIME one animation frame covers. The same for every setting,
#: and that is the whole point: it makes a second of watching worth a fixed
#: number of years, so a run that takes four times as long in the model takes
#: four times as long to watch.
#:
#: One frame used to be one drift-controlled step, which taught the opposite.
#: The controller holds the visible CHANGE per frame constant, so it hands a
#: slow-weathering run more years per frame -- measured at 5 cm, 1.84 kyr per
#: frame at 0 degrees C against 0.19 kyr at 30 -- and a cold section reached
#: 90 % dissolved in 149 frames where a warm one needed 333. Cold took 4.2x
#: longer in model time and less than half the real time. The animation was
#: teaching the reverse of the model. (Those five figures describe the old
#: behaviour and were measured before E_a and delta_H_r were sourced; they are
#: kept as the record of why this changed, not as current numbers.)
#:
#: Accuracy is unaffected: the frame sub-steps as c_drift_max demands, so this
#: sets the pace and the controller still sets the step.
#:
#: 1 kyr per frame. The pace rose to this when tau was derived from the
#: mineralogy rather than calibrated, which made the model sevenfold slower.
#: It went to 2 kyr briefly when correcting the matrix transport slowed it
#: another threefold, and came back: 2 kyr put the warm end over the frame
#: budget at 46 ms, and keeping the whole temperature range in budget is worth
#: more than halving the watch times, because the temperature comparison is
#: what the pace exists to protect.
#:
#: The reasoning that first chose 250 yr, kept because it is the argument that
#: sets the ceiling on any pace, and 1 kyr has to clear the same bar. The
#: frame budget is 33 ms. Local
#: cost per frame, and the same figure scaled by the 3.0x this app measured
#: slower under Pyodide in Chrome (Show to 50 kyr: 0.96 s against 0.32 s):
#:
#: At 1 kyr/frame, median over 250 frames, and scaled by the 3.0x this app
#: measured slower under Pyodide in Chrome:
#:
#:      5 cm,  12 C      2.97 ms  ->   8.9 ms     inside 33 ms
#:      5 cm,  30 C     11.78     ->  35.3        AT the budget
#:      2.5 cm, 12 C    55.16     -> 165.5        over
#:      2 cm,  12 C    137.65     -> 412.9        over
#:
#: The warm end now sits ON the budget rather than inside it -- two
#: measurements of the same thing gave 32 and 35 ms against 33 -- so 30 C will
#: drop the occasional frame. That is the price of converging the flow
#: tolerance, and it is affordable in a way it was not before artesian's
#: animator began yielding every frame: a dropped frame now stretches the run
#: instead of freezing the controls.
#:
#: 5 cm keeps time across the whole temperature range and the finer grids do
#: not, which has been true at every pace. At 2 kyr the warm end came to 46 ms
#: and did not: a longer frame gathers more sub-steps, and at 30 C the rock
#: changes fastest, so the fast end is where a longer frame breaks first. Superseded
#: measurements at 250 and 500 yr per frame:
#:
#:                        500 yr/frame        250 yr/frame
#:      0 C, default      1.8 ms   5 ms       2.1 ms   6 ms
#:     12 C, default      3.4      10         2.2      7
#:     30 C, default     11.9      36  OVER   6.4     19
#:
#: At 500 the warm end is already over budget in a browser, so it drops frames
#: and its wall-clock stretches -- compressing the very difference this is for.
#: At 250 the whole temperature range keeps up. Measured in Chrome on the
#: built page, against the 7.50 kyr per real second this asks for:
#:
#:      0 C   7.52 kyr/s     36 s to 90 % dissolved
#:     12 C   7.48           19 s
#:     30 C   7.04            9 s
#:
#: Uniform within 6 %, so what a reader feels is the model time, which is the
#: whole claim. With the sourced kinetics that buys, at 90 % dissolved:
#:
#:      0 C  10044 kyr / 30 kyr/s = 335 s of watching
#:     12 C   3713        / 30     = 124 s
#:     30 C   1025        / 30     =  34 s
#:
#: Those are long -- nearly six minutes to watch a cold section reach 90 %.
#: Run is pausable and Show exists for exactly this, and the alternative was
#: dropping frames at the warm end, which corrupts the comparison rather than
#: merely lengthening it.
#:
#: an 7.8x spread, and longer runs than before: the model slowed sevenfold
#: when tau was derived and the pace rose only fourfold, so watching costs
#: about 1.7x what it did. That is the price of the timescale being a
#: prediction rather than a calibration.
#:
#: What it costs: twice as many frames as 500, and the slowest corner on offer
#: (3 m joints, 0.05 m/yr, 0 degrees C) now passes 2000 kyr before reaching
#: 90 %, so more than four minutes of watching. That is the honest price of a proportional pace -- rock 30x
#: slower takes 30x longer to watch -- and Run is pausable.
#:
#: The limit, stated plainly: on a machine too slow to hold 30 fps the pace
#: stops being set by this number and starts being set by compute, and since a
#: warm run costs more per model year, the felt difference shrinks. Total
#: compute to 90 % is 1.19 s at 0 degrees C against 1.54 s at 30, so on a slow
#: enough laptop the two look equally long. Nothing here can fix that; it is
#: the frame budget, not the choice of pace.
YEARS_PER_FRAME = 1000.0

#: Tighter than the model's own default of 0.03, for two reasons that happen
#: to agree. One frame is one step, so the budget sets how long the animation
#: lasts: 0.03 is 58 frames, under two seconds at 30 fps, which is over before
#: a reader has focused on it. 0.01 is 180 frames, about six seconds. It also
#: cuts the error by roughly three, since the error is close to linear in this.
#: A step costs 2.9 ms here, so six seconds of animation is 0.5 s of arithmetic
#: and the frame budget is nowhere near threatened.
C_DRIFT_MAX = 0.01

#: How far the rock may weather before the head is re-solved. Looser than the
#: model's converged 0.01, because the feedback is expensive: re-solving the
#: head is what makes weathering open the rock and draw more water in, and it
#: costs about ten times as much as holding the flow fixed. Measured over
#: 200 kyr on the 3 m section at 5 cm, against the converged answer:
#:
#:     flow_tolerance   ms/frame   200 kyr   max|dM|
#:              0.01       11.9      3.45 s    0
#:              0.05        4.6      1.31 s    0.027
#:              0.10        2.6      0.73 s    0.059
#:
#: 0.05 keeps a frame inside the 33 ms animation budget with room for a
#: browser being several times slower, and 0.027 on a field in [0, 1] is far
#: inside the uncertainty on the conductivities themselves, which span an
#: order of magnitude in the measurements they come from.
# Converged from 0.05 on 2026-09-05. At 0.05 the head lagged the rock enough
# to move the answer by max|dM| = 0.381 against a 0.005 reference -- visibly
# wrong pictures. 0.02 gives 0.141 and 0.01 gives 0.045, but 0.01 puts the
# warm end well over the frame budget at 61 ms. 0.02 is the tightest tolerance
# that keeps the whole temperature range near it.
#
# It was 0.05 because it had to be: re-solving the head was half the run at
# 0.01. Warm-started CG cut the flow factorisations from 850 to 21 over
# 2000 kyr, which is what made converging it affordable. See
# Weathering._solve_head.
FLOW_TOLERANCE = 0.02

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
# WHICH REACTION. Two assignments live in this one app: the in-class activity
# runs feldspar dissolution, which is the textbook case -- an Arrhenius rate
# constant, a solubility ceiling, and water that stops working once it is full
# -- and the problem set runs biotite oxidation, where the solute is a
# REACTANT that runs out instead of a product that fills up.
#
# It is a radio group and not a slider because it is not a parameter. Nothing
# about it is continuous: it changes which equation the model is solving, and
# every label on the right-hand figure changes with it.
DRIVER_LABELS = {"Feldspar dissolution": "dissolution",
                 "Biotite oxidation": "oxidation"}
driver = pn.widgets.RadioButtonGroup(
    name="Reaction", options=list(DRIVER_LABELS), value="Feldspar dissolution",
    button_type="default", sizing_mode="stretch_width",
    max_width=2 * SLIDER_WIDTH)

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
    m.set_driver(DRIVER_LABELS[driver.value])
    m.set_infiltration(infiltration.value / YEAR)
    m.set_temperature(temperature.value + 273.15)
    m.c_drift_max = C_DRIFT_MAX
    m.flow_tolerance = FLOW_TOLERANCE
    m.dt_max = YEARS_PER_FRAME * YEAR
    m.initialize()
    m.c = m.solve_solute(m.reaction_coefficient)
    return net, m


# `sim`, never `state`: panel exports pn.state, and shadowing it fails silently.
sim = {}


def step():
    """Advance one frame -- a fixed span of MODEL TIME, for every setting.

    Not one step: one step is a fixed amount of visible change, which makes
    slow rock race and fast rock crawl. Sub-stepping inside the frame keeps
    the accuracy control in charge of the step while this stays in charge of
    the pace. See YEARS_PER_FRAME.

    No end. Run runs until it is paused, and the time selector belongs to
    :func:`show_result` alone -- a stopping point on the animation would make
    what you are looking at depend on which button you pressed to get there.
    """
    m = sim["model"]
    target = m.t + YEARS_PER_FRAME * YEAR
    while m.t < target - 1e-9 * YEAR:
        m.update(dt_limit=target - m.t)
    _redraw()


def _retitle():
    """Point the right-hand colour bar at whichever reaction is running.

    The field is the same array either way -- 1 - M -- and it does not mean
    the same thing, so the label is not decoration. Dissolving, it is mass
    that has left the rock; oxidising, it is iron that has rusted in place
    without leaving.
    """
    bar_right.title = EXTENT_LABEL[DRIVER_LABELS[driver.value]]


def do_reset():
    """Rebuild from the sliders. Every slider here is structural."""
    _retitle()
    net, m = _build()
    sim["net"], sim["model"] = net, m
    _joints()
    _redraw()


async def show_result(event=None):
    """Jump straight to the state at the chosen time, without animating.

    ASYNC, and the ``await`` below is the whole reason. Show rebuilds from
    fresh rock and integrates the entire interval, which at the far end of the
    slider is minutes: 15000 kyr costs about 56 s here and three times that in
    a browser, and 2 cm is slower again. Run at least redraws every frame, so
    a reader can see it working; Show used to return nothing at all until it
    was finished, which is indistinguishable from a hung page.
    
    So the figures are put in their loading state, and then control is handed
    back to the event loop for a moment. Without that yield Panel would send
    the loading state and the finished result in the same update and the
    reader would see neither. The model runs in a web worker, so the browser's
    main thread stays free and the indicator actually turns.

    Always from fresh rock, even when the model has not yet reached that time
    and could simply be advanced. Watching it evolve is one question; asking
    what the rock looks like at 50 kyr is another, and the answer to the second
    should not depend on what was pressed before it. Rebuilding makes this a
    function of the sliders alone, which is what lets two settings be compared:
    same time, same answer, every time.
    """
    run.value = False                          # stop animating, if it was
    figures.loading = True
    jump.disabled = True
    await asyncio.sleep(0.05)                  # let the indicator reach the page
    try:
        do_reset()
        m = sim["model"]
        target = at_time.value * 1e3 * YEAR
        while m.t < target - 1e-9 * YEAR:
            m.update(dt_limit=target - m.t)
        _redraw()
    finally:
        figures.loading = False
        jump.disabled = False


def _joints():
    """The joint traces, as segments in metres."""
    net = sim["net"]
    seg = np.array([[p0[0], p0[1], p1[0], p1[1]]
                    for p0, p1 in net.segments]) if net.segments else \
        np.zeros((0, 4))
    for src in (joints_left, joints_right):
        src.data = {"x0": seg[:, 0], "y0": seg[:, 1],
                    "x1": seg[:, 2], "y1": seg[:, 3]}


#: Colour-bar range for the water speed, as log10 of metres per year.
#:
#: Five decades, 0.0001 to 10 m/yr. Across every slider setting the speed runs
#: from 5e-5 m/yr in stagnant rock at the lowest infiltration to 20 m/yr in a
#: joint at the highest, and a range covering all of that would be seven
#: decades -- which washes the picture out, since at any ONE setting the
#: section only spans about four. So the ends clip a little at the extremes of
#: the infiltration slider, and only inside joints, which are already the
#: darkest thing on the map and lose nothing by saturating.
SPEED_LOG_LOW, SPEED_LOG_HIGH = -4.0, 1.0


def _speed_field(m):
    """Water speed in metres per year, logarithmically.

    In real units, not as a multiple of the infiltration rate: a reader can
    ask whether a metre a year is fast for groundwater, and cannot ask that of
    a ratio.

    Logarithmic because it spans four orders of magnitude at any one setting
    -- the matrix starts near a thousandth of what falls on the surface while
    the joints carry twenty times it -- and on a linear scale everything but
    the joints would be white, which is a picture of the joint network rather
    than of the flow.

    Fixed limits, not per-frame, so two settings can be compared: the same
    reason Show rebuilds from fresh rock.
    """
    return np.log10(np.maximum(m.darcy_speed, 1e-30) * YEAR)


def _redraw():
    m = sim["model"]
    speed.data = {"image": [_speed_field(m)]}
    dissolved.data = {"image": [m.dissolved_fraction]}
    fig_left.title.text = "How fast the water is moving"
    fig_right.title.text = "What is left of the rock"
    # Time and a mean, and nothing that needs a threshold. This used to read
    # "grus X %, corestone Y %", which was two claims the model cannot make.
    # Both rested on cut-offs in dissolved fraction that were chosen and never
    # justified, and neither word means a fraction dissolved: the sequence
    # fresh rock - saprock - saprolite - grus is defined by fabric and
    # mineralogy, and a corestone is a SHAPE, a rounded block surrounded by
    # weathered rock, not a cell that happens to be under a cut-off. Deep
    # intact bedrock counted as corestone.
    # The mean over cells IS the fraction of the section's soluble phase that
    # has gone, since every cell starts with the same amount, so this needs no
    # "mean" qualifier to be exact.
    readout.object = (
        "**%.0f kyr** &nbsp;·&nbsp; **%.0f %%** %s"
        % (m.t / YEAR / 1e3, 100 * m.dissolved_fraction.mean(),
           EXTENT_LABEL[m.driver]))


# ---- figures ----------------------------------------------------------------
# Depth increases DOWNWARD, so the y range runs from LZ at the bottom of the
# axis to 0 at the top. Row 0 of the array is the ground surface, and bokeh
# draws row 0 at the anchor and later rows at increasing y, which with a
# reversed range puts the surface at the top where it belongs -- no flip.
speed = ColumnDataSource(data={"image": [np.zeros((2, 2))]})
dissolved = ColumnDataSource(data={"image": [np.zeros((2, 2))]})
joints_left = ColumnDataSource(data={"x0": [], "y0": [], "x1": [], "y1": []})
joints_right = ColumnDataSource(data={"x0": [], "y0": [], "x1": [], "y1": []})


def _panel(source, joints, palette, label, labels, low=0.0, high=1.0,
           ticks=(0.0, 0.25, 0.5, 0.75, 1.0)):
    """
    One map of the section, with its colour bar and joint traces.

    ``labels`` relabels chosen ticks, because a bare number on a colour bar
    says nothing. For the speed it undoes the logarithm and gives metres per
    year; for the dissolved fraction it names the ends of the QUANTITY --
    "none" and "all" of the soluble phase. Not a rock type: calling the dark
    end "grus" was a petrological claim the model cannot support, since the
    sequence fresh rock, saprock, saprolite, grus is a matter of fabric and
    mineralogy and not of how much of one phase has gone.
    """
    fig = figure(width=FIG_W, height=FIG_H,
                 x_axis_label="Distance [m]", y_axis_label="Depth [m]",
                 x_range=Range1d(0, LX), y_range=Range1d(LZ, 0),
                 tools="", toolbar_location=None)
    mapper = LinearColorMapper(palette=palette, low=low, high=high)
    fig.image(image="image", x=0, y=0, dw=LX, dh=LZ, source=source,
              color_mapper=mapper)
    fig.segment("x0", "y0", "x1", "y1", source=joints,
                color="#2a2a2a", line_width=1, alpha=0.45)
    bar = ColorBar(color_mapper=mapper, width=8, title=label,
                   label_standoff=6, padding=4,
                   ticker=FixedTicker(ticks=list(ticks)),
                   major_label_overrides=dict(labels))
    fig.add_layout(bar, "right")
    # Half the design width each, since the two sit side by side. Without a
    # bound, responsive() defaults to 1200 PER FIGURE, so the row is entitled
    # to 2400 -- far wider than the app is laid out for.
    responsive(fig, aspect_ratio=float(FIG_W) / FIG_H,
               max_width=DESIGN_WIDTH // 2)
    # The BAR comes back too, because the right-hand one is retitled when the
    # reaction changes. Returned rather than dug out of the figure: bokeh
    # models do not take arbitrary attributes, and hiding it in ``tags``
    # worked but read as a trick.
    return fig, bar


# Palettes reversed so that 0 is pale and 1 is saturated: bokeh's 256-step
# ramps run dark to light.
# Water, so blue. The bar is labelled in metres per year, undoing the
# logarithm, so the numbers on it are speeds and not exponents.
fig_left, _ = _panel(speed, joints_left, Blues256[::-1],
                  "water speed [m/yr]",
                  {-4.0: "0.0001", -3.0: "0.001", -2.0: "0.01",
                   -1.0: "0.1", 0.0: "1", 1.0: "10"},
                  low=SPEED_LOG_LOW, high=SPEED_LOG_HIGH,
                  ticks=(-4.0, -3.0, -2.0, -1.0, 0.0, 1.0))
#: What the right-hand field MEANS, which is not the same in the two modes.
#: Dissolving, M is the soluble phase remaining and 1 - M is mass that has
#: left the rock. Oxidising, M is unoxidised Fe(II) and 1 - M is iron that has
#: rusted IN PLACE -- Goodfellow et al. (2016) put it as "major changes in
#: rock properties can occur with only minor element leaching". The same
#: picture, and not the same claim, so the label has to follow the driver.
EXTENT_LABEL = {"dissolution": "soluble phase dissolved",
                "oxidation": "biotite iron oxidised"}

fig_right, bar_right = _panel(
    dissolved, joints_right, Oranges256[::-1],
    EXTENT_LABEL["dissolution"], {0.0: "none", 1.0: "all"})

#: The two figures, named so :func:`show_result` can put them in their
#: loading state while it computes.
figures = pn.Row(fig_left, fig_right, sizing_mode="stretch_width",
                 max_width=DESIGN_WIDTH)

readout = pn.pane.Markdown("", sizing_mode="stretch_width",
                           margin=(0, 10, 0, 10))
#: Which moment to jump to. It belongs to Show, not to Run: Run animates from
#: fresh rock until there is nothing left to watch, and this asks a different
#: question -- what does the rock look like at 50 kyr? -- answered directly.
#: It does not reset the model, because it is not a property of the rock.
at_time = pn.widgets.IntSlider(
    name="View results at [kyr]", start=250, end=int(END_KYR), step=250,
    value=2000, sizing_mode="stretch_width", max_width=260)

run = animator(step)
reset = reset_button(do_reset, name="Fresh rock")
jump = pn.widgets.Button(name="Show", button_type="success", width=90)
jump.on_click(show_result)

# Every slider rebuilds: the joint geometry is the initial condition, and the
# infiltration rate sets a flow field that is solved once and held. None of
# them is a forcing that can be turned while the rock evolves, so changing one
# starts the clock again rather than pretending otherwise.
# ...and so does the driver, which is more than a rebuild: it changes which
# equation is being solved, so the cached operator and factorisation go with
# it. set_driver does that; this only has to start the clock again.
for w in (angle, spacing, infiltration, temperature, driver):
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
        "round inward, or set **View results at** and press **Show** to go "
        "straight there. Left is where the water goes; right is how far the "
        "reaction has got.",
        margin=(0, 10, 5, 10), sizing_mode="stretch_width"),
    # One row per way of driving the model: watch it happen, or ask what it
    # looks like at one moment. The state readout rides with the first, where
    # there is room for it, rather than taking a line of its own.
    pn.Row(run, reset, readout,
           sizing_mode="stretch_width", max_width=DESIGN_WIDTH),
    pn.Row(at_time, jump, align="end",
           sizing_mode="stretch_width", max_width=DESIGN_WIDTH),
    # The reaction gets its own line and sits above the sliders, because it
    # is not one of them: the sliders move parameters, this changes the
    # equation. It is also the line that switches between the two
    # assignments -- feldspar in class, biotite for the problem set.
    pn.Row(pn.pane.Markdown("**Reaction**", margin=(5, 8, 0, 10),
                            width=80),
           driver, pn.Spacer(sizing_mode="stretch_width"),
           sizing_mode="stretch_width", max_width=DESIGN_WIDTH),
    pn.Row(angle, spacing, infiltration, temperature, cell,
           sizing_mode="stretch_width", max_width=DESIGN_WIDTH),
    figures,
    # Centred, not jammed left. The cap means the app can be narrower than the
    # frame -- whenever the embedding page has not scaled the frame to the
    # design width -- and left-aligned that reads as a broken layout with a
    # slab of empty space beside it rather than as a demo.
    sizing_mode="stretch_width", max_width=DESIGN_WIDTH, align="center",
).servable(title="corestone – fracture-controlled granite weathering")

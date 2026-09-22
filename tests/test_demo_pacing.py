"""
The animation's pace, which is a teaching claim as much as a numerical one.

One frame must cover a fixed span of MODEL TIME, the same at every setting.
When a frame was instead one drift-controlled step, the controller held the
visible CHANGE per frame constant and so handed slow-weathering rock more
years per frame -- 1.84 kyr per frame at 0 degrees C against 0.19 kyr at 30,
measured at 5 cm. A cold section then reached 90 % dissolved in 149 frames
where a warm one needed 333: cold takes 4.2x longer in the model and less
than half the real time on the screen. The demo taught the reverse of the
model, which is worse than teaching nothing.
"""

import os
import sys

import numpy as np
import pytest

from corestone import YEAR

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir,
                                "interactive_demo"))
demo = pytest.importorskip("corestone_panel")


@pytest.mark.parametrize("tC", [0.0, 12.0, 30.0])
def test_a_frame_covers_the_same_model_time_at_every_temperature(tC):
    """Only the 30 C case bites: at 0 and 12 the drift control asks for a
    longer step than a frame and is clamped to the frame, so one step is one
    frame either way. That clamp is half the fix -- it is what stops cold rock
    sprinting -- and the 30 C case covers the other half, where a frame has to
    gather several steps to cover its span."""
    demo.temperature.value = tC
    demo.do_reset()
    m = demo.sim["model"]
    for _ in range(3):
        before = m.t
        demo.step()
        assert m.t - before == pytest.approx(demo.YEARS_PER_FRAME * YEAR,
                                             rel=1e-9)


def test_the_frame_respects_the_drift_budget():
    """
    The frame sets the pace; the drift control still sets the step, and this
    is the property that says so.

    It used to assert that a hot, wet frame takes MORE THAN ONE sub-step,
    which was true of the dissolution driver -- seven at 30 C and 1.00 m/yr --
    and is false of oxidation, which takes one step per frame at every setting
    the demo offers. That is not the accuracy control being switched off; it
    is a 1 kyr frame being short enough that one step stays inside the budget,
    and the model saying so. Measured, sub-steps in a 1 kyr frame at
    omega_drift_max = 0.01:

        driver        30 C 1.00   30 C 0.30   0 C 0.30
        dissolution       7           2          1
        oxidation         1           1          1

    So the invariant is asserted directly instead of through a proxy that only
    held for one of the two reactions: however many steps a frame takes, the
    drift each one produced is inside the budget.
    """
    demo.temperature.value = 30.0
    demo.rainfall.value = 1.00
    demo.do_reset()
    m = demo.sim["model"]
    drifts = []
    real = m.update
    def watched(*a, **k):
        out = real(*a, **k)
        drifts.append(m._drift)
        return out
    m.update = watched
    demo.step()
    demo.rainfall.value = 0.30                # leave the sliders as found
    assert drifts, "the frame took no step at all"
    assert all(d is not None and d <= m.omega_drift_max * 1.001 for d in drifts), \
        drifts


def test_the_reaction_control_switches_the_model_and_relabels_the_figure():
    """
    The demo carries two assignments: feldspar dissolution for the in-class
    activity, biotite oxidation for the problem set. The control that switches
    them has to change BOTH the equation and the label, because the
    right-hand field is the same array either way -- 1 - M -- and it does not
    mean the same thing.

    Dissolving, it is mass that has left the rock. Oxidising, it is iron that
    has rusted IN PLACE without leaving; Goodfellow et al. (2016) put that as
    "major changes in rock properties can occur with only minor element
    leaching". A picture relabelled wrongly would teach the second as the
    first.
    """
    for label, expected in demo.DRIVER_LABELS.items():
        demo.driver.value = label
        assert demo.sim["model"].driver == expected, label
        assert demo.bar_right.title == demo.EXTENT_LABEL[expected], label
    demo.driver.value = "Feldspar dissolution"     # leave it as found
    assert demo.sim["model"].driver == "dissolution"


def test_the_demo_opens_on_the_in_class_activity():
    """Feldspar, because that is the reaction the exercise teaches first.
    Pinned so that it cannot drift; see the library's own default test."""
    assert demo.driver.value == "Feldspar dissolution"
    assert demo.DRIVER_LABELS[demo.driver.value] == "dissolution"


def test_the_spacing_slider_reaches_no_joints_at_all():
    """
    The endpoint that shows what the joints were for.

    With no joints the rock can take only what its own matrix passes, so the
    surface ponds and most of the rain runs off. That setting is only
    meaningful because the surface can refuse water: without the ponding
    boundary the model forces the full rainfall through intact granite, which
    needs a hydraulic gradient of 23 and reports 67 m of head at the land
    surface.
    """
    # The unfractured case is INFINITE spacing, and sits at the far left of a
    # slider that runs descending, so the sequence stays monotonic.
    opts = list(demo.spacing.options.items())
    assert opts[0][1] == demo.NO_JOINTS, opts[0]
    assert opts[0][1] == float("inf")
    assert [v for _, v in opts[1:]] == sorted(
        [v for _, v in opts[1:]], reverse=True), opts

    demo.spacing.value = demo.NO_JOINTS
    m = demo.sim["model"]
    assert not m.network.link_v.any() and not m.network.link_h.any()
    assert m.ponded.all()
    assert m.infiltration < 0.1 * m.rainfall
    assert m.infiltration == pytest.approx(m.k_matrix_at_T, rel=0.05)

    # ...and any joints at all take every drop, because one 100 um joint
    # carries about twenty-three times the rain on a 3 m section.
    demo.spacing.value = 1.0
    m = demo.sim["model"]
    assert not m.ponded.any()
    assert m.infiltration == pytest.approx(m.rainfall, rel=1e-9)


def test_the_slider_is_named_for_what_arrives_not_what_enters():
    """A teaching page should not call the prescribed rate 'infiltration' when
    the model can now deliver less than it."""
    assert "Rainfall" in demo.rainfall.name
    assert not hasattr(demo, "infiltration")


def test_the_unjointed_case_is_solved_as_one_column():
    """
    With no joints the rock and the rain are both uniform across the section,
    so the solution is too. Solving one column imposes a symmetry the problem
    genuinely has, and the figures repeat it across the width.

    THE ALTERNATIVE IS WORSE, which is why this exists. Solved as a grid, the
    reactive-infiltration instability amplifies floating-point rounding error
    -- the only asymmetry a uniform section has -- and the model has nothing to
    select a wavelength with, diffusion through intact rock being ~7e-14 m2/s.
    Measured, the fingers are 2.1 cells wide at 5 cm cells and 13.3 at 2.5 cm:
    a pattern that follows the mesh rather than the rock.
    """
    demo.spacing.value = demo.NO_JOINTS
    m = demo.sim["model"]
    assert m.nx == 1, "no joints should be solved as a single column"
    assert not m.network.periodic_x, "one column has no seam to wrap"

    demo._redraw()
    img = demo.dissolved.data["image"][0]
    assert img.shape == (m.nz, demo._cells(demo.cell.value))
    assert np.ptp(img, axis=1).max() == 0.0, "every column is that column"

    # ...and a jointed run is still solved as a grid.
    demo.spacing.value = 1.0
    m = demo.sim["model"]
    assert m.nx == demo._cells(demo.cell.value)
    assert m.network.periodic_x

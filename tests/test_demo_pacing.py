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
    c_drift_max = 0.01:

        driver        30 C 1.00   30 C 0.30   0 C 0.30
        dissolution       7           2          1
        oxidation         1           1          1

    So the invariant is asserted directly instead of through a proxy that only
    held for one of the two reactions: however many steps a frame takes, the
    drift each one produced is inside the budget.
    """
    demo.temperature.value = 30.0
    demo.infiltration.value = 1.00
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
    demo.infiltration.value = 0.30                # leave the sliders as found
    assert drifts, "the frame took no step at all"
    assert all(d is not None and d <= m.c_drift_max * 1.001 for d in drifts), \
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

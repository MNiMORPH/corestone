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


def test_the_frame_still_sub_steps_for_accuracy():
    """The frame sets the pace; the drift control still sets the step. At the
    warm end a frame needs several sub-steps, and taking it in one would be
    the accuracy control quietly switched off."""
    demo.temperature.value = 30.0
    demo.infiltration.value = 1.00
    demo.do_reset()
    m = demo.sim["model"]
    calls = []
    real = m.update
    m.update = lambda *a, **k: (calls.append(1), real(*a, **k))[1]
    demo.step()
    demo.infiltration.value = 0.30                # leave the sliders as found
    assert len(calls) > 1, "a hot, wet frame should not be a single step"

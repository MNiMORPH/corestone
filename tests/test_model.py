"""
The skeleton runs and its lifecycle behaves.

REPLACE these with tests of the actual physics. Keep the shape: one test per
claim, and a name that states the claim.
"""

import numpy as np

from corestone import Model


def _spike(n=11):
    x = np.arange(0., float(n), 1.)
    z = np.zeros(n)
    z[n // 2] = 1.
    return x, z


def test_initialize_records_the_initial_state():
    m = Model()
    x, z = _spike()
    m.initialize(x, z)
    assert m.t == 0.
    assert len(m.z_out) == 1
    assert np.array_equal(m.z_out[0], z)


def test_run_advances_time_and_records_every_step():
    m = Model()
    m.set_k(1.)
    m.set_dt(0.1)
    m.initialize(*_spike())
    m.run(5)
    m.finalize()
    assert m.t == 5 * 0.1
    assert m.z_out.shape == (6, 11)          # initial state plus five steps
    assert np.allclose(m.t_out, np.arange(6) * 0.1)


def test_diffusion_lowers_the_peak_and_holds_the_ends():
    """The fixed-value ends do not move; the interior spike spreads."""
    m = Model()
    m.set_k(1.)
    m.set_dt(0.1)
    x, z = _spike()
    m.initialize(x, z)
    m.run(5)
    m.finalize()
    assert m.z_out[-1][5] < m.z_out[0][5]     # peak decayed
    assert m.z_out[-1][4] > 0.                # spread to its neighbours
    assert m.z_out[-1][0] == 0.               # ends held
    assert m.z_out[-1][-1] == 0.

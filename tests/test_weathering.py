"""
The affinity term does what the design says, and the flow loses no water.
"""

import numpy as np
import pytest

from corestone import FractureNetwork, Weathering, orthogonal_grid, YEAR


def _model(nz=75, nx=100, dx=0.20, spacing=1.5):
    net = FractureNetwork(nz, nx, dx).seed(sets=orthogonal_grid(spacing),
                                           rng=np.random.default_rng(12345))
    return Weathering(net)


def test_the_rock_starts_fresh_and_the_water_starts_clean():
    m = _model().initialize()
    assert np.all(m.dissolved_fraction == 0.0)
    assert np.all(m.c == 0.0)
    assert m.t == 0.0


def test_the_flow_routing_conserves_water():
    """
    Every cell hands on all of its water, and the edge columns have no
    off-grid receiver, so each row carries the whole infiltration.
    """
    m = _model().initialize()
    top = m.q[0, :].sum()
    for iz in range(m.nz):
        assert m.q[iz, :].sum() == pytest.approx(top, rel=1e-12)


def test_water_enters_fresh_and_saturates_with_depth():
    m = _model().run(years=50e3)
    assert m.c[0, :].max() == 0.0                    # rain arrives undersaturated
    assert np.median(m.c[-1, :]) > np.median(m.c[5, :])
    assert m.c.max() <= 1.0 + 1e-12
    assert np.allclose(m.affinity, 1.0 - m.c)


def test_raising_the_temperature_shortens_the_equilibration_length():
    """
    The counterintuitive core of the model: a hotter system does LESS
    weathering per metre of flow path, because the water saturates sooner.
    """
    m = _model()
    m.set_temperature(275.0)
    cold = m.equilibration_length
    m.set_temperature(305.0)
    hot = m.equilibration_length
    assert hot < cold
    assert m.equilibration_length > 0.0


def test_hotter_leaves_more_corestone_not_less():
    """
    Ten times the rate constant, and MORE rock survives -- the weathering
    concentrates at the joints instead of spreading. Warm is not weathered.
    """
    cold, hot = _model(), _model()
    cold.set_temperature(275.0)
    hot.set_temperature(305.0)
    cold.run(years=100e3)
    hot.run(years=100e3)
    assert hot.equilibration_length < cold.equilibration_length
    assert hot.is_corestone.mean() > cold.is_corestone.mean()


def test_grus_forms_at_the_joints_and_corestones_away_from_them():
    """The claim the whole model rests on."""
    m = _model().run(years=100e3)
    d = m.network.distance_to_fracture()
    assert m.is_grus.any() and m.is_corestone.any()
    assert np.median(d[m.is_grus]) < np.median(d[m.is_corestone])


def test_weathering_only_ever_advances():
    m = _model().initialize()
    previous = m.dissolved_fraction.sum()
    for _ in range(5):
        m.update()
        now = m.dissolved_fraction.sum()
        assert now >= previous
        previous = now
    assert m.t > 0.0


def test_no_rain_means_no_weathering():
    m = _model()
    m.set_infiltration(0.0)
    m.run(years=100e3)
    assert np.all(m.dissolved_fraction == 0.0)

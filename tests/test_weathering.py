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


def test_the_flow_solution_conserves_water():
    """
    Steady flow with no sources below the surface: whatever infiltrates must
    cross every horizontal plane and leave at the base.
    """
    m = _model().initialize()
    inflow = m.infiltration * m.dx * m.nx
    for iz in range(m.nz - 1):
        assert m.q_v[iz, :].sum() == pytest.approx(inflow, rel=1e-9)


def test_flow_is_downward_everywhere_so_rows_can_be_swept():
    """
    A flow driven by the gradient of a potential cannot circulate. The solute
    step sweeps rows in order and relies on this; it is checked, not assumed.
    """
    m = _model().initialize()
    assert (m.q_v < 0.0).sum() == 0


def test_the_horizontal_joints_carry_water():
    """
    The reason for solving a flow equation rather than routing water downhill.
    A gravity cascade can enter a subhorizontal joint but never travel along
    one, which left the whole horizontal set inert.

    The contrast grows with resolution -- 9x at dx = 0.20 m and 34x at
    0.05 m -- because a coarse grid marks a larger fraction of links as
    fractured and dilutes it. The threshold here is set for the coarse case.
    """
    m = _model().initialize()
    net = m.network
    on_joint = np.abs(m.q_h[net.link_h]).mean()
    intact = np.abs(m.q_h[~net.link_h]).mean()
    assert on_joint > 5.0 * intact


def test_water_enters_fresh_and_saturates_with_depth():
    """
    ``c`` is the concentration leaving a cell, not entering it, so the top row
    is small rather than exactly zero: rain arrives fresh and picks up a little
    on its way through.

    Concentration is NOT monotonic with depth. The top rows have given up most
    of their soluble phase, so L_eq there is long and the water leaves them
    nearly fresh; the matrix saturates within a metre or two below that. The
    bottom row runs lower again, which is not yet explained -- see the open
    defect in design/06.
    """
    m = _model().run(years=50e3)
    assert m.c[0, :].max() < 0.05                    # rain arrives fresh
    assert np.median(m.c[5, :]) > 0.9                # matrix saturates quickly
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


def test_temperature_barely_matters_because_transport_limits_it():
    """
    Shortening the equilibration length 28-fold changes the total weathering by
    under a tenth. Once L_eq is far shorter than the joint spacing, the water
    saturates within a cell of the joint whatever the rate constant, so the
    rock dissolved is set by how much water arrives and how much it can carry
    -- not by how fast the rock dissolves. That is the transport-limited
    (high-Damkohler) regime.

    This replaces an earlier test asserting that hotter leaves MORE corestone.
    That held under the gravity cascade, where a shorter L_eq visibly
    concentrated the weathering at the joints. Under Darcy flow the water is
    already concentrated at the joints, so the effect is small and not
    monotonic -- it reverses between 275 and 285 K. The earlier result was
    partly an artifact of the flow model, and the test failed honestly when the
    flow was corrected.
    """
    cold, hot = _model(), _model()
    cold.set_temperature(275.0)
    hot.set_temperature(315.0)
    cold.run(years=100e3)
    hot.run(years=100e3)

    assert cold.equilibration_length / hot.equilibration_length > 20.0
    change = abs(hot.dissolved_fraction.mean() - cold.dissolved_fraction.mean())
    assert change / cold.dissolved_fraction.mean() < 0.15


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
    """
    No infiltration, no flux, no dissolution. The residual is the head solver's
    roundoff -- the exact solution is a uniform head and identically zero flux.
    """
    m = _model()
    m.set_infiltration(0.0)
    m.run(years=100e3)
    assert np.abs(m.q_v).max() < 1e-15
    assert m.dissolved_fraction.max() < 1e-8


def _periodic_model(nz=None, nx=None, dx=0.10, spacing=1.5):
    # dx must divide the spacing a whole number of times for a periodic
    # tiling; 1.5 / 0.20 = 7.5 does not, and periodic_grid_shape now says so.
    from corestone import periodic_grid_shape
    if nz is None:
        nz, nx = periodic_grid_shape(20.0, 15.0, dx, spacing)
    net = FractureNetwork(nz, nx, dx, periodic_x=True).seed(
        sets=orthogonal_grid(spacing), rng=np.random.default_rng(12345))
    return Weathering(net)


def test_periodic_walls_make_every_block_weather_alike():
    """
    A no-flow wall is not neutral. It forces the lateral flow to vanish there,
    and with subhorizontal joints that manufactures a domain-scale circulation:
    the centre block weathered a third as much as the blocks two in from the
    walls, and the effect grew with the width of the section rather than
    staying near the edges. Wrapping the section onto itself removes the walls,
    and every block then behaves identically.
    """
    m = _periodic_model().run(years=100e3)
    net = m.network
    j = np.nonzero(net.link_v[m.nz // 2, :])[0]
    X = m.dissolved_fraction
    blocks = np.array([X[:, a + 1:b].mean() for a, b in zip(j[:-1], j[1:])])
    assert len(blocks) >= 5
    spread = (blocks.max() - blocks.min()) / blocks.mean()
    assert spread < 0.01


def test_the_periodic_network_tiles_across_the_seam():
    """The wrap gap is one spacing like every other, not a doubled joint."""
    m = _periodic_model()
    net = m.network
    j = np.nonzero(net.link_v[net.nz // 2, :])[0]
    gaps = np.diff(j)
    assert len(set(gaps)) == 1                       # uniform inside
    assert net.nx - j[-1] + j[0] == gaps[0]          # and across the seam
    # The subhorizontal joints reach both walls, so the seam link conducts.
    assert net.link_wrap.sum() > 0


def test_periodic_flow_still_conserves_water():
    m = _periodic_model().initialize()
    inflow = m.infiltration * m.dx * m.nx
    for iz in range(m.nz - 1):
        assert m.q_v[iz, :].sum() == pytest.approx(inflow, rel=1e-8)

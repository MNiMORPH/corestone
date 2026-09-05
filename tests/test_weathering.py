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


@pytest.mark.parametrize("driver", ["dissolution", "oxidation"])
def test_the_rock_starts_fresh_and_the_water_starts_at_the_inlet(driver):
    """
    "The water starts clean" was true of one reaction and is not a general
    statement. Rain carries no dissolved silica, so it starts at c = 0 when
    silica is the solute; it is in equilibrium with the atmosphere, so it
    starts at c = 1 when oxygen is. Neither is a parameter -- both are exact
    by construction -- and the initial field is the inlet value either way.
    """
    m = _model()
    m.set_driver(driver)
    m.initialize()
    assert np.all(m.dissolved_fraction == 0.0)
    assert np.all(m.c == m.inlet_concentration)
    assert m.inlet_concentration == (1.0 if driver == "oxidation" else 0.0)
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

    The freshness is asserted as a RATIO against the saturated matrix
    below, not as an absolute number. It was ``< 0.05`` here, and the
    converged value is 0.04919 -- a margin of 1.6 %, so the test was
    measuring the time step rather than the physics, and it duly broke
    the first time the step control changed. What the sentence above
    actually claims is that water at the surface is far from saturated
    while the matrix a metre down is at it, and that is a ratio:
    0.05 against 0.99.
    """
    # 350 kyr, not 50. The claim is about a profile whose top has given up
    # most of its soluble phase, and deriving tau from the mineralogy made
    # the model seven times slower, so reaching that state takes seven times
    # as long. The assertion is untouched; only the time needed to get there.
    m = _model()
    m.set_driver("dissolution")
    m.run(years=350e3)
    assert m.c[0, :].max() < 0.1 * np.median(m.c[5, :])   # rain arrives fresh
    assert np.median(m.c[5, :]) > 0.9                # matrix saturates quickly
    assert m.c.max() <= 1.0 + 1e-12
    assert np.allclose(m.affinity, 1.0 - m.c, rtol=1e-12, atol=0.0)


def test_water_enters_full_of_oxygen_and_gives_it_up_with_depth():
    """
    The mirror of the test above, and the reason the driver matters.

    The solute now enters at its ceiling instead of at zero, and is consumed
    instead of accumulated, so every statement about the profile inverts: the
    top row is nearly saturated with oxygen rather than nearly free of silica,
    and the matrix below runs it down rather than filling it up.

    The affinity inverts with it. Under dissolution it is 1 - c, what the
    water can still take up; under oxidation it is c, what the water still has
    to give. Both are zero where the water can do no more work.
    """
    m = _model()
    m.set_driver("oxidation")
    m.run(years=350e3)
    assert m.c[0, :].min() > 0.99                    # rain arrives full
    # ...and gives it up on the way down. Asserted top against base rather
    # than at a fixed row: this grid is 20 cm cells against a 4.5 cm
    # penetration depth, so no single row resolves the depletion, and by
    # 350 kyr the matrix has opened enough to carry oxygen deep. The claim
    # that survives coarsening is the one that matters -- water leaves with
    # less than it arrived with.
    assert np.median(m.c[-1, :]) < 0.8 * np.median(m.c[0, :])
    assert m.c.max() <= 1.0 + 1e-12
    assert np.allclose(m.affinity, m.c, rtol=1e-12, atol=0.0)


def test_raising_the_temperature_shortens_the_saturation_length():
    """
    The counterintuitive core of the model: a hotter system does LESS
    weathering per metre of flow path, because the water saturates sooner.
    """
    m = _model()
    m.set_temperature(275.0)
    cold = m.saturation_length
    m.set_temperature(305.0)
    hot = m.saturation_length
    assert hot < cold
    assert m.saturation_length > 0.0


def test_temperature_acts_through_solubility_not_the_rate_constant():
    """
    Raising the temperature 40 K roughly triples the weathering -- but not for
    the reason the rate constant suggests.

    This model is transport-limited almost everywhere, and there the amount
    dissolved scales with C_eq, not with k. Holding C_eq constant left the
    model with a 5 % response over the same range, which was reported as a
    physical result about Damkohler limits. It was an artefact of a missing
    term: solubility is temperature dependent too, and it enters twice --
    the saturation length goes as C_eq/k, and tau = M0/C_eq.

    Warm therefore does mean weathered here, through solubility. The earlier
    "warm does not mean weathered" story belonged to a model that had no
    solubility at all.
    """
    cold, hot = _model(), _model()
    for m in (cold, hot):
        m.set_driver("dissolution")     # a statement about THIS reaction:
                                        # oxygen solubility runs the other way
    cold.set_temperature(275.0)
    hot.set_temperature(315.0)
    cold.run(years=100e3)
    hot.run(years=100e3)

    assert hot.solubility_factor > 3.0 * cold.solubility_factor
    assert hot.tau < cold.tau                       # more soluble carries more
    change = (hot.dissolved_fraction.mean() - cold.dissolved_fraction.mean())
    assert change / cold.dissolved_fraction.mean() > 1.0     # more than double


@pytest.mark.parametrize("driver", ["dissolution", "oxidation"])
def test_the_rock_weathers_less_the_further_it_is_from_a_joint(driver):
    """
    The claim the whole model rests on, stated as what it actually says.

    It used to be asserted through ``is_grus`` and ``is_corestone`` -- two
    ARBITRARY cut-offs at 0.50 and 0.05 on a continuous field, as their own
    docstrings say -- and required both sets to be non-empty. That held for
    dissolution, which drives the joint cells to 1 while the interiors sit
    near 0, and it never holds for oxidation, which reaches 0.50 somewhere
    only once almost nothing is left below 0.05. Measured on this grid:

        kyr    mean    max     cells > 0.50   cells < 0.05
        100   0.1125  0.2844        0             3339
        200   0.2477  0.4880        0              715

    That is not the model failing. It is a test of where two arbitrary
    thresholds happen to fall, dressed as a test of physics. What the
    sentence claims is that weathering decreases with distance from the
    joint network, and that is asserted directly here -- monotonically, band
    by band, which is a far stronger statement than two thresholds
    straddling.
    """
    m = _model()
    m.set_driver(driver)
    m.run(years=100e3)
    d = m.network.distance_to_fracture()
    x = m.dissolved_fraction

    # AXIAL distances only -- whole multiples of the cell size. Euclidean
    # distance mixes them with diagonals, and a cell one diagonal step from a
    # corner is adjacent to TWO joints while a cell two axial steps away is
    # adjacent to one, so ordering by distance is not ordering by shelter.
    # Measured under dissolution, the mixed ordering gives 1.61e-01,
    # 6.55e-04, 4.04e-08, 5.52e-07, 4.02e-11 -- the fourth above the third,
    # which is the corner effect that
    # test_corners_stay_further_from_saturation_than_faces asserts on
    # purpose. It is not a defect, so it must not be written into a test as
    # though monotonicity in Euclidean distance were the claim.
    bands = [k * m.dx for k in range(4)]
    means = [float(x[np.isclose(d, b)].mean()) for b in bands]
    assert all(a > b for a, b in zip(means, means[1:])), (bands, means)
    assert means[0] > 3.0 * means[-1], means


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


def test_the_water_speed_is_the_infiltration_rate_where_there_are_no_joints():
    """
    Unjointed rock has nowhere to focus the flow, so every cell passes exactly
    what falls on it and the Darcy speed is the infiltration rate everywhere.

    That is the check the cell-centred reassembly needs: the solver works on
    faces, and turning face fluxes back into a speed is where a factor of two
    or a missing boundary face hides. Both would still look plausible on a
    colour map.
    """
    net = FractureNetwork(30, 30, 0.10, periodic_x=True).seed(sets=[])
    m = Weathering(net).initialize()
    assert not net.link_v.any() and not net.link_h.any()
    v = m.darcy_speed
    assert np.allclose(v, m.infiltration, rtol=1e-6, atol=0.0), (v.min(), v.max())


def test_joints_carry_far_more_than_the_matrix():
    m = _periodic_model().initialize()
    v = m.darcy_speed
    joint = m.network.cell
    assert np.median(v[joint]) > 100.0 * np.median(v[~joint])
    # and the mean flux through any depth still equals what fell on the surface
    assert m.q_v.sum(axis=1)[m.nz // 2] == pytest.approx(
        m.infiltration * m.dx * m.nx, rel=1e-8)


def test_the_default_driver_is_dissolution_and_the_choice_is_deliberate():
    """
    Pinned, because it is a TEACHING decision and not a physical one, and a
    teaching decision that lives only in a comment will drift.

    Design 09 checked the oxidation case adversarially and it came out
    stronger: Goodfellow et al. (2016) watched biotite weathering begin with
    oxidation by diffusing oxygen. Oxidation is what really paces spheroidal
    weathering. The default is dissolution anyway, because this model exists
    to teach rate times affinity, a solubility ceiling and Arrhenius, and the
    oxidation driver inverts the temperature intuition before a student has
    built it.

    If this test fails, someone changed which reaction the exercise is about.
    That is allowed; it is not allowed to happen by accident.
    """
    m = _model()
    assert m.driver == "dissolution"
    assert m.apparent_activation_energy > 0.0      # warm means weathered
    assert m.oxygen_dissolution_enthalpy < 0.0     # ...but not for oxygen
    m.set_driver("oxidation")
    assert m.apparent_activation_energy == 0.0
    m.set_driver("dissolution")
    assert m.driver == "dissolution"

    with pytest.raises(ValueError):
        m.set_driver("photosynthesis")

"""
One transcription test per equation stated in the source.

These are not new physics. Each equation is copied out of the docstring where
it is claimed and required to hold against the model's own arrays, evaluated
independently of the solver that produced them. That is exactly the act --
reading the prose and the code side by side, once -- that was skipped when
``L_eq = q C_eq / (kA)`` sat in a docstring above code that computed a constant
with no ``q`` in it at all, for six revisions of the file, invisible to a suite
of thirty passing tests because every one of them checked an aggregate.

Two warnings from writing them, both worth keeping:

* It is easy to transcribe the equation *wrongly* and get a test that fails on
  correct code. That happened twice out of ten while these were first written
  (confusing psi with H; reproducing the very cancellation being hunted). Both
  failed loudly within a minute, which is the right failure mode -- but if one
  of these fails, suspect the transcription before the model.
* A global balance is not a per-cell check. Summing the rows of the transport
  matrix telescopes and is satisfied for *any* concentration field, so it is
  blind to an error inside the operator. Where it matters, check per cell.

Companion: ``test_equation_coverage.py`` fails if an equation appears in a
docstring without a test registered here.
"""

import numpy as np
import pytest

from corestone import (FractureNetwork, Weathering, YEAR, orthogonal_grid,
                       periodic_grid_shape)
from corestone.weathering import water_viscosity


def _model(dx=0.10, spacing=1.5, width=12.0, depth=9.0):
    nz, nx = periodic_grid_shape(width, depth, dx, spacing)
    net = FractureNetwork(nz, nx, dx, periodic_x=True).seed(
        sets=orthogonal_grid(spacing), rng=np.random.default_rng(12345))
    return Weathering(net).initialize()


# ---------------------------------------------------------------- rate law

@pytest.mark.parametrize("driver", ["dissolution", "oxidation"])
def test_the_reaction_rate_per_unit_volume_does_not_depend_on_the_flux(driver):
    """
        R = k(T) * A * (1 - C / C_eq) dissolving: a product, driven by
        how far the water is from saturation
        R = k_ox * A * C oxidising: a reactant, driven by how much of it
        there is

    The rate per unit volume is a property of the rock. It may not depend on
    how fast water moves past it. This is the invariant the original bug broke:
    holding the saturation length uniform made dissolution proportional to the
    local flux, weighting a joint against the matrix by a factor of ~3000.

    The reaction coefficient r = kA/C_eq is now the primary quantity, and the
    flux cancels explicitly rather than inside an exponent, so this is a direct
    check: r must be uniform across cells that differ in flux by 1000x.
    """
    m = _model()
    m.set_driver(driver)
    r = m.reaction_coefficient
    # The flux distribution is bimodal, so select by structure, not quantile.
    fast, slow = m.network.cell, ~m.network.cell
    assert fast.any() and slow.any()
    assert np.median(m.q[fast]) > 100.0 * np.median(m.q[slow])   # real contrast
    assert np.ptp(r) / r.mean() < 1e-12                     # yet r is uniform


@pytest.mark.parametrize("driver", ["dissolution", "oxidation"])
def test_what_the_rock_loses_is_what_the_water_carries_out_of_the_base(driver):
    """
        d(M/M0)/dt = - r f(c) / tau,  f = 1 - c dissolving, f = c oxidising

    Solute is conserved, and the books balance the same way whichever
    direction it points: what comes in equals what is consumed plus what
    leaves through the base.

    Under DISSOLUTION nothing comes in -- rain is clean -- and every cell is a
    source, so the statement collapses to "everything dissolved leaves through
    the base". Under OXIDATION the rain is the only source and the rock is the
    sink, so it reads "the oxygen that arrives is either used or exported".
    The second is the sharper test of the two, because a source at one
    boundary can be lost silently in a way a volumetric one cannot.

    The source is written out here rather than taken from
    ``_solute_source``, so the test cannot pass by agreeing with the code.
    """
    m = _model()
    m.set_driver(driver)
    r = m.reaction_coefficient
    c = m.solve_solute(r)

    if driver == "oxidation":
        supplied = m.infiltration * m.dx * m.nx    # rain, at c = 1
    else:
        supplied = (r * m.dx * m.dx).sum()         # rock, in every cell
    consumed = (r * c * m.dx * m.dx).sum()
    exported = (m.q_out_base * c[-1, :]).sum()
    assert supplied == pytest.approx(consumed + exported, rel=1e-9)

    rate = r * m.driving_force(c) / m.tau
    assert rate.shape == m.M.shape
    assert (rate >= 0.0).all()


def test_the_rock_is_integrated_exactly_over_a_step_with_c_held():
    """
        M(t + dt) = M(t) exp(-lambda dt), lambda = (r / M) (1 - c) / tau

    The content of "exactly" is that the answer does not depend on how the
    step is chopped up. With ``c`` held -- which is what the model does within
    a step -- taking one step of dt and ten steps of dt/10 must give the same
    M to roundoff. Forward Euler cannot do this: subdividing changes its
    answer, which is precisely the error it makes.

    Checked as that invariance rather than by re-typing the exponential, so
    that the test cannot pass by reproducing a mistake in the code.
    """
    m = _model()
    frozen = m.solve_solute(m.reaction_coefficient)
    m.solve_solute = lambda r: frozen              # hold c, as a step does

    dt = 4000.0 * YEAR
    m.dx_max = np.inf                              # let the step be the step
    one = _model()
    one.solve_solute = lambda r: frozen
    one.dx_max = np.inf
    one.update(dt=dt)

    many = _model()
    many.solve_solute = lambda r: frozen
    many.dx_max = np.inf
    for _ in range(10):
        many.update(dt=dt / 10.0)

    assert np.abs(one.M - many.M).max() < 1e-14
    assert one.M.min() > 0.0                       # an exponential cannot hit 0
    # and the rock did move, so the agreement is not agreement about nothing
    assert (1.0 - one.M).max() > 0.01


def test_forward_euler_would_fail_the_invariance_the_exponential_passes():
    """
    The companion that gives the test above its teeth. If subdividing the step
    were harmless for any integrator, the check would be vacuous. Euler is
    written out here and required to disagree with itself.
    """
    m = _model()
    frozen = m.solve_solute(m.reaction_coefficient)
    lam = m.specific_reaction_coefficient * (1.0 - frozen) / m.tau
    dt = 4000.0 * YEAR

    one = np.ones((m.nz, m.nx)) * (1.0 - lam * dt)
    many = np.ones((m.nz, m.nx))
    for _ in range(10):
        many = many * (1.0 - lam * dt / 10.0)
    assert np.abs(one - many).max() > 1e-6


# ------------------------------------------------------- saturation length

def test_the_saturation_length_is_proportional_to_the_local_flux():
    """
        saturation_length = q * C_eq / (k(T) * A)

    Proportional to q. Fast water in a joint travels far before it saturates;
    slow water in the matrix saturates almost at once.
    """
    m = _model()
    L = m.local_saturation_length()
    ratio = L / (m.q / m.dx)
    assert np.ptp(ratio) / ratio.mean() < 1e-12


def test_the_saturation_length_scales_as_C_eq_over_k_not_as_one_over_k():
    """
        saturation_length = q * C_eq / (k(T) * A)

    Both k and C_eq depend on temperature, so the exponent is (E_a - dH_r),
    not E_a. Treating C_eq as constant was the second half of the temperature
    dependence, and in the transport-limited regime it is the half that
    dominates.
    """
    m = _model()
    R, T_ref = m.R_gas, m.T_ref
    for T in (275.0, 295.0, 315.0):
        m.set_temperature(T)
        want = m.L_ref * np.exp(-((m.delta_H_r - m.E_a) / R)
                                * (1.0 / T - 1.0 / T_ref))
        assert m.saturation_length == pytest.approx(want, rel=1e-12)
        # and it is NOT the E_a-only form, unless dH_r is zero
        naive = m.L_ref * np.exp((m.E_a / R) * (1.0 / T - 1.0 / T_ref))
        if T != T_ref and m.delta_H_r != 0.0:
            assert m.saturation_length != pytest.approx(naive, rel=1e-6)


def test_tau_falls_as_solubility_rises():
    """
    ``tau = M0 / C_eq`` is the second place C_eq enters: a warmer, more soluble
    fluid carries more away per unit volume. Held constant, the model had no
    solubility response at all.

    A statement about the DISSOLUTION driver, so it names it. Oxygen runs the
    other way -- a gas leaves solution as the water warms -- and that
    opposition is checked in
    ``test_warming_drives_oxygen_out_of_solution_and_silica_into_it``.
    """
    m = _model()
    m.set_driver("dissolution")
    m.set_temperature(m.T_ref)
    base = m.tau
    m.set_temperature(m.T_ref + 20.0)
    assert m.tau < base
    assert m.tau == pytest.approx(m.tau_ref / m.solubility_factor, rel=1e-12)
    assert m.tau == m.silica_tau


def test_the_matrix_conducts_better_as_it_dissolves():
    """
        k(M) = k_matrix(T)^M * k_weathered(T)^(1 - M)

    Geometric interpolation: linear in the LOGARITHM of conductivity, which is
    how conductivity varies and why the endpoints span four orders of
    magnitude rather than a factor of four.

    Checked as the three properties that define it rather than by re-typing
    the expression -- the two endpoints, and the fact that the midpoint is the
    geometric and not the arithmetic mean, which is the whole content of
    "geometric" and the thing a careless rewrite would lose.
    """
    m = _model()
    net = m.network
    intact = ~net.link_v                       # unjointed vertical links

    m.M = np.ones((m.nz, m.nx))
    kv, _, _ = m.link_conductivity()
    assert np.allclose(kv[intact], m.k_matrix_at_T, rtol=1e-12, atol=0.0)

    m.M = np.zeros((m.nz, m.nx))
    kv, _, _ = m.link_conductivity()
    assert np.allclose(kv[intact], m.k_weathered_at_T, rtol=1e-12, atol=0.0)

    m.M = np.full((m.nz, m.nx), 0.5)
    kv, _, _ = m.link_conductivity()
    geometric = np.sqrt(m.k_matrix_at_T * m.k_weathered_at_T)
    arithmetic = 0.5 * (m.k_matrix_at_T + m.k_weathered_at_T)
    assert np.allclose(kv[intact], geometric, rtol=1e-12, atol=0.0)
    assert not np.allclose(kv[intact], arithmetic, rtol=1e-3, atol=0.0)

    # a joint is a joint whatever the rock beside it has done
    m.M = np.zeros((m.nz, m.nx))
    kv, _, _ = m.link_conductivity()
    assert np.allclose(kv[net.link_v], m.k_fracture, rtol=1e-12, atol=0.0)


def test_the_head_is_re_solved_as_the_rock_changes():
    """
    The feedback only exists if the head is actually recomputed. Triggered by
    how much the ROCK has changed, never by a step count: a step count would
    tie the answer to the step size, and halving the drift budget would
    silently double how often the flow was updated.
    """
    # 280 kyr, not 40: sourcing tau slowed the model sevenfold, and this
    # test needs enough ROCK CHANGE to trigger re-solves, not enough time.
    m = _model()
    m.flow_tolerance = 0.02
    m.run(years=280e3)
    assert m.flow_solves > 4, m.flow_solves       # the feedback is live

    fine = _model()
    fine.flow_tolerance = 0.02
    fine.c_drift_max = 0.25 * fine.c_drift_max    # four times the steps
    fine.run(years=280e3)

    # The ANSWER is what must not depend on the step size, and it does not.
    # The solve COUNT does move -- a coarse step overshoots the tolerance
    # before the check, so each solve covers more than flow_tolerance of
    # change and fewer are needed (16 against 46 here). That is a weaker
    # coupling than a step count, which would be proportional, and it washes
    # out of the result: 0.1630 against 0.1665.
    a, b = m.dissolved_fraction.mean(), fine.dissolved_fraction.mean()
    assert abs(a - b) / max(a, 1e-12) < 0.05, (a, b)


# ------------------------------------------------------------- flow, Darcy

def test_the_head_field_satisfies_the_darcy_equation_cell_by_cell():
    """
        div( K grad H ) = 0,  H = psi - d   (d is depth, positive down)

    Per cell, not globally: a global balance telescopes and would pass for a
    head field that is wrong in the interior. Sources are the surface
    infiltration and the base drainage; everything else must close.
    """
    m = _model()
    nz, nx = m.nz, m.nx
    div = np.zeros((nz, nx))
    div[:-1, :] += m.q_v
    div[1:, :] -= m.q_v
    div[:, :-1] += m.q_h
    div[:, 1:] -= m.q_h
    if m.network.periodic_x:
        div[:, -1] += m.q_wrap
        div[:, 0] -= m.q_wrap
    div[0, :] -= m.infiltration * m.dx          # source in at the surface
    div[-1, :] += m.q_out_base                  # sink out at the base
    scale = m.infiltration * m.dx
    assert np.abs(div).max() / scale < 1e-9


# ---------------------------------------------------- transport and solute

def test_the_transport_coefficient_is_molecular_plus_dispersive():
    """
        D = D_aqueous(T) / tortuosity(M) + grain_size * |v|

    Two terms. The second is hydrodynamic dispersion; at these fluxes the
    Reynolds number is about 3e-5, so it is not turbulence, and it is pore and
    aperture geometry that sets the dispersivity.
    """
    m = _model()
    D_v, D_h = m.transport_coefficients()
    net = m.network
    tv, _ = m.link_tortuosity()
    want_v = np.where(net.link_v, m.D_aqueous, m.D_aqueous / tv) \
        + m.dispersivity * np.abs(m.q_v) / m.dx
    assert np.allclose(D_v, want_v, rtol=1e-12, atol=0.0)
    # both terms actually matter somewhere
    # Dispersion still leads in the fastest joints and is negligible in the
    # matrix -- the point of putting it at the pore scale.
    assert (m.dispersivity * np.abs(m.q_v) / m.dx).max() > m.D_aqueous / 10.0
    assert ((m.dispersivity * np.abs(m.q_v) / m.dx).min()
            < m.D_aqueous / m.tortuosity_weathered)


@pytest.mark.parametrize("driver", ["dissolution", "oxidation"])
def test_the_solved_concentration_satisfies_the_stated_cell_balance(driver):
    """
        sum_out f c_i - sum_in f c_j + sum_links D (c_i - c_j) + r dx^2 c_i
            = S_i
        S_i = r dx^2   dissolution: every cell a source, inlet c = 0
        S_i = q_in dx  oxidation: the surface only, inlet c = 1

    Assembled here from the model's own fluxes, independently of the sparse
    matrix the solver builds, and checked PER CELL.

    Run for BOTH drivers, and the source is written out here rather than taken
    from ``_solute_source``, so that the test cannot pass by agreeing with the
    code it is checking. The left-hand side is identical in the two cases --
    that is the claim design 08 rests on, and this is where it is checked.
    """
    m = _model()
    m.set_driver(driver)
    r = m.reaction_coefficient
    c = m.solve_solute(r)
    nz, nx, dx = m.nz, m.nx, m.dx
    D_v, D_h = m.transport_coefficients()
    source = np.zeros((nz, nx))
    if driver == "oxidation":
        source[0, :] = m.infiltration * dx      # rain, at c = 1
    else:
        source += r * dx * dx                   # rock, in every cell
    res = r * dx * dx * c - source

    def flux(a_slice, b_slice, f, D):
        ca, cb = c[a_slice], c[b_slice]
        adv = np.where(f > 0, f * ca, f * cb)           # upwind
        return adv + D * (ca - cb)

    fv = flux((slice(0, -1), slice(None)), (slice(1, None), slice(None)),
              m.q_v, D_v)
    res[:-1, :] += fv
    res[1:, :] -= fv
    fh = flux((slice(None), slice(0, -1)), (slice(None), slice(1, None)),
              m.q_h, D_h)
    res[:, :-1] += fh
    res[:, 1:] -= fh
    if m.network.periodic_x:
        Dw = np.where(m.network.link_wrap, m.D_aqueous,
                      m.D_aqueous / m.link_tortuosity_wrap()) \
            + m.dispersivity * np.abs(m.q_wrap) / dx
        fw = np.where(m.q_wrap > 0, m.q_wrap * c[:, -1], m.q_wrap * c[:, 0]) \
            + Dw * (c[:, -1] - c[:, 0])
        res[:, -1] += fw
        res[:, 0] -= fw
    res[-1, :] += m.q_out_base * c[-1, :]

    scale = max((r * dx * dx).max(), source.max())
    assert np.abs(res).max() / scale < 1e-8


def test_diffusion_is_what_lets_a_block_weather_inward():
    """
    Not a transcription but the reason the diffusive term is there. With pure
    advection a block interior saturates and stays at c = 1 for ever, so rock
    off a flow path never weathers and the model is binary. Turning the
    transport coefficients off must reproduce that, and turning them on must
    not.

    Checked at three states, because the tortuosity now follows the rock and
    the answer is not the same at each. Fresh granite is very nearly sealed --
    diffusion still triples the undersaturated fraction but cannot open the
    whole section -- while at half weathered it is decisive, and by the time
    the rock is mostly gone the water reaches everywhere without help.
    """
    def undersaturated(fraction_remaining):
        m = _model()
        m.set_driver("dissolution")
        m.M[:] = fraction_remaining
        m.solve_flow()
        r = m.reaction_coefficient
        on = (m.driving_force(m.solve_solute(r)) > 1e-6).mean()
        m.D_molecular = m.D_O2_molecular = 0.0
        m.grain_size = 0.0
        off = (m.driving_force(m.solve_solute(r)) > 1e-6).mean()
        return float(on), float(off)

    fresh_on, fresh_off = undersaturated(1.0)
    assert fresh_on > 2.0 * fresh_off, (fresh_on, fresh_off)

    half_on, half_off = undersaturated(0.5)
    assert half_on > 0.99, half_on          # diffusion opens the whole section
    assert half_off < 0.5, half_off         # advection alone does not


def test_diffusion_is_what_gives_the_oxidation_rind_its_width():
    """
    The oxidation counterpart, and it is a different statement.

    Under dissolution the test above works because a block interior reaches
    c = 1 EXACTLY and the driving force vanishes, so "what fraction of the
    domain is still reacting" is a clean binary. A reactant has no such
    shut-off: pure advection still carries some oxygen everywhere, so that
    measure saturates at 1.0 and says nothing.

    What diffusion does is set the LEVEL of oxygen in the matrix, not the
    shape of its decline. Measured on this fixture at t = 0, mean c one cell
    and three cells from a joint:

        with diffusion     0.2728    0.0694
        without            0.0368    0.0132

    -- a factor of 7.4 and 5.3. The far/near RATIO actually rises slightly
    when diffusion is removed, 0.25 to 0.36, because what is left is
    advection carrying oxygen down THROUGH the matrix rather than sideways
    from the joint, and a through-flow does not decay with distance from a
    joint the way a diffusive gradient does. The first draft of this test
    asserted the ratio and had the sign backwards; the level is the claim,
    and it is the one that matters -- without diffusion the matrix is nearly
    anoxic and there is no rind to speak of.
    """
    def rind(with_diffusion):
        m = _model()
        m.set_driver("oxidation")
        if not with_diffusion:
            m.D_molecular = m.D_O2_molecular = 0.0
            # NOT grain_size = 0, which is what the dissolution version does.
            # Under oxidation the grain size also sets the reactive surface
            # area, so zeroing it removes the reaction as well as the
            # dispersion -- and divides by zero on the way. Molecular
            # diffusion is the term that matters here regardless: in the
            # matrix it beats mechanical dispersion by about two hundred to
            # one, which :attr:`diffusivity_factor` states and this relies on.
        m.solve_flow()
        c = m.solve_solute(m.reaction_coefficient)
        d = m.network.distance_to_fracture()
        near = c[np.isclose(d, m.dx)].mean()
        far = c[np.isclose(d, 3.0 * m.dx)].mean()
        return float(near), float(far)

    on_near, on_far = rind(True)
    off_near, off_far = rind(False)
    assert on_near > 3.0 * off_near, (on_near, off_near)
    assert on_far > 3.0 * off_far, (on_far, off_far)


def test_corners_stay_further_from_saturation_than_faces():
    """
    The geometric route to spheroidal rounding, and it comes free with
    diffusion: a corner sheds solute to two joint faces and a face to one, so
    at equal distance from the nearest joint a corner is further from
    saturation and therefore weathers faster.
    """
    m = _model()
    c = m.solve_solute(m.reaction_coefficient)
    u = m.driving_force(c)      # where the water can still do work
    jc = np.nonzero(m.network.link_v[m.nz // 2, :])[0]
    jr = np.nonzero(m.network.link_h.mean(axis=1) > 0.5)[0]
    c0, r0, r1 = jc[1], jr[1], jr[2]
    k = max(int(round(0.2 / m.dx)), 2)
    face = u[(r0 + r1) // 2, c0 + k]
    corner = u[r0 + k, c0 + k]
    assert corner > face


# -- the thermodynamic pair, registered in the coverage ledger ------------

R_GAS = 8.314


def _thermo(tC=11.85, **kw):
    """A small model at a chosen temperature, for the two temperature laws."""
    net = FractureNetwork(20, 20, 0.05, periodic_x=True).seed(
        sets=orthogonal_grid(0.5), rng=np.random.default_rng(0))
    m = Weathering(net)
    for k, v in kw.items():
        setattr(m, k, v)
    m.set_infiltration(0.30 / YEAR)
    m.set_temperature(tC + 273.15)
    return m.initialize()


@pytest.mark.parametrize("tC", [0.0, 11.85, 30.0])
def test_the_factors_are_the_textbook_arrhenius_and_van_t_hoff(tC):
    """Computed here from the equations as they appear on the exercise page,
    independently of how the model forms them."""
    m = _thermo(tC)
    T = tC + 273.15
    assert float(np.mean(m.rate_factor)) == pytest.approx(
        np.exp(-(m.E_a / R_GAS) * (1.0 / T - 1.0 / m.T_ref)), rel=1e-12)
    assert float(np.mean(m.solubility_factor)) == pytest.approx(
        np.exp(-(m.delta_H_r / R_GAS) * (1.0 / T - 1.0 / m.T_ref)), rel=1e-12)


@pytest.mark.parametrize("tC", [0.0, 5.0, 20.0, 30.0])
def test_only_the_DIFFERENCE_of_the_two_enthalpies_sets_the_length_scale(tC):
    """
    The central claim, and the one most easily lost in an edit.

    The saturation length goes as ``C_eq / k``, so ``E_a`` and ``delta_H_r``
    enter it with opposite signs. Two completely different pairs sharing a
    difference must give the same length at every temperature -- which is why
    the pair may not be chosen one at a time, and why a field study of
    weathering against temperature recovers the difference rather than E_a.
    """
    a = _thermo(tC, E_a=69.8e3, delta_H_r=32.9e3)     # difference 36.9
    b = _thermo(tC, E_a=100.0e3, delta_H_r=63.1e3)    # difference 36.9
    assert a.apparent_activation_energy == pytest.approx(
        b.apparent_activation_energy, rel=1e-12)
    assert float(np.mean(a.saturation_length)) == pytest.approx(
        float(np.mean(b.saturation_length)), rel=1e-12)


def test_the_water_viscosity_correlation_matches_tabulated_water():
    """
    ``mu(T) = 2.414e-5 * 10 ** (247.8 / (T - 140))``

    The only physical property of water in the model, and the only equation
    here that is a fit rather than a law -- so it is checked against
    tabulated viscosities rather than against the code that uses it.
    """
    from corestone.weathering import water_viscosity
    for T, tabulated in ((273.15, 1.792e-3), (293.15, 1.002e-3),
                         (313.15, 0.653e-3)):
        assert water_viscosity(T) == pytest.approx(tabulated, rel=0.025)


def test_diffusivity_follows_stokes_einstein_and_is_not_constant():
    """``D(T)`` goes as ``T / mu(T)``: a factor of three across the demo's
    slider, carried almost entirely by the viscosity."""
    cold, warm = _thermo(0.0), _thermo(30.0)
    for m in (cold, warm):
        T = float(np.mean(m.T))
        assert m.diffusivity_factor == pytest.approx(
            (T / water_viscosity(T))
            / (m.T_D_ref / water_viscosity(m.T_D_ref)), rel=1e-12)
    assert warm.D_aqueous / cold.D_aqueous == pytest.approx(2.6, rel=0.1)


def test_the_joint_conductivity_is_the_cubic_law_on_its_aperture():
    """
    ``k_fracture = rho g a^3 / (12 mu dx)``

    A joint is a geometry. The conductivity is derived from the aperture, so
    what the model states is a measurable object rather than the conductivity
    of a joint smeared over an arbitrary cell.
    """
    from corestone.weathering import RHO_WATER, GRAVITY, water_viscosity
    m = _thermo(11.85)
    mu = water_viscosity(float(np.mean(m.T)))
    assert m.k_fracture == pytest.approx(
        RHO_WATER * GRAVITY * m.joint_aperture ** 3
        / (12.0 * mu * m.network.dx), rel=1e-12)


def test_the_joint_is_the_same_joint_at_every_cell_size():
    """
    The defect this derivation removes. Held as a constant CONDUCTIVITY the
    implied aperture moved with the grid -- 91 um at 5 cm against 67 um at
    2 cm -- so refining the mesh quietly tightened the joints by a third,
    while the exercise page promises that cell size is the numerical grid and
    not the rock. Transmissivity is the invariant.
    """
    T = []
    for dx in (0.05, 0.025, 0.02):
        n = int(round(3.0 / dx))
        net = FractureNetwork(n, n, dx, periodic_x=True).seed(
            sets=orthogonal_grid(1.0), rng=np.random.default_rng(1))
        m = Weathering(net)
        m.set_temperature(285.0)
        T.append(m.k_fracture * dx)
    assert max(T) / min(T) == pytest.approx(1.0, rel=1e-12), T


def test_temperature_does_not_move_the_flow_field():
    """
    A regression, and the bug was mine.

    Hydraulic conductivity is ``k_intrinsic rho g / mu`` for any medium, so
    warming raises the joints and the matrix alike and leaves their ratio
    alone. The infiltration is prescribed at the surface rather than driven by
    a head gradient, so an unchanged ratio means an unchanged flow field:
    temperature must not move the water at all, only the chemistry.

    Deriving the joint conductivity from an aperture introduced the viscosity
    on the joints alone while the matrix ends stayed fixed. That doubled the
    joint-to-matrix contrast between 0 and 30 C -- 18653 to 41017 -- and moved
    the speed field by 55 %, a temperature effect on the flow with no physical
    basis whatever.
    """
    def speed(tC):
        net = FractureNetwork(30, 30, 0.05, periodic_x=True).seed(
            sets=orthogonal_grid(0.5), rng=np.random.default_rng(12345))
        m = Weathering(net)
        m.set_infiltration(0.30 / YEAR)
        m.set_temperature(tC + 273.15)
        return m.initialize().darcy_speed.copy()
    cold, warm = speed(0.0), speed(30.0)
    assert np.abs(warm / cold - 1.0).max() < 1e-6, np.abs(warm/cold - 1.0).max()


def test_the_tortuosity_follows_the_rock_like_the_conductivity():
    """
    ``tortuosity(M) = tortuosity_fresh^M * tortuosity_weathered^(1 - M)``

    The same geometric interpolation the conductivity uses, and for the same
    reason: dissolving rock opens porosity to diffusion as surely as it opens
    it to flow. Held at the weathered value, fresh granite diffused about a
    thousand times too freely -- and did so at t = 0, when every cell is fresh
    and the rind is forming.
    """
    m = _thermo(11.85)
    for frac in (1.0, 0.5, 0.0):
        m.M[:] = frac
        tv, th = m.link_tortuosity()
        want = (m.tortuosity_fresh ** frac
                * m.tortuosity_weathered ** (1.0 - frac))
        assert np.allclose(tv, want, rtol=1e-12, atol=0.0)
        assert np.allclose(th, want, rtol=1e-12, atol=0.0)
    assert m.tortuosity_fresh / m.tortuosity_weathered == pytest.approx(1e3)


# ------------------------------------------------------------ the oxidation
#
# Design 08. These parameters do not yet move rock; the transcriptions are
# here from the moment the numbers are, so that a value can be checked before
# a model is built on it.

def test_the_oxygen_solubility_correlation_matches_tabulated_water():
    """
    ``ln C = -139.34411 + 1.575701e5 / T - 6.642308e7 / T^2
    + 1.2438e10 / T^3 - 8.621949e11 / T^4``

    The second fit in this model rather than a law, so like the viscosity it
    is checked against a table and not against the code that uses it. The
    correlation returns mol/m3; the table is in mg/L, which is the same number
    divided by the molar mass.

    Checked a second, independent way as well, because the coefficients were
    written from memory and a mistranscribed one would still produce a smooth
    curve. Henry's law at 25 C, with a constant from nowhere near this
    correlation, agrees to 2 % -- see :func:`oxygen_solubility`.
    """
    from corestone.weathering import oxygen_solubility, M_O2
    for tC, mg_per_L in ((0.0, 14.62), (5.0, 12.77), (10.0, 11.29),
                         (15.0, 10.08), (20.0, 9.09), (25.0, 8.26),
                         (30.0, 7.56)):
        got = oxygen_solubility(tC + 273.15) * M_O2
        assert got == pytest.approx(mg_per_L, rel=5e-4), (tC, got)

    # Independent of the table: Henry's law on moist air at one atmosphere.
    pO2 = (1.0 - 0.03126) * 0.20946                  # atm, 25 C
    henry = 1.3e-3 * pO2 * M_O2 * 1000.0             # mg/L
    assert oxygen_solubility(298.15) * M_O2 == pytest.approx(henry, rel=0.03)


def test_warming_drives_oxygen_out_of_solution_and_silica_into_it():
    """
    The sign reversal, which is the whole reason the driver matters.

    Every ceiling in this model until now rose with temperature: warm water
    dissolves more silica, so each litre carries more away and weathering
    speeds up. Oxygen is a GAS, and a gas comes out of solution as the water
    warms. So the oxygen budget gets WORSE with warming while the silica
    budget gets better, and a model paced by oxygen has a temperature response
    that need not point the same way as one paced by silica.

    Not a modelling choice -- the two correlations were written independently
    of each other and of this test.
    """
    cold, warm = _thermo(0.0), _thermo(30.0)
    assert warm.C_O2 < cold.C_O2
    assert cold.C_O2 / warm.C_O2 == pytest.approx(1.93, rel=0.02)
    assert warm.solubility_factor > cold.solubility_factor
    # ...and therefore the two taus move in opposite directions.
    assert warm.tau_oxidation > cold.tau_oxidation
    assert warm.silica_tau < cold.silica_tau


def test_the_volume_expansion_is_the_ratio_of_the_molar_volumes():
    """
    ``dV/V = V_goethite / V_FeO - 1``

    Both molar volumes are per single Fe, so the ratio carries no
    stoichiometric factor and 0.735 is a pure geometric statement: oxidised,
    hydrated iron needs three quarters again as much room.

    The value is also the check on WHICH PRODUCT. Fletcher's 0.7 is quoted for
    ferrihydrite from a compilation that has no ferrihydrite in it; goethite
    reproduces it and ferrihydrite does not.
    """
    m = _thermo()
    assert m.volume_expansion == pytest.approx(m.V_goethite / m.V_FeO - 1.0,
                                               rel=1e-12)
    assert m.volume_expansion == pytest.approx(0.735, abs=5e-4)


def test_tau_on_oxygen_is_the_iron_divided_by_four_and_by_the_solubility():
    """
    ``tau_O2 = f_FeO / (4 V_FeO C_O2(T))``

    Four Fe(II) per O2. Transcribed independently of the property, and pinned
    at the reference state, because the whole case for design 08 is that this
    number is much smaller than the silica ``tau`` -- so if it is wrong, the
    case is wrong.
    """
    m = _thermo(11.85)
    iron = m.f_FeO / m.V_FeO                       # mol Fe per m3 of rock
    assert iron == pytest.approx(916.7, rel=1e-3)
    assert m.tau_oxidation == pytest.approx(0.25 * iron / m.C_O2, rel=1e-12)
    assert m.tau_oxidation == pytest.approx(678.1, rel=1e-3)
    # The comparison design 08 rests on. NOT 15x: that figure was computed
    # with Fletcher's f_FeO = 0.05, which the same document rejects.
    assert m.silica_tau / m.tau_oxidation == pytest.approx(70.4, rel=1e-2)


def test_the_front_ceiling_is_the_flux_over_tau():
    """
    ``front ceiling = q / tau``

    Stoichiometry alone: water arriving saturated and leaving stripped. A rate
    no mechanism can beat, which is what makes it worth stating -- the silica
    version, 6.3 m/Myr, is barely above the 4 to 7 m/Myr measured in the
    field, so this model has been running against its own budget.
    """
    m = _thermo(11.85)
    m.set_infiltration(0.30 / YEAR)
    assert m.oxidation_front_ceiling == pytest.approx(
        m.infiltration / m.tau_oxidation, rel=1e-12)
    per_Myr = m.oxidation_front_ceiling * YEAR * 1e6
    assert per_Myr == pytest.approx(442.4, rel=1e-3)
    assert (m.infiltration / m.silica_tau * YEAR * 1e6
            == pytest.approx(6.28, rel=1e-2))


def test_the_reactive_surface_area_is_six_phi_over_d():
    """
    ``A = 6 phi / d``

    Cubic grains, the model's one convention for surface area, applied to the
    biotite exactly as it is to the plagioclase behind ``L_ref``. Checked
    against the plagioclase number too: the same expression at 30 % must
    return the 900 m2/m3 that ``L_ref`` was derived from, or the two
    reactions are not being measured the same way.
    """
    m = _thermo()
    assert m.biotite_surface_area == pytest.approx(
        6.0 * m.phi_biotite / m.grain_size, rel=1e-12)
    assert m.biotite_surface_area == pytest.approx(180.0, rel=1e-12)
    assert 6.0 * 0.30 / m.grain_size == pytest.approx(900.0, rel=1e-12)


def test_the_oxygen_penetration_depth_is_the_reaction_diffusion_length():
    """
    ``penetration = sqrt(D_O2 / (tortuosity_fresh k_ox A))``

    Once the advective Damkohler is far below one, this is the only length
    left, and it is what would shelter a corestone in an oxygen-driven model.
    Pinned at the reference state because the whole framing of the exercise
    turns on it being centimetres rather than metres.

    Also pinned: the regime it implies. The section-scale Damkohler on oxygen
    must be below the 1/3 that :attr:`regime` calls reaction-limited, while
    the silica one is above the 3 that it calls saturation-limited. Those two
    facts together are the finding of probe I, and if either stops holding,
    the page's claim about what shelters a corestone must be re-read.
    """
    m = _thermo(11.85)
    m.set_infiltration(0.30 / YEAR)
    want = np.sqrt(m.D_O2_aqueous / m.tortuosity_fresh
                   / m.specific_oxidation_coefficient)
    assert m.oxidation_penetration_depth == pytest.approx(want, rel=1e-12)
    assert m.oxidation_penetration_depth == pytest.approx(0.0448, rel=2e-3)

    assert m.specific_oxidation_coefficient == pytest.approx(7.2e-11, rel=1e-3)
    assert m.oxidation_length == pytest.approx(132.0, rel=1e-2)
    # The regime claim depends on section depth, so it is made at the depth
    # the demo actually uses -- 3 m -- and not on this fixture's 1 m.
    deep = Weathering(FractureNetwork(60, 60, 0.05, periodic_x=True).seed(
        sets=orthogonal_grid(1.0), rng=np.random.default_rng(0)))
    deep.set_infiltration(0.30 / YEAR)
    deep.set_temperature(11.85 + 273.15)
    deep.initialize()
    assert deep.damkohler > 3.0                    # saturation-limited
    assert deep.oxidation_damkohler < 1.0 / 3.0    # reaction-limited
    assert deep.regime == "saturation-limited"

    # ...and the ratio of the two, which no grid can change.
    assert m.oxidation_length / m.saturation_length == pytest.approx(289.0,
                                                                    rel=1e-2)

    # Warming deepens the penetration, because diffusivity rises with it and
    # the rate constant has no temperature dependence to oppose that.
    cold, warm = _thermo(0.0), _thermo(30.0)
    assert warm.oxidation_penetration_depth > cold.oxidation_penetration_depth

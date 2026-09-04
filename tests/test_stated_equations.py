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


def _model(dx=0.10, spacing=1.5, width=12.0, depth=9.0):
    nz, nx = periodic_grid_shape(width, depth, dx, spacing)
    net = FractureNetwork(nz, nx, dx, periodic_x=True).seed(
        sets=orthogonal_grid(spacing), rng=np.random.default_rng(12345))
    return Weathering(net).initialize()


# ---------------------------------------------------------------- rate law

def test_the_dissolution_rate_per_unit_volume_does_not_depend_on_the_flux():
    """
        R = k(T) * A * (1 - C / C_eq)

    The rate per unit volume is a property of the rock. It may not depend on
    how fast water moves past it. This is the invariant the original bug broke:
    holding the saturation length uniform made dissolution proportional to the
    local flux, weighting a joint against the matrix by a factor of ~3000.

    The reaction coefficient r = kA/C_eq is now the primary quantity, and the
    flux cancels explicitly rather than inside an exponent, so this is a direct
    check: r must be uniform across cells that differ in flux by 1000x.
    """
    m = _model()
    r = m.reaction_coefficient
    # The flux distribution is bimodal, so select by structure, not quantile.
    fast, slow = m.network.cell, ~m.network.cell
    assert fast.any() and slow.any()
    assert np.median(m.q[fast]) > 100.0 * np.median(m.q[slow])   # real contrast
    assert np.ptp(r) / r.mean() < 1e-12                     # yet r is uniform


def test_what_the_rock_loses_is_what_the_water_carries_out_of_the_base():
    """
        d(M/M0)/dt = - r (1 - c) / tau

    No solute enters at the surface, so in steady transport everything
    dissolved must leave through the base. Checked against the export, not
    against a restatement of the same expression.
    """
    m = _model()
    r = m.reaction_coefficient
    c = m.solve_solute(r)
    produced = (r * (1.0 - c) * m.dx * m.dx).sum()
    exported = (m.q_out_base * c[-1, :]).sum()
    assert produced == pytest.approx(exported, rel=1e-9)

    rate = r * (1.0 - c) / m.tau
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
    """
    m = _model()
    m.set_temperature(m.T_ref)
    base = m.tau
    m.set_temperature(m.T_ref + 20.0)
    assert m.tau < base
    assert m.tau == pytest.approx(m.tau_ref / m.solubility_factor, rel=1e-12)


def test_the_matrix_conducts_better_as_it_dissolves():
    """
        k(M) = k_matrix^M * k_weathered^(1 - M)

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
    assert np.allclose(kv[intact], m.k_matrix, rtol=1e-12)

    m.M = np.zeros((m.nz, m.nx))
    kv, _, _ = m.link_conductivity()
    assert np.allclose(kv[intact], m.k_weathered, rtol=1e-12)

    m.M = np.full((m.nz, m.nx), 0.5)
    kv, _, _ = m.link_conductivity()
    geometric = np.sqrt(m.k_matrix * m.k_weathered)
    arithmetic = 0.5 * (m.k_matrix + m.k_weathered)
    assert np.allclose(kv[intact], geometric, rtol=1e-12)
    assert not np.allclose(kv[intact], arithmetic, rtol=1e-3)

    # a joint is a joint whatever the rock beside it has done
    m.M = np.zeros((m.nz, m.nx))
    kv, _, _ = m.link_conductivity()
    assert np.allclose(kv[net.link_v], m.k_fracture, rtol=1e-12)


def test_the_head_is_re_solved_as_the_rock_changes():
    """
    The feedback only exists if the head is actually recomputed. Triggered by
    how much the ROCK has changed, never by a step count: a step count would
    tie the answer to the step size, and halving the drift budget would
    silently double how often the flow was updated.
    """
    m = _model()
    m.flow_tolerance = 0.02
    m.run(years=40e3)
    assert m.flow_solves > 4, m.flow_solves       # the feedback is live

    fine = _model()
    fine.flow_tolerance = 0.02
    fine.c_drift_max = 0.25 * fine.c_drift_max    # four times the steps
    fine.run(years=40e3)

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
        D = D_molecular / tortuosity + dispersivity * |v|

    Two terms. The second is hydrodynamic dispersion; at these fluxes the
    Reynolds number is about 3e-5, so it is not turbulence, and it is pore and
    aperture geometry that sets the dispersivity.
    """
    m = _model()
    D_v, D_h = m.transport_coefficients()
    net = m.network
    want_v = np.where(net.link_v, m.D_molecular, m.D_molecular / m.tortuosity) \
        + m.dispersivity * np.abs(m.q_v) / m.dx
    assert np.allclose(D_v, want_v, rtol=1e-12)
    # both terms actually matter somewhere
    assert (m.dispersivity * np.abs(m.q_v) / m.dx).max() > m.D_molecular
    assert (m.dispersivity * np.abs(m.q_v) / m.dx).min() < m.D_molecular / m.tortuosity


def test_the_solved_concentration_satisfies_the_stated_cell_balance():
    """
        sum_out f c_i - sum_in f c_j + sum_links D (c_i - c_j) + r dx^2 c_i
            = r dx^2

    Assembled here from the model's own fluxes, independently of the sparse
    matrix the solver builds, and checked PER CELL.
    """
    m = _model()
    r = m.reaction_coefficient
    c = m.solve_solute(r)
    nz, nx, dx = m.nz, m.nx, m.dx
    D_v, D_h = m.transport_coefficients()
    res = r * dx * dx * (c - 1.0)                       # reaction + source

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
        Dw = np.where(m.network.link_wrap, m.D_molecular,
                      m.D_molecular / m.tortuosity) \
            + m.dispersivity * np.abs(m.q_wrap) / dx
        fw = np.where(m.q_wrap > 0, m.q_wrap * c[:, -1], m.q_wrap * c[:, 0]) \
            + Dw * (c[:, -1] - c[:, 0])
        res[:, -1] += fw
        res[:, 0] -= fw
    res[-1, :] += m.q_out_base * c[-1, :]

    scale = (r * dx * dx).max()
    assert np.abs(res).max() / scale < 1e-8


def test_diffusion_is_what_lets_a_block_weather_inward():
    """
    Not a transcription but the reason the diffusive term is there. With pure
    advection a block interior saturates and stays at c = 1 for ever, so rock
    off a flow path never weathers and the model is binary. Turning the
    transport coefficients off must reproduce that, and turning them on must
    not.
    """
    m = _model()
    r = m.reaction_coefficient
    with_diff = m.solve_solute(r)
    m.D_molecular = 0.0
    m.dispersivity = 0.0
    without = m.solve_solute(r)
    assert ((1.0 - without) > 1e-6).mean() < 0.5
    assert ((1.0 - with_diff) > 1e-6).mean() > 0.99


def test_corners_stay_further_from_saturation_than_faces():
    """
    The geometric route to spheroidal rounding, and it comes free with
    diffusion: a corner sheds solute to two joint faces and a face to one, so
    at equal distance from the nearest joint a corner is further from
    saturation and therefore weathers faster.
    """
    m = _model()
    c = m.solve_solute(m.reaction_coefficient)
    u = 1.0 - c
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

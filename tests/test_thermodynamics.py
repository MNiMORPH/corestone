"""
The thermodynamic identities, stated as tests so they cannot quietly stop
being true.

These are not regression tests for bugs: each one is a statement about the
model's thermodynamics that a student should be able to check, written in a
form the machine also checks. If one of them fails, the physics changed.
"""

import numpy as np
import pytest

from corestone import FractureNetwork, Weathering, orthogonal_grid, YEAR

R_GAS = 8.314


def model(tC=11.85, **kw):
    net = FractureNetwork(20, 20, 0.05, periodic_x=True).seed(
        sets=orthogonal_grid(0.5), rng=np.random.default_rng(0))
    m = Weathering(net)
    for k, v in kw.items():
        setattr(m, k, v)
    m.set_infiltration(0.30 / YEAR)
    m.set_temperature(tC + 273.15)
    m.initialize()
    return m


def test_both_temperature_factors_are_exactly_one_at_the_reference():
    """The reference state is a definition, not an approximation. This is why
    the demo's temperature slider opens where it does."""
    m = model()
    m.set_temperature(m.T_ref)
    assert float(np.mean(m.rate_factor)) == pytest.approx(1.0, abs=1e-12)
    assert float(np.mean(m.solubility_factor)) == pytest.approx(1.0, abs=1e-12)
    assert float(np.mean(m.saturation_length)) == pytest.approx(m.L_ref, rel=1e-12)


@pytest.mark.parametrize("tC", [0.0, 11.85, 30.0])
def test_the_factors_are_the_textbook_arrhenius_and_van_t_hoff(tC):
    """Computed here from the equations as they appear on the exercise page,
    independently of how the model forms them."""
    m = model(tC)
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
    a = model(tC, E_a=69.8e3, delta_H_r=32.9e3)     # difference 36.9
    b = model(tC, E_a=100.0e3, delta_H_r=63.1e3)    # difference 36.9
    assert a.apparent_activation_energy == pytest.approx(
        b.apparent_activation_energy, rel=1e-12)
    assert float(np.mean(a.saturation_length)) == pytest.approx(
        float(np.mean(b.saturation_length)), rel=1e-12)


def test_a_larger_reaction_enthalpy_than_activation_energy_reverses_temperature():
    """
    Warming does not have to accelerate weathering, and the model must not
    quietly assume it does.

    If the ceiling on the solute rises with temperature faster than the rate
    constant does, the saturation length GROWS with warming and water travels
    further before it stops working. That is what a kaolinite-buffered
    reading of calcic plagioclase gives, and it is a real regime; this model
    is not in it, but nothing here forbids it.
    """
    normal = [float(np.mean(model(t, E_a=69.8e3, delta_H_r=32.9e3)
                            .saturation_length)) for t in (0.0, 30.0)]
    assert normal[1] < normal[0], "warming should shorten it in this model"
    reversed_ = [float(np.mean(model(t, E_a=17.8e3, delta_H_r=100.0e3)
                               .saturation_length)) for t in (0.0, 30.0)]
    assert reversed_[1] > reversed_[0], "delta_H_r > E_a must lengthen it"


def test_the_temperature_dependence_does_not_depend_on_the_calibration():
    """``L_ref`` sets the scale and nothing else. Doubling it doubles every
    length and changes no ratio -- which is what lets it be a free choice
    without the thermodynamics being wrong."""
    a, b = model(30.0, L_ref=0.50), model(30.0, L_ref=1.00)
    ra = float(np.mean(a.saturation_length)) / a.L_ref
    rb = float(np.mean(b.saturation_length)) / b.L_ref
    assert ra == pytest.approx(rb, rel=1e-12)


def test_the_damkohler_number_counts_e_foldings_across_the_section():
    m = model()
    depth = m.network.nz * m.network.dx
    assert float(np.mean(m.damkohler)) == pytest.approx(
        depth / float(np.mean(m.saturation_length)), rel=1e-12)


def test_the_regime_belongs_to_the_SECTION_and_not_to_the_model():
    """
    Written the other way round first, asserting the model is always
    saturation-limited, and it failed -- correctly. The Damkohler number is
    depth divided by a length, so it is a property of how much rock you are
    looking at, not of the rock.

    The 3 m section the exercise ships gives Da = 6 and sits firmly in the
    transport limit, which is what makes corestones: water saturates before
    it has crossed, so block interiors are sheltered. Take the same granite,
    the same water and the same temperature and look at only 1 m of it, and
    Da = 2 -- the same physics, a weaker limit, and less shelter.
    """
    def da(cells, dx=0.05):
        net = FractureNetwork(cells, cells, dx, periodic_x=True).seed(
            sets=orthogonal_grid(0.5), rng=np.random.default_rng(0))
        m = Weathering(net)
        m.set_infiltration(0.30 / YEAR); m.set_temperature(285.0); m.initialize()
        return float(np.mean(m.damkohler)), m.regime

    da3, regime3 = da(60)                       # the exercise's 3 m section
    da1, regime1 = da(20)                       # 1 m of the same rock
    assert da3 == pytest.approx(6.0, rel=1e-9)
    assert regime3 == "saturation-limited"
    assert da1 == pytest.approx(2.0, rel=1e-9)
    assert regime1 == "mixed"
    assert da3 == pytest.approx(3.0 * da1, rel=1e-9)


def test_warming_pushes_further_into_the_transport_limit():
    """Only because E_a exceeds delta_H_r here. It is not a general truth
    about weathering; see the reversal test above."""
    cold, warm = model(0.0).damkohler, model(30.0).damkohler
    assert warm > cold


def test_the_thermo_report_states_every_governing_number():
    m = model()
    text = m.thermo_report()
    for token in ("E_a", "delta_H_r", "E_a - delta_H_r", "Damkohler",
                  "saturation length", "regime", "T_ref"):
        assert token in text, token

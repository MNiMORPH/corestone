"""
The thermodynamic identities, stated as tests so they cannot quietly stop
being true.

The two that transcribe a displayed equation live in
``test_stated_equations.py``, where the coverage ledger requires them; the
rest are here.

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
    # 3.0 m / L_ref, with L_ref now derived from grain size rather than
    # calibrated, so 0.457 m and Da = 6.56 where the round 0.50 gave 6.00.
    assert da3 == pytest.approx(3.0 / 0.457, rel=1e-3)
    assert regime3 == "saturation-limited"
    assert da1 == pytest.approx(1.0 / 0.457, rel=1e-3)
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
                  "saturation length", "T_ref", "driver"):
        assert token in text, token
    # ...and the regime by name, taken from the model rather than spelled
    # here: this fixture is a 1 m section, so it is "mixed", not the 3 m
    # demo's "saturation-limited".
    assert m.regime in text, (m.regime, text)


def test_the_thermo_report_puts_the_two_budgets_side_by_side():
    """
    The comparison design 08 turns on, printed where a student can read it
    rather than buried in a document. Silica caps the front at 6.3 m/Myr,
    which is barely above the 4-7 m/Myr the field measures, so the model has
    been running against its own stoichiometry; oxygen caps it at 442.
    """
    m = model()
    m.set_infiltration(0.30 / YEAR)
    text = m.thermo_report()
    for token in ("OXIDATION --", "DISSOLUTION --", "<== DRIVING",
                  "C_O2", "tau_O2", "tau, silica", "oxidation length",
                  "O2 penetration", "reaction-limited", "front ceiling"):
        assert token in text, token

    # The regression this exists for: when tau became driver-aware, the line
    # labelled "tau, silica" started printing the OXYGEN value, and so did
    # both front ceilings. A report that mislabels a number is worse than no
    # report, so the two are pinned apart.
    ox, diss = text.split("DISSOLUTION --")
    assert "678" in ox and "442" in ox, ox
    assert "47744" in diss and "6.28" in diss, diss
    assert "47744" not in ox and "678 " not in diss, text
    assert "SLOWS the oxidation" in text and "SPEEDS the dissolution" in text


def test_the_oxidation_drivers_whole_temperature_response_is_the_gas_law():
    """
    No Arrhenius term is in the oxidation rate constant, and none is measured
    -- verified across 55 local PDFs, of which nine mention activation energy
    or Arrhenius and NONE in an oxidation context. (The one apparent hit is
    Svante Arrhenius cited as an author, 1954, which is why one reads the
    sentence.) Hogg & Meads (1975), a dedicated Mossbauer kinetics study of
    exactly this reaction, contains "activation energ", "Arrhenius", "kJ" and
    "kcal" zero times in 4303 words.

    So the driver still has a temperature response, and it is entirely the van
    't Hoff enthalpy of dissolving oxygen. The sign is NEGATIVE, because a gas
    leaves solution as water warms, and there is no activation energy on the
    other side to cancel it -- so cold rock oxidises faster.

    Checked analytically here rather than by three runs to 90 %: tau_O2 goes
    as 1/C_O2 and carries the whole of it. The runs agree -- 653, 869 and
    1238 kyr at 0, 11.85 and 30 C, an apparent -14.6 kJ/mol against this
    quantity's -14.5 -- and are too slow to be a unit test.
    """
    m = model()
    m.set_driver("oxidation")
    dH = m.oxygen_dissolution_enthalpy
    assert dH < 0.0
    assert dH / 1e3 == pytest.approx(-14.5, abs=0.2)

    # ...and tau_O2 carries exactly that enthalpy, which is what makes it the
    # model's whole response: r has no temperature dependence at all.
    temps = np.array([273.15, 288.15, 303.15])
    taus = []
    for T in temps:
        m.set_temperature(T)
        taus.append(m.tau_oxidation)
        assert m.specific_oxidation_coefficient == \
            pytest.approx(m.k_oxidation * m.biotite_surface_area, rel=1e-12)
    slope = np.polyfit(1.0 / temps, np.log(1.0 / np.array(taus)), 1)[0]
    # 6 %, not 2 %: the solubility correlation is a five-term polynomial in
    # 1/T, not a straight line in van 't Hoff coordinates, so the effective
    # enthalpy depends on the interval fitted. Three points over 273-303 K
    # give -15.2 kJ/mol where the property's 41 points over 273-313 give
    # -14.5. Both are right; a single van 't Hoff enthalpy is the
    # approximation, and the tolerance should say so rather than hide it by
    # fitting the same interval the code does.
    assert -slope * R_GAS == pytest.approx(dH, rel=0.06)

    # The oxidation LENGTH, by contrast, does not move at all.
    lengths = []
    for T in temps:
        m.set_temperature(T)
        lengths.append(m.oxidation_length)
    assert max(lengths) == pytest.approx(min(lengths), rel=1e-12)
    assert m.apparent_activation_energy == 0.0


def test_the_two_drivers_disagree_about_whether_warm_means_weathered():
    """
    The sign reversal, stated as the comparison a student would make. Under
    dissolution the saturation length shortens with warming and each litre
    carries more away; under oxidation each litre carries LESS oxygen and
    nothing speeds up to compensate.
    """
    cold, warm = model(0.0), model(30.0)
    for m in (cold, warm):
        m.set_driver("oxidation")
    assert warm.tau > cold.tau                 # warm water brings less oxygen
    for m in (cold, warm):
        m.set_driver("dissolution")
    assert warm.tau < cold.tau                 # warm water carries more silica

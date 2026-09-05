#! /usr/bin/python3

"""
Dissolve granite along its joints, and see what is left.

The model is one equation. Dissolution runs at an Arrhenius rate constant
multiplied by how far the pore water is from equilibrium,

    R = k(T) * A * (1 - C / C_eq)

so water that has equilibrated stops weathering rock, however soluble the rock
and however warm the water. Fresh water descends the joints; the joints
therefore decide where weathering happens, and rock the water never reaches --
or reaches already saturated -- survives as a corestone.

Working in normalised concentration ``c = C / C_eq`` removes the need to assert
a solubility. The scale it sets is the **saturation length**

    saturation_length = q * C_eq / (k(T) * A)

the e-folding length of the approach to saturation -- *not* a distance at which
equilibrium is reached, because there is no equilibrium here. ``c`` approaches
1 asymptotically and never arrives.

THE THERMODYNAMICS, IN FULL, BECAUSE IT IS THE PART MOST OFTEN GOT WRONG
------------------------------------------------------------------------

Temperature enters twice, with opposite effects on the length scale, and the
second one is the larger here. Write both in the textbook form:

    k(T)    = k_0 exp(-E_a      / R_g T)          Arrhenius, on the RATE
    C_eq(T) = C_0 exp(-dH_r     / R_g T)          van 't Hoff, on the CEILING

Warming speeds the reaction. It also raises the ceiling -- which does not make
the rock dissolve faster where it stands, it lets each litre of water carry
more away before it stops working. Those are different things and they act on
different terms.

Now form the saturation length, the only length this chemistry has:

    L = q C_eq / (k A)
      = (q C_0 / (k_0 A)) exp( -(dH_r - E_a) / R_g T )
      = L_ref exp( +(E_a - dH_r) / R_g ( 1/T - 1/T_ref ) )      *

``E_a`` and ``dH_r`` enter with OPPOSITE SIGNS and only their difference
survives. That difference, ``E_a - dH_r``, is the apparent activation energy
of weathering -- it is what a field study measuring weathering rate against
temperature actually recovers, and it is NOT ``E_a``. Here it is
69.8 - 32.9 = 36.9 kJ/mol, so the length scale is about half as
temperature-sensitive as the rate constant alone would suggest. If ``dH_r``
exceeded ``E_a``, warming would LENGTHEN L and slow the weathering down. That
is a real regime and this model does not forbid it; see
:attr:`Weathering.apparent_activation_energy`.

Two things follow that are easy to miss:

**k_0 and C_0 are gone.** Step * absorbed them into ``L_ref``, and every rate
in the code is a ratio to the reference state, so the absolute pre-exponential
factors cancel and are never evaluated. This model has a NORMALISATION where a
research model would need a thermodynamics. That is why ``L_ref`` can be
calibrated freely without the chemistry being wrong -- it sets the scale and
nothing else -- and why ``C_eq`` never appears alone anywhere below.

**C_eq still enters a second time, and not through L.** The rock must supply
``tau = M0 / C_eq`` volumes of saturated water per volume dissolved, so a
warmer, more soluble fluid needs fewer of them. That term carries ``dH_r``
alone, with no ``E_a`` to cancel against it. So temperature is not one dial:
it moves where weathering happens (through L) and how much water the job takes
(through tau), and the two do not have to point the same way.

Which limit the section is in is one dimensionless number, the Damkohler group
``Da = depth / L``, counting the e-foldings of saturation a parcel undergoes
on the way down. At the reference state it is 6: water reaches the base within
exp(-6) of saturation, so the section is firmly SATURATION-LIMITED, and that is
precisely what shelters a block interior and makes a corestone. In the
opposite limit water crosses barely touched, every block dissolves at the same
rate, and nothing is sheltered. Corestones are a saturation-limited phenomenon.

(Physical chemistry calls that limit *transport-limited*, for the transport of
solute. The word is avoided here: in geomorphology it means an erosion rate
set by the capacity to move sediment, and nothing in this model moves
sediment.)

NOTHING IN THE CHEMISTRY OR THE FLOW IS FITTED, AND THE TIMESCALE IS A RESULT
-----------------------------------------------------------------------------

Every parameter in the reaction and the flow now comes from a measurement or
from geometry: the rate constant and its activation energy from Palandri &
Kharaka (2004) for oligoclase; the solute ceiling and its enthalpy from quartz
saturation; the matrix conductivities from Goodfellow et al. (2016); the joint
conductivity from a 100 um aperture through the cubic law; the diffusivity of
dissolved silica from Rebreanu et al. (2008), scaled by Stokes-Einstein; and
tau and the saturation length from the mineralogy and a 2 mm grain size.

Which makes the weathering timescale a PREDICTION, and it can be checked. The
default settings -- 1 m joints, 0.30 m/yr, 12 C -- take 3713 kyr to dissolve
90 % of a 3 m section, a weathering front of about 0.81 m/Myr. Measured
temperate granite regoliths give 7 m/Myr at Panola and 4 m/Myr at Davis Run
(White et al. 2001); tropical Rio Icacos runs 43-45 m/Myr, which is the right
direction for a warmer, wetter site. So this runs five to nine times SLOWER
than the temperate field rate, with nothing tuned to it.

The gap is worth stating precisely, because it is not evenly distributed
across the inputs. It sits in the reactive surface area, which is here the
GEOMETRIC area of 2 mm grains, 900 m2/m3, while BET for granite is 3e5 to
3e6 -- two to three orders higher. Closing the rate gap needs a factor of
five, which is still five hundred times below BET. So the discrepancy lives
inside a range the field has not resolved (White & Brantley 2003), and it is
reported rather than removed: choosing a surface area to make the rate come
out is the one move that would make this number meaningless.

Two earlier versions got closer and were wrong to. The calibrated model ran at
17.9 m/Myr, three times too FAST, with nothing checking it. Deriving tau
brought it to 2.5, which looked like agreement within a factor of two -- but
30 % of the section was part-dissolved at once, so it was not advancing a
front at all, and dividing 3 m by t90 was not measuring one. Only with the
matrix transport corrected, and the part-dissolved zone down to 4.6 %, is the
number a front rate at all. It is a worse match and a better measurement.

Solute moves by advection **and by diffusion**:

    div(q c) - div(D grad c) = r (1 - c),      r = k A / C_eq

Without the diffusive term, rock off a flow path never weathers at all: a block
interior saturates and then sits at ``c = 1`` for ever, and the model gives
joints entirely dissolved beside blocks entirely untouched with nothing in
between. Diffusive export of solute toward a flushed joint keeps the interior
undersaturated, which is what forms a weathering rind -- and, because a corner
sheds solute to two faces and a face to one, what rounds corestones.

Flow is steady Darcy flow, solved for the hydraulic head,

    div( K grad H ) = 0,        H = psi - d        (d is depth, positive down)

with infiltration prescribed at the surface, a fixed head at the base and
no-flow sides. Conductance lives on the **links** between cells, which is where
the fracture network lives: a fractured link conducts, an intact one barely
does. Lateral flow along a subhorizontal joint is then not a special case -- it
is what the head field does when a low-resistance path exists.

An earlier version routed water down a gravity cascade instead, each cell
handing its water to the three cells below. That could enter a horizontal joint
but never travel along one, so the whole horizontal set was inert: deleting it
changed the weathering by 0.2 percentage points. The cascade was cheaper, but
the thing it could not represent was exactly the thing the joints are for.

The conductance EVOLVES: dissolving rock opens connected porosity, so the
matrix conducts better as it weathers, water is drawn into the weathered zone,
and it dissolves faster still. That feedback is what makes a weathering profile
-- shallow blocks destroyed, deeper ones surviving -- rather than a section
that weathers uniformly at every depth. The head is therefore re-solved as the
rock changes, every ``flow_interval`` steps; see :meth:`link_conductivity`.

(Two claims stood here until 2026-09-03 and both were stale. The conductance
was fixed and the head solved once, which is the simplification this replaces.
The solute was also described as swept row by row with a tridiagonal solve
inside each row -- that went when diffusion arrived and the whole field became
one sparse advection-diffusion-reaction solve.)

**Every parameter here is a placeholder.** None is measured. They are tabulated
in the design document, and no number from this module should be used as a
result.
"""

import inspect

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

#: What SciPy calls the relative tolerance of an iterative solve. It was ``tol``
#: until 1.12, ``rtol`` from 1.12, and ``tol`` was REMOVED in 1.14. Pyodide
#: ships 1.14, so a hardcoded ``tol=`` runs fine against an older SciPy on a
#: workstation and raises TypeError in the browser -- and since the model runs
#: inside a web worker there, the traceback never reaches the page console. The
#: demo simply refused to advance, silently. Resolved once, from the signature.
_RTOL = ("rtol" if "rtol" in inspect.signature(spl.bicgstab).parameters
         else "tol")

#: Seconds in a Julian year. Time is seconds internally; years appear only at
#: the input and output edges, and this is where the conversion is named.
YEAR = 365.25 * 24 * 3600.0

#: Column ordering for every sparse factorisation here. Both matrices -- the
#: conductance Laplacian for the head and the transport operator for the solute
#: -- are STRUCTURALLY SYMMETRIC, because every link contributes both (i, j) and
#: (j, i). That is the precondition for minimum degree on ``A + A.T``, and it is
#: asserted in the tests rather than assumed, since an upwind stencil that
#: reached only one way would silently break it. SuperLU's default is COLAMD,
#: which orders for ``A.T A`` and is the right choice only when the structure is
#: unsymmetric; here it roughly doubles the fill. Measured on the 22,650-cell
#: section: nnz(L + U) 2,070,424 -> 1,060,262, factorisation 302 -> 170 ms,
#: back-substitution 3.16 -> 2.74 ms.
ORDERING = "MMD_AT_PLUS_A"


#: Water density [kg/m3] and gravity [m/s2], for the cubic law.
RHO_WATER = 1000.0
GRAVITY = 9.81


def water_viscosity(T):
    """
    Dynamic viscosity of liquid water [Pa s], for temperature ``T`` in K.

    The standard engineering correlation,

        mu(T) = 2.414e-5 * 10 ** (247.8 / (T - 140))

    checked here against tabulated values: -2.2 % at 273.15 K, -0.0 % at
    293.15 K and -0.2 % at 313.15 K, so within a couple of per cent across
    every temperature this model is used at. It is in the code rather than in
    a table because it is the only place a *physical property of water*
    enters, and a student should be able to see that it is not a fit to
    anything in this model.
    """
    return 2.414e-5 * 10.0 ** (247.8 / (T - 140.0))

#: Molar mass of dioxygen [g/mol]. Present because every water-chemistry table
#: of dissolved oxygen is in mg/L and this model is in mol/m3, which are the
#: same units divided by this number.
M_O2 = 32.0


def oxygen_solubility(T):
    """
    Dissolved oxygen in fresh water in equilibrium with the atmosphere
    [mol/m3], for temperature ``T`` in K at one atmosphere.

    The standard freshwater saturation correlation, which is stated in mg/L:

        ln C = -139.34411 + 1.575701e5 / T - 6.642308e7 / T^2
        + 1.2438e10 / T^3 - 8.621949e11 / T^4

    **Warming LOWERS this**, which is the opposite of what temperature does to
    every other ceiling in this model. A gas is driven out of solution as the
    water warms: 0.457 mol/m3 at 0 C against 0.236 at 30 C, a factor of 1.93
    the wrong way. So where a warmer, more soluble silica fluid carries more
    away per litre, warmer water carries LESS oxygen in, and the two
    temperature effects on the weathering rate point in opposite directions.
    That is a physical fact about oxygen and not a modelling choice.

    Provenance, marked honestly because the two halves are not equally sure.
    The coefficients are the standard ones (Benson & Krause, in the form used
    by APHA Standard Methods 4500-O and the USGS DOTABLES tool) and are
    written here from memory: THE ATTRIBUTION HAS NOT BEEN CHECKED AGAINST THE
    PRIMARY. What has been checked, in the session that added this, is that
    they are not mistranscribed, two independent ways.

    Against the tabulated freshwater saturation values, 0 to 30 C: this
    returns 14.621, 12.771, 11.288, 10.084, 9.092, 8.263 and 7.559 mg/L at 0,
    5, 10, 15, 20, 25 and 30 C, every one within 0.04 % of the table.

    And against Henry's law, with a constant that has nothing to do with this
    correlation. Moist air at one atmosphere and 25 C carries an oxygen
    partial pressure of 0.20291 atm; at a Henry constant of 1.3e-3 mol per
    litre per atmosphere that is 8.44 mg/L, against this correlation's 8.263.
    Two per cent apart, on a constant good to two significant figures.

    A five-term polynomial in 1/T that lands on a table at seven temperatures
    and on Henry's law independently is the right polynomial. Whose it is, is
    the part still to confirm.
    """
    return np.exp(-139.34411 + 1.575701e5 / T - 6.642308e7 / T ** 2
                  + 1.2438e10 / T ** 3 - 8.621949e11 / T ** 4) / M_O2


class Weathering(object):
    """
    Weathering of a jointed rock section, driven by infiltrating rainwater.

    Takes a seeded :class:`~corestone.fractures.FractureNetwork` and evolves
    the fraction of the soluble mineral phase that has dissolved in each cell.
    """

    def __init__(self, network):
        self.network = network            # the seeded joint network
        self.dx = network.dx              # cell size [m]
        self.nz = network.nz              # rows, increasing downward
        self.nx = network.nx              # columns, increasing rightward

        # ---- parameters. ALL PLACEHOLDERS; see design/02-teaching-scope.md
        self.infiltration = 0.30 / YEAR   # recharge at the surface [m/s]
        # A JOINT IS A GEOMETRY, NOT A CONDUCTIVITY. What is set here is the
        # aperture; the conductivity follows from it by the cubic law, and so
        # does its dependence on cell size. See :attr:`k_fracture`.
        #
        # 100 um is the HYDRAULIC aperture -- the one the cubic law wants --
        # and not the millimetre-plus opening you can see at an outcrop. The
        # walls touch at asperities, the surfaces are rough, and near-surface
        # joints carry clay and grus, so the hydraulic aperture runs one to
        # two orders below the mechanical one. Bracketed from both sides:
        #
        #   Rukavickova et al. (2021) measured fractured granitoid at
        #   1e-8 to 1e-7 m/s, borehole scale; inverted through the cubic law
        #   for one joint set at 0.5-2 m spacing that is 20-67 um. Those
        #   boreholes are deep, where joints are held shut.
        #
        #   Laboratory tension fractures in granite run from 250 um down to
        #   4 um as normal stress rises to 20 MPa, so the unstressed end is
        #   100-250 um.
        #
        # This model is the top three metres, where the normal stress is
        # essentially nil, so it belongs at the open end of the deep range and
        # the low end of the unstressed one. 100 um is both.
        #
        # PRE-WEATHERING, and it never changes. Spheroidal weathering does not
        # widen a joint into a void: the wall rock rots to grus and the grus
        # stays where it is. What opens is the MATRIX beside the joint, which
        # is k(M) and is modelled -- and that alone drops the joint's
        # advantage from 20000x to 1.8x over a run. Aperture growth belongs to
        # a later stage, once grus has been flushed out and there is a cavity;
        # that is the channelization this model excludes.
        self.joint_aperture = 100.0e-6    # hydraulic aperture of a joint [m]
        # The two ends of the matrix, and the only two numbers in this model
        # taken from a measurement rather than invented. Goodfellow et al.
        # (2016), JGR Earth Surface 121, 1410-1435, measured the hydraulic
        # conductivity of granodiorite MATRIX across a range of weathering
        # grades: 9e-9 to 8e-8 cm/s in the parent rock and 9e-5 to 9e-4 cm/s
        # in the most weathered samples, an increase of three to four orders
        # of magnitude. These are the mid-points of those two ranges, in m/s.
        # Called k_weathered rather than k_grus because that is what the
        # measurement is -- the most weathered samples in a granodiorite
        # suite. Grus is a particular material with a particular fabric, and
        # naming a conductivity after it claims more than the number carries.
        # AT T_K_ref. A hydraulic conductivity is k_intrinsic rho g / mu, so
        # a measured one belongs to the temperature it was measured at;
        # Goodfellow et al.'s permeameter work is laboratory temperature.
        # Corrected to the working temperature by viscosity_factor, exactly
        # as the joint is through the cubic law -- if only one of them
        # carried the viscosity the joint-to-matrix contrast would depend on
        # temperature, which is an artefact and not a physics.
        self.T_K_ref = 293.15             # temperature of the measured
                                          # conductivities [K]
        self.k_matrix = 5.0e-10           # intact granite [m/s]
        self.k_weathered = 5.0e-6         # fully weathered matrix [m/s]
        self.flow_tolerance = 0.01        # re-solve the head once the rock has
                                          # changed this much anywhere. NOT a
                                          # step count: that would tie the
                                          # answer to the step size, which the
                                          # drift control exists to prevent.
                                          # CONVERGED, not chosen: 0.2, 0.1,
                                          # 0.05, 0.02, 0.01 give max|dM| of
                                          # 1.9e-1, 9.4e-2, 4.4e-2, 9.4e-3 and
                                          # 0 against 0.002. A tolerance is a
                                          # numerical quantity and should be
                                          # converged; the uncertainty belongs
                                          # in the conductivities.
        # DERIVED. L = q C_eq / (k A), and every term is sourced:
        #   q    0.30 m/yr infiltration = 9.51e-9 m/s
        #   C_eq 0.10 mol Si/m3, quartz saturation
        #   k    10^-11.84 mol m-2 s-1, oligoclase neutral (Palandri &
        #        Kharaka 2004), times 1.6 mol Si per formula unit
        #   A    reactive surface area, 900 m2/m3 -- GEOMETRIC, for 2 mm
        #        cubic grains at 30 % plagioclase, which is an ordinary
        #        coarse granite
        # giving L_ref = 0.457 m, where the old calibration guessed 0.50.
        #
        # A is the one real choice and it is the deep one. BET surface area
        # of granite is 0.1-1 m2/g, i.e. 3e5-3e6 m2/m3, two to three orders
        # ABOVE the geometric figure, and using it would give a saturation
        # length of about a millimetre. That gap is the long-standing
        # laboratory-versus-field rate discrepancy (White & Brantley 2003):
        # laboratory rate constants on BET areas overpredict field weathering
        # by two to three orders of magnitude. Pairing a laboratory k with a
        # GEOMETRIC area is the standard way to land near field behaviour,
        # and the validation in the module docstring is the check that it
        # does.
        self.L_ref = 0.457                # saturation length at T_ref, mean
                                          # infiltration, fresh rock [m]
        self.T_ref = 285.0                # reference temperature [K]
        # The soluble phase is PLAGIOCLASE -- oligoclase, the An10-30 the
        # feldspar of a granite usually is. Not because it is the most
        # abundant phase (by the IUGS definition granite is 10-65 % of its
        # feldspar as plagioclase, so K-feldspar often wins on volume) but
        # because it sets the rate: at 25 C and neutral pH it dissolves 3.7x
        # faster than K-feldspar, and given equal surface areas it still
        # carries 79 % of the dissolution. Taking the resistant phase's
        # kinetics would describe a rock paced by what survives, and what
        # makes grus is the phase that goes.
        self.E_a = 69.8e3                 # activation energy [J/mol]:
                                          # oligoclase, neutral mechanism,
                                          # Palandri & Kharaka (2004) Table 13
                                          # (K-feldspar, Table 15, is 38.0e3)
        # Enthalpy of the equilibrium that stops the reaction [J/mol], van 't
        # Hoff. NOT a property of the dissolving mineral: it is a statement
        # about WHAT saturates. Quartz is treated as inert here and silica
        # saturates readily, so the ceiling on the solute is silica, and this
        # is quartz's dissolution enthalpy.
        #
        # The alternative reading -- solution buffered by kaolinite -- gives a
        # NEGATIVE value for calcic plagioclase, computed from llnl.dat as
        # +24.0 kJ/mol for albite to kaolinite and -151.3 for anorthite,
        # so -11 at An20 and -46 at An40. That flips the sign of the
        # temperature effect on C_eq: at An40 the section takes 124 kyr to
        # reach 90 % dissolved at 0 C and 308 kyr at 30 C, i.e. warming the
        # rock slows it down. A real regime, and not this model's.
        self.delta_H_r = 32.9e3           # quartz, llnl.dat at 25 C
        # Both matter twice over. (E_a - delta_H_r) alone sets how the
        # saturation length moves with temperature, and C_eq enters tau
        # separately, so the pair has to come from one consistent story
        # rather than being picked one at a time.
        self.R_gas = 8.314                # gas constant [J/mol/K]
        # DERIVED, not calibrated. tau = M0 / C_eq, and both come from
        # elsewhere in this file:
        #
        #   M0: oligoclase An20 has M = 265.4 g/mol at 2640 kg/m3, so a molar
        #       volume of 100.5 cm3/mol. At 30 % of the rock by volume that is
        #       2984 mol/m3 of plagioclase. Dissolving incongruently to
        #       kaolinite releases 2 Si per albite and 0 per anorthite -- the
        #       anorthite silicon stays in the clay -- so 1.6 Si per formula
        #       unit at An20, giving M0 = 4774 mol Si per m3 of rock.
        #
        #   C_eq: quartz saturation, 1.0e-4 mol/kg, i.e. 0.10 mol Si/m3. The
        #       same ceiling that sets delta_H_r, so the two agree.
        #
        # tau = 4774 / 0.10 = 47744. The placeholder was 6700, seven times
        # too few volumes, and the model therefore weathered seven times too
        # fast. See the note on validation in the module docstring.
        self.tau_ref = 47744.0            # M0/C_eq at T_ref: volumes of
                                          # saturated water per volume of rock
        # MEASURED, and of the right species. The solute here is silica --
        # C_eq is quartz saturation -- and the diffusion coefficient of
        # dissolved silica is (1.02 +/- 0.02)e-9 m2/s at 25 C (Rebreanu,
        # Vanderborght & Chou, 2008, Marine Chemistry 112, 230-233),
        # confirming Wollast & Garrels (1971) at (1.00 +/- 0.05)e-9. This was
        # carried for a long time as "order of magnitude for a dissolved ion",
        # which undersold it: it is the measurement, for the ion in question.
        self.D_molecular = 1.0e-9         # aqueous diffusivity of dissolved
        self.T_D_ref = 298.15             # silica [m2/s] AT this temperature
                                          # [K]; scaled to the working
                                          # temperature by Stokes-Einstein,
                                          # see diffusivity_factor.
        # Tortuosity divides the aqueous diffusivity to give the matrix
        # value. 10 puts D_matrix at 1e-10 m2/s, which is right for WEATHERED
        # material: saprolite at porosity ~0.3 and tortuosity ~3 gives
        # D_eff/D_0 ~ 0.1. It is emphatically NOT right for intact granite,
        # where through-diffusion experiments give D_eff/D_0 of 1e-4 to 1e-5,
        # a tortuosity nearer 1e4.
        #
        # KNOWN INCONSISTENCY, stated rather than hidden: the conductivity
        # evolves with M and this does not. Dissolving rock opens porosity to
        # diffusion exactly as it opens it to flow, so a full treatment would
        # interpolate tortuosity between the two ends the way
        # link_conductivity interpolates k. Fresh rock therefore diffuses far
        # too freely here, which flatters the early rind.
        # Tortuosity is the cost of the detour. A molecule diffusing through
        # rock cannot travel in a straight line: it follows the pore network
        # around every grain, so it covers far more distance than the
        # separation it achieves. Diffusivity in rock is the free-water value
        # divided by this. The two ends are the two rocks.
        self.tortuosity_fresh = 1.0e4     # intact crystalline rock, which has
                                          # almost no connected porosity:
                                          # measured 2e-14 to 1.3e-12 m2/s
                                          # against a free-water 1e-9, so
                                          # 1e3 to 1e5, centre 1e4
        self.tortuosity_weathered = 10.0  # saprolite, which is full of holes:
                                          # near 30 % porosity, D_eff/D_0 ~ 0.1

        # Longitudinal dispersivity scales with the transport distance:
        # Gelhar, Welty & Rehfeldt (1992) put it near a tenth of the scale
        # observed, with orders of magnitude of scatter and their
        # high-reliability data at the low end. Over a 3 m section a tenth
        # would be 0.3 m; 0.05 m is a fiftieth, inside the scatter and
        # deliberately conservative, since dispersion this large would smear
        # the rind the model exists to show.
        # Grain size, which does two jobs and must do them consistently. It
        # sets the reactive surface area behind L_ref above (900 m2/m3 for
        # 2 mm cubes at 30 % plagioclase) and it sets the dispersivity below.
        self.grain_size = 2.0e-3          # mean grain diameter [m]
        # Two ARBITRARY cut-offs on a continuous field, kept for convenience
        # and named badly. Neither word is a fraction dissolved: fresh rock,
        # saprock, saprolite and grus are distinguished by fabric and
        # mineralogy, and a corestone is a SHAPE -- a rounded block surrounded
        # by weathered rock -- so intact bedrock at depth satisfies x_core
        # without being a corestone at all. They were reported as percentages
        # in the demo, which made two claims this model cannot make; that has
        # been removed. Set them yourself if you want them to mean something.
        self.x_grus = 0.50                # ARBITRARY: "mostly dissolved"
        self.x_core = 0.05                # ARBITRARY: "barely touched"
        # No f_inert. Quartz used to be carried here as a 0.30 fraction that
        # "never dissolves" -- a knob, set by hand, and read by nothing in the
        # code. It is not a knob. C_eq is quartz saturation, ~6 mg/L as SiO2
        # at 25 C (llnl.dat: log K -3.999; chalcedony would be 11.2 and
        # amorphous silica 116.2), so quartz sits exactly AT the ceiling this
        # model normalises to and its driving force (1 - C/C_eq) is zero by
        # construction. Quartz is inert here as a CONSEQUENCE of the
        # saturation choice, not as an assumption laid on top of it, and the
        # same silica ceiling is what stops the feldspar.
        self.c_drift_max = 0.03           # THE step control: how far c may move
                                          # while it is held across a step.
                                          # This is the model's one time-step
                                          # approximation, so it is the thing
                                          # bounded. Error is very nearly
                                          # linear in it -- halve it, halve the
                                          # error. A CHOSEN error budget; see
                                          # update() for what each value costs.
        self.dt_max = 50000.0 * YEAR      # ceilings, not the control. They
        self.dx_max = 0.50                # bound a step that the drift control
                                          # has no history for (the first) or
                                          # would let run away (a field that
                                          # has stopped moving).
        self.dt_growth = 2.0              # most a step may lengthen at once
        self.dt_min = 1.0 * YEAR          # floor, so a pathological drift
                                          # cannot spin the controller for ever
        self.krylov_tol = 1.0e-10         # convergence of the reused solve
        self.max_krylov_iterations = 15   # past this, refactorise instead

        # ---- state
        self.T = self.T_ref               # temperature [K]
        self.t = 0.0                      # model time [s]
        self._tort = None                 # link tortuosity, refreshed with
                                          # the head; see solve_flow
        self.M = None                     # soluble mineral remaining, M/M0
        self.c = None                     # normalised concentration C/C_eq
                                          # LEAVING each cell, not entering
        self.H = None                     # hydraulic head [m]
        self.q = None                     # through-flux per cell [m2/s]
        self.q_v = None                   # flux on vertical links, down [m2/s]
        self.q_out_base = None            # flux leaving through the base
        self.q_h = None                   # flux on horizontal links, right
        self._in_above = None             # inflow from the row above
        self._in_left = None              # inflow from the left neighbour
        self._in_right = None             # inflow from the right neighbour
        self._T = None                    # constant part of the solute operator
        self._T_key = None                # the values it was built from
        self._A = None                    # the step matrix, reused in place
        self._diag = None                 # where its diagonal sits in .data
        self._flow_lu = None              # cached head factorisation, kept
        self._H_prev = None               # across solves as a preconditioner
        self.head_tol = 1e-12             # CG residual tolerance on the head
        self.max_head_iterations = 20     # above this, refactorise instead
        self._lu = None                   # cached factorisation, reused
        self._x = None                    # last solute solution, unclipped
        self._dt = None                   # the step the drift control chose
        self._c_held = None               # the c actually held over a step
        self._drift = None                # the drift the last step produced
        self._M_flow = None               # M when the head was last solved
        self.flow_solves = 0              # how many times the head was solved
        self.rejected_steps = 0           # steps retried for overrunning
        self.factorisations = 0           # how many times it was rebuilt

    # ---- parameter setters (one per parameter; units in the docstring)

    def set_temperature(self, value):
        """Temperature [K]. Higher shrinks the equilibration length."""
        self.T = value

    def set_infiltration(self, value):
        """Recharge at the ground surface [m/s]."""
        self.infiltration = value

    def set_saturation_length(self, value):
        """Saturation length at the reference temperature [m]."""
        self.L_ref = value

    # ---- the physics

    @property
    def rate_factor(self):
        """``k(T) / k_ref``, Arrhenius."""
        return np.exp(-(self.E_a / self.R_gas)
                      * (1.0 / self.T - 1.0 / self.T_ref))

    @property
    def solubility_factor(self):
        """
        ``C_eq(T) / C_eq_ref``, van 't Hoff.

        Solubility is temperature dependent too, and treating it as constant
        was not a small error: in the saturation-limited regime -- which is most
        of this model -- the amount dissolved scales with ``C_eq`` and not with
        the rate constant at all, so holding it fixed removed the half of the
        temperature dependence that dominates.
        """
        return np.exp(-(self.delta_H_r / self.R_gas)
                      * (1.0 / self.T - 1.0 / self.T_ref))

    @property
    def saturation_length(self):
        """
        The e-folding length of the approach to saturation [m], for water at
        the mean infiltration rate through fresh rock.

            saturation_length = q * C_eq / (k(T) * A)

        **Not** a distance at which equilibrium is reached. There is no
        equilibrium: ``c`` approaches 1 asymptotically, and after n of these
        lengths the undersaturation is ``exp(-n)`` of what it was. The system
        has an asymptote, and the model has a normalisation rather than a
        thermodynamics -- ``C_eq`` never appears alone.

        Because it goes as ``C_eq / k``, its temperature dependence is set by
        ``(E_a - delta_H_r)``, not by ``E_a`` alone.
        """
        return self.L_ref * self.solubility_factor / self.rate_factor

    # ---- the thermodynamics, named
    #
    # These four exist to be READ. None is needed by the solver: every one is
    # a quantity the model already computes implicitly, given a name and a
    # docstring so that a student can ask the model what regime it is in and
    # get an answer instead of inferring it from behaviour.

    @property
    def k_fracture(self):
        """
        Hydraulic conductivity of a jointed link [m/s], from the cubic law.

            k_fracture = rho g a^3 / (12 mu dx)

        The "cubic law" is not a law of its own. It is the Navier-Stokes
        solution for steady laminar flow between two parallel plates -- plane
        Poiseuille flow -- integrated across the gap, which puts discharge in
        proportion to the cube of the separation. The name and its
        applicability to real, rough, deformable rock fractures come from
        Witherspoon, Wang, Iwai & Gale (1980), Water Resources Research 16,
        1016-1024, who showed it holds down to apertures of a few microns
        provided the aperture is the HYDRAULIC one rather than the visible gap.

        A planar joint of aperture ``a`` carries transmissivity
        ``rho g a^3 / 12 mu``; smearing it across a cell of width ``dx`` gives
        a conductivity. Derived rather than set, for two reasons.

        It makes the joint a measurable object. An aperture can be measured
        and is in the literature; the conductivity of a joint smeared over an
        arbitrary cell is a modelling artefact and is in nobody's table.

        And it makes the cell size mean what the exercise says it means. Held
        as a CONSTANT conductivity, the implied aperture moved with the grid
        -- 91 um at 5 cm against 67 um at 2 cm -- so choosing a finer grid
        quietly tightened the joints by a third, while the page promised that
        cell size is the numerical grid and not the rock. Transmissivity is
        the invariant, so the conductivity has to scale as 1/dx, and now does.

        Temperature enters through the viscosity, and the MATRIX ends carry
        it too -- see :attr:`k_matrix_at_T`. That matters more than it looks.
        Hydraulic conductivity is ``k_intrinsic rho g / mu`` for any medium,
        so warming raises joints and matrix alike and leaves their ratio
        alone; with the infiltration prescribed at the surface rather than
        driven by a head gradient, an unchanged ratio means an unchanged flow
        field. Letting only the joint carry the viscosity, which this model
        did briefly, doubled the joint-to-matrix contrast between 0 and 30 C
        and moved the speed field by 55 % -- a temperature effect on the flow
        that has no physical basis and was purely an asymmetry in the code.
        """
        mu = water_viscosity(float(np.mean(self.T)))
        return (RHO_WATER * GRAVITY * self.joint_aperture ** 3
                / (12.0 * mu * self.network.dx))

    @property
    def dispersivity(self):
        """
        Longitudinal dispersivity [m]: one grain diameter.

        Water threading a porous medium does not travel at one speed. It takes
        many paths, of different lengths, through pores of different widths,
        and a solute front smears out because of it. The dispersivity is the
        distance over which those paths differ from one another -- so at the
        pore scale it is the grain size, and it multiplies the local velocity
        to give a spreading coefficient with the units of a diffusivity.

        The scale matters as much as the number. A dispersivity measured in
        the field is much larger, because it is standing in for heterogeneity
        the modeller did not draw; here the heterogeneity IS drawn, cell by
        cell, as the joint network. Design 07 has the arithmetic.
        """
        return self.grain_size

    def link_tortuosity(self):
        """
        Matrix tortuosity on each link: how far a diffusing molecule must
        detour around grains, and therefore how much slower diffusion is here
        than in open water.

            tortuosity(M) = tortuosity_fresh^M * tortuosity_weathered^(1 - M)

        It follows the rock because dissolution is what opens the detour up.
        Intact granite has almost no connected pore network and a molecule can
        barely cross it; as the soluble phase goes, the holes it leaves join up
        and the path straightens out, by three orders of magnitude between the
        two ends. Geometric interpolation on the mean ``M`` of the two cells a
        link joins -- the same form as :meth:`link_conductivity`, because it is
        the same porosity opening to diffusion that opens to flow.

        Returns ``(vertical, horizontal)``. Joint links are handled by the
        caller and are never tortuous: an open aperture is not a maze.
        """
        M = np.clip(self.M, 0.0, 1.0)
        lo, hi = np.log(self.tortuosity_fresh), np.log(self.tortuosity_weathered)

        def tort(m):
            return np.exp(m * lo + (1.0 - m) * hi)

        return (tort(0.5 * (M[:-1, :] + M[1:, :])),
                tort(0.5 * (M[:, :-1] + M[:, 1:])))

    def link_tortuosity_wrap(self):
        """Matrix tortuosity on the periodic seam links; see
        :meth:`link_tortuosity`."""
        M = np.clip(self.M, 0.0, 1.0)
        lo, hi = np.log(self.tortuosity_fresh), np.log(self.tortuosity_weathered)
        m = 0.5 * (M[:, -1] + M[:, 0])
        return np.exp(m * lo + (1.0 - m) * hi)

    @property
    def viscosity_factor(self):
        """
        ``mu(T_K_ref) / mu(T)``: what warming does to any hydraulic
        conductivity.

        ``K = k_intrinsic rho g / mu``, so a conductivity measured at one
        temperature applies at another only after this correction. Water is
        2.2 times as viscous at 0 C as at 30 C, so it is not a small one.
        """
        return (water_viscosity(self.T_K_ref)
                / water_viscosity(float(np.mean(self.T))))

    @property
    def k_matrix_at_T(self):
        """Intact-matrix conductivity at the working temperature [m/s]."""
        return self.k_matrix * self.viscosity_factor

    @property
    def k_weathered_at_T(self):
        """Fully weathered matrix conductivity at the working temperature."""
        return self.k_weathered * self.viscosity_factor

    @property
    def diffusivity_factor(self):
        """
        ``D(T) / D(T_D_ref)``, Stokes-Einstein.

        A diffusing ion is dragged by the water around it, so its diffusivity
        goes as ``T / mu(T)`` -- and over the range of a terrestrial climate
        it is the viscosity that moves, not the absolute temperature. Water is
        2.6 times as viscous at 275 K as at 315 K, so ``D`` spans a factor of
        three across the temperatures this demo offers.

        Left constant, this was the one temperature-dependent quantity in the
        model that was not allowed to depend on temperature. It matters where
        it would be easiest to assume it does not: in the JOINTS mechanical
        dispersion is about 95 times molecular diffusion and this is
        irrelevant, but in the MATRIX the ratio is 0.005 -- molecular
        diffusion beats dispersion two hundred to one -- and the matrix is
        where the weathering rind forms and where corners round. This is the
        term that carries solute out of a block interior.

        Note what it does NOT touch. Viscosity also enters the hydraulic
        conductivity, so warm water should flow more freely -- but the
        infiltration here is prescribed at the surface rather than driven by a
        head gradient, so scaling every conductivity together rescales the
        head and leaves the flow field identical (verified: a 1.6x change
        moves the speed field by 7e-8). The effect is not missing from this
        model; it is unable to act, and that is a property of the boundary
        condition.
        """
        T = np.mean(self.T)
        return ((T / water_viscosity(T))
                / (self.T_D_ref / water_viscosity(self.T_D_ref)))

    @property
    def D_aqueous(self):
        """Molecular diffusivity at the working temperature [m2/s]."""
        return self.D_molecular * self.diffusivity_factor

    @property
    def apparent_activation_energy(self):
        """
        ``E_a - delta_H_r`` [J/mol]: what temperature actually does here.

        The single most misread thing in this model. Warming does two opposing
        things at once -- it speeds the reaction (Arrhenius, ``E_a``) and it
        raises the solubility (van 't Hoff, ``delta_H_r``), and the second
        does not make the rock dissolve faster in place, it lets each litre
        of water carry more away before it stops. The saturation length goes
        as ``C_eq / k``, so those two enter it with OPPOSITE signs and only
        the difference survives.

        That difference is the apparent activation energy of the weathering
        *length scale*, and it is what a field study measuring weathering
        against temperature would recover -- not ``E_a``. With the values
        here, 69.8 - 32.9 = 36.9 kJ/mol, so the length scale is about half as
        temperature-sensitive as the rate constant alone would suggest.

        It can be zero, or negative. If ``delta_H_r`` exceeded ``E_a`` --
        which happens if the ceiling on the solute is set by a reaction whose
        enthalpy is large, and it is what a kaolinite-buffered reading of
        calcic plagioclase gives -- then warming would LENGTHEN the
        saturation length and slow the weathering down. Nothing in the model
        forbids it, and the sign is a consequence of the chemistry, not an
        assumption.
        """
        return self.E_a - self.delta_H_r

    @property
    def damkohler(self):
        """
        Section depth divided by the saturation length [-].

        The dimensionless group that decides which of the two limits this
        model is in, and therefore what the pictures mean. It counts the
        e-foldings of the approach to saturation that a parcel of water
        undergoes on its way down through the section, at the mean
        infiltration and through fresh rock.

        ``Da >> 1`` -- SATURATION-LIMITED. The water saturates long before it
        runs out of rock, so what limits weathering is how much solute each
        litre can carry away, and dissolution happens where fresh water
        arrives. This is the regime that makes corestones, and it is where
        this model sits by construction: about 6 at the reference state, so
        water leaving the base is within ``exp(-6)``, a quarter of a percent,
        of saturation.

        ``Da << 1`` -- REACTION-LIMITED. Water crosses the whole section
        barely touched, weathering is set by the rate constant everywhere at
        once, and the section dissolves uniformly. No corestones: there is
        nothing to shelter a block interior from water that is everywhere
        undersaturated.

        A NOTE ON THE NAME. Physical chemistry calls the ``Da >> 1`` limit
        *transport-limited*, meaning the transport of solute. That name is
        not usable here: in geomorphology transport-limited means a landscape
        whose erosion is set by the capacity to move sediment, and the
        companion exercise on the same teaching site is about exactly that.
        Nothing in this model transports sediment. What runs out is the
        water's capacity to hold solute, so it is called saturation-limited
        here, and a reader meeting *transport-limited* in a chemistry text
        should know it is the same limit.

        The number is the SECTION-scale one. Two others matter and are worth
        forming by hand: the joints carry a higher flux, so their local
        saturation length is longer and they flush deeper than this suggests;
        and the ratio of the JOINT SPACING to the saturation length is what
        decides whether a block interior can be sheltered at all.
        """
        return self.network.nz * self.network.dx / self.saturation_length

    @property
    def regime(self):
        """Name of the limit the model is currently in; see :attr:`damkohler`."""
        da = self.damkohler
        if da > 3.0:
            return "saturation-limited"
        if da < 1.0 / 3.0:
            return "reaction-limited"
        return "mixed"

    def thermo_report(self):
        """
        A readable statement of the thermodynamic state. Returns the text.

        Written to be pasted into a lab report or read aloud in class: every
        number that governs the temperature behaviour, with the arithmetic
        that produced it visible rather than asserted.
        """
        T = float(np.mean(self.T))
        lines = [
            "thermodynamic state",
            "  T                       %8.2f K   (%.2f C)" % (T, T - 273.15),
            "  T_ref                   %8.2f K   (%.2f C)"
            % (self.T_ref, self.T_ref - 273.15),
            "",
            "  E_a                     %8.1f kJ/mol  Arrhenius, on the rate"
            % (self.E_a / 1e3),
            "  delta_H_r               %8.1f kJ/mol  van 't Hoff, on the ceiling"
            % (self.delta_H_r / 1e3),
            "  E_a - delta_H_r         %8.1f kJ/mol  what the length scale feels"
            % (self.apparent_activation_energy / 1e3),
            "",
            "  k(T) / k(T_ref)         %8.3f       reaction this much faster"
            % float(np.mean(self.rate_factor)),
            "  C_eq(T) / C_eq(T_ref)   %8.3f       each litre carries this much more"
            % float(np.mean(self.solubility_factor)),
            "  saturation length       %8.3f m     = L_ref * C_eq-factor / k-factor"
            % float(np.mean(self.saturation_length)),
            "",
            "  Damkohler (section)     %8.2f       depth / saturation length"
            % float(np.mean(self.damkohler)),
            "  regime                  %8s" % self.regime,
        ]
        return "\n".join(lines)

    @property
    def specific_reaction_coefficient(self):
        """
        ``r / M`` [1/s]: the reaction coefficient per unit mineral remaining.

        A SCALAR, and that is the point of naming it. ``r`` falls with the
        soluble mineral because the reactive surface area does, and it falls
        exactly in proportion -- so the quotient carries no ``M`` at all, and
        the ratio can be formed without ever dividing by a mineral content that
        is on its way to zero. It is what lets the rock be integrated exactly
        rather than tangentially; see :meth:`update`.
        """
        return ((self.infiltration / self.L_ref)
                * self.rate_factor / self.solubility_factor)

    @property
    def reaction_coefficient(self):
        """
        ``r = k A / C_eq`` [1/s]: the rate at which undersaturation is consumed.

        This is the flux-independent form. Since ``r = q / saturation_length``
        and the saturation length is itself proportional to ``q``, the flux
        cancels **explicitly** here rather than implicitly inside an exponent.
        It falls with the soluble mineral remaining, because the reactive
        surface area does.
        """
        return self.specific_reaction_coefficient * np.maximum(self.M, 0.0)

    @property
    def tau(self):
        """
        ``M0 / C_eq``: volumes of saturated water needed per volume of rock.

        Falls as solubility rises, so a warmer and more soluble fluid carries
        more away per unit volume. This is the second place ``C_eq`` enters.
        """
        return self.tau_ref / self.solubility_factor

    def local_saturation_length(self):
        """The saturation length cell by cell [m]: ``q / r``."""
        r = self.reaction_coefficient
        return (self.q / self.dx) / np.maximum(r, 1e-300)

    def transport_coefficients(self):
        """
        Solute transport coefficient on every link [m2/s].

            D = D_aqueous(T) / tortuosity(M)  +  grain_size * |v|

        Molecular diffusion, at the working temperature, plus hydrodynamic
        (mechanical) dispersion. The
        second is sometimes called turbulent diffusion and the operator is the
        same, but at these fluxes the Reynolds number is around 3e-5 -- laminar
        by five orders of magnitude. The distinction matters only because it is
        pore and aperture geometry that sets the dispersivity, not eddies.

        Without this term the model is pure advection, and rock that is not on
        a flow path never weathers at all: a block interior saturates and then
        sits at ``c = 1`` for ever. Diffusive export of solute toward a flushed
        joint is what keeps the interior undersaturated, and is therefore what
        lets a weathering rind form. It is also the geometric route to
        spheroidal rounding, since a corner sheds solute to two faces and an
        edge to one.
        """
        D = self.D_aqueous
        tv, th = self._tort if self._tort is not None else self.link_tortuosity()
        dm_v = np.where(self.network.link_v, D, D / tv)
        dm_h = np.where(self.network.link_h, D, D / th)
        return (dm_v + self.dispersivity * np.abs(self.q_v) / self.dx,
                dm_h + self.dispersivity * np.abs(self.q_h) / self.dx)

    def _transport_operator(self):
        """
        The part of the solute operator that never changes: advection,
        diffusion and the outflow through the base.

        Built once. The flow field is static, so only the reaction term on the
        diagonal moves from step to step -- which is what makes a cached
        factorisation worth keeping.
        """
        # Cache on the values that built it. A bare `is not None` check made
        # changing D_molecular or the dispersivity silently do nothing, which a
        # test caught -- the same shape of defect as a docstring drifting from
        # its code, and one a cache invites.
        tv, th = self._tort if self._tort is not None else self.link_tortuosity()
        key = (self.D_aqueous, self.dispersivity,
               float(np.sum(tv)), float(np.sum(th)),
               id(self.q_v), float(np.sum(self.q_v)), float(np.sum(self.q_h)))
        if self._T is not None and self._T_key == key:
            return self._T
        self._lu = None                    # the factorisation went with it
        nz, nx, dx = self.nz, self.nx, self.dx
        idx = np.arange(nz * nx).reshape(nz, nx)
        D_v, D_h = self.transport_coefficients()
        rows, cols, vals = [], [], []

        def add(i, j, v):
            rows.append(np.asarray(i).ravel())
            cols.append(np.asarray(j).ravel())
            vals.append(np.asarray(v, dtype=float).ravel()
                        * np.ones(np.asarray(i).size))

        for a_, b_, f in ((idx[:-1, :], idx[1:, :], self.q_v),
                          (idx[:, :-1], idx[:, 1:], self.q_h)):
            fwd, rev = np.maximum(f, 0.0), np.maximum(-f, 0.0)
            add(a_, a_, fwd);  add(b_, a_, -fwd)
            add(b_, b_, rev);  add(a_, b_, -rev)
        for a_, b_, D in ((idx[:-1, :], idx[1:, :], D_v),
                          (idx[:, :-1], idx[:, 1:], D_h)):
            add(a_, a_, D);  add(a_, b_, -D)
            add(b_, b_, D);  add(b_, a_, -D)
        if self.network.periodic_x:
            fwd = np.maximum(self.q_wrap, 0.0)
            rev = np.maximum(-self.q_wrap, 0.0)
            Dw = np.where(self.network.link_wrap, self.D_aqueous,
                          self.D_aqueous / self.link_tortuosity_wrap()) \
                 + self.dispersivity * np.abs(self.q_wrap) / dx
            l, rgt = idx[:, -1], idx[:, 0]
            add(l, l, fwd + Dw);     add(rgt, l, -(fwd + Dw))
            add(rgt, rgt, rev + Dw); add(l, rgt, -(rev + Dw))
        add(idx[-1, :], idx[-1, :], self.q_out_base)

        self._T = sp.coo_matrix(
            (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
            shape=(nz * nx, nz * nx)).tocsc()
        self._T_key = key

        # Where the diagonal sits inside the CSC value array, so the reaction
        # term can be written straight into it every step instead of building a
        # second matrix and adding. Every cell has a diagonal entry -- every
        # cell has at least one link, and a link puts a coefficient on both of
        # its cells' diagonals -- and that is asserted here rather than trusted,
        # because a missing one would silently drop that cell's reaction.
        n = nz * nx
        col_of = np.repeat(np.arange(n), np.diff(self._T.indptr))
        self._diag = np.nonzero(self._T.indices == col_of)[0]
        if self._diag.size != n:
            raise AssertionError("%d of %d cells have no diagonal entry"
                                 % (n - self._diag.size, n))
        self._A = self._T.copy()
        return self._T

    def _step_matrix(self, r):
        """
        The operator for one step: transport, plus the reaction on the diagonal.

        Only the diagonal moves from step to step, so this writes into a matrix
        that is allocated once rather than adding two sparse matrices and
        converting the result. The values are identical -- the same two floats
        added in the same order -- and it costs 0.105 ms in place against
        1.025 ms via ``sp.diags`` at 22,650 cells.
        """
        T = self._transport_operator()
        A = self._A
        A.data[:] = T.data
        A.data[self._diag] += np.broadcast_to(
            r * self.dx * self.dx, (self.nz, self.nx)).ravel()
        return A

    def solve_solute(self, r):
        """
        Steady advection-diffusion-reaction for the normalised concentration.

            sum_out f c_i - sum_in f c_j + sum_links D (c_i - c_j)
                + r dx^2 c_i  =  r dx^2

        One sparse solve for the whole field. Diffusion is not directional, so
        the row-by-row sweep that pure advection allowed is gone -- and with it
        the cyclic Sherman-Morrison seam solve and the ``Q * beta * (1 - c)``
        pickup, which cancelled catastrophically once ``beta`` grew past about
        1e5. The reaction term here is linear and bounded, so neither hazard
        remains.

        The solve STARTS FROM THE PREVIOUS STEP'S FIELD. The rock changes by at
        most ``dx_max`` of its mineral content in a step, so the concentration
        moves very little, and the old field is a far better guess than the
        zero the solver would otherwise assume. This costs no accuracy: the
        convergence test is on the residual ``||b - A x||``, which does not
        know where the iteration started.

        Water entering at the surface carries ``c = 0``, so it contributes
        nothing to the inflow sum: undersaturation enters only through the
        reaction term on the right.

        Cost. Only the diagonal changes between steps, so the LU factorisation
        is cached and reused as a preconditioner: at 22,650 cells a
        back-substitution is 2.7 ms against 170 ms to refactorise.

        WHY IT IS REFRESHED, corrected. This docstring used to say the
        preconditioner is rebuilt when the iteration count says so, and called
        that self-tuning. It is not what happens. Raising
        ``max_krylov_iterations`` from 15 to a million changes nothing at all:
        the same 15 factorisations, the same answer to the last bit. Every
        refresh but the first is triggered by ``info = -10`` -- a BiCGSTAB
        BREAKDOWN, after two iterations, on 14 of 108 steps at dx = 0.02. The
        algorithm divides by a quantity that has gone to zero, which happens
        here because the preconditioner is very good rather than because it is
        stale. The cap is a backstop that has never fired at its current value.

        Lowering the cap so that it DOES fire is worth about 15 %, and that is
        an untaken decision rather than an oversight: it is a tuned number,
        the optimum moves between 4 and 5 from case to case, and it buys
        little. Measured, min of three alternating runs, against cap = 15:

            cap     app 3 m   fine 3 m   45 deg   wide
              4       x1.28      x0.86    x1.09   x1.38
              5       x1.28      x1.14    x1.04   x1.15
              8       x1.17      x0.87    x0.94   x1.04

        Defect correction -- ``x <- x + LU^-1 (b - A x)``, which cannot break
        down and needs one back-substitution per iteration instead of two --
        was implemented and measured, and it is a WASH: x1.06, x0.97, x1.17 on
        the three cases, alternated in one process. It is not here because it
        was not better, and the 3.2x that first appeared to favour it was
        machine load drifting between sequential runs, not the solver.
        """
        dx = self.dx
        A = self._step_matrix(r)
        b = np.broadcast_to(r * dx * dx, (self.nz, self.nx)).ravel().copy()

        def direct():
            self._lu = spl.splu(A, permc_spec=ORDERING)
            self.factorisations += 1
            return self._lu.solve(b)

        if self._lu is None:
            x = direct()
        else:
            n_it = [0]
            P = spl.LinearOperator(A.shape, matvec=self._lu.solve)
            x, info = spl.bicgstab(
                A, b, x0=self._x, M=P, atol=0.0,
                callback=lambda xk: n_it.__setitem__(0, n_it[0] + 1),
                **{_RTOL: self.krylov_tol})
            if info != 0 or n_it[0] > self.max_krylov_iterations:
                x = direct()                       # preconditioner has gone stale
        # Kept unclipped, and kept whether or not it was clipped on the way
        # out, because it is the next step's starting guess and not an answer.
        self._x = x
        return np.clip(x, 0.0, 1.0).reshape(self.nz, self.nx)

    def link_conductivity(self):
        """
        Hydraulic conductivity on every link: ``(vertical, horizontal, wrap)``.

        The matrix CONDUCTS BETTER AS IT DISSOLVES, interpolated geometrically
        -- linearly in the logarithm, which is how conductivity varies --
        between intact granite and fully dissolved rock:

            k(M) = k_matrix(T)^M * k_weathered(T)^(1 - M)

        on the mean of the two cells a link joins. A jointed link keeps
        ``k_fracture``: an open joint is an open joint whatever the rock beside
        it has done.

        This is the model's one positive feedback, and it is not decoration.
        Dissolving rock opens connected porosity, so water is drawn into the
        weathered zone, which dissolves it faster still -- the mechanism behind
        wormholes in limestone and reactive-infiltration fingers generally.
        With the conductance held fixed, as it was until now, the section
        weathers almost uniformly with depth. With the feedback the shallow
        blocks are destroyed while the deeper ones survive and taper, because
        water opens the rock it passes through on the way down and arrives
        saturated below. That is a weathering PROFILE, which is what a
        saprolite actually looks like, and the fixed version cannot produce it
        at all (prototypes/probe_h_evolving_permeability.py).

        The two endpoints are measured, not chosen. Goodfellow et al. (2016)
        report granodiorite matrix conductivity rising three to four orders of
        magnitude across weathering grades, from 9e-9 to 8e-8 cm/s in parent
        rock to 9e-5 to 9e-4 cm/s in the most weathered samples; ``k_matrix``
        and ``k_weathered`` are the mid-points of those ranges. Everything else in
        this module is still a placeholder.
        """
        net = self.network
        lo, hi = np.log(self.k_matrix_at_T), np.log(self.k_weathered_at_T)

        def k_of(m):
            return np.exp(m * lo + (1.0 - m) * hi)

        M = np.clip(self.M, 0.0, 1.0) if self.M is not None \
            else np.ones((self.nz, self.nx))
        kv = np.where(net.link_v, self.k_fracture,
                      k_of(0.5 * (M[:-1, :] + M[1:, :])))
        kh = np.where(net.link_h, self.k_fracture,
                      k_of(0.5 * (M[:, :-1] + M[:, 1:])))
        if net.periodic_x:
            kw = np.where(net.link_wrap, self.k_fracture,
                          k_of(0.5 * (M[:, -1] + M[:, 0])))
        else:
            kw = np.zeros(self.nz)
        return kv, kh, kw

    def flow_operator(self):
        """
        The conductance matrix for the head, and the right-hand side.

        Assembled entirely as triplets, base boundary included. COO sums
        duplicate entries, so the boundary conductance is simply one more
        contribution to the diagonal -- there is no need to reach into an
        assembled matrix and modify a row, which needed a LIL round trip and
        cost 44 ms of a 180 ms initialize() at 22,650 cells.

        Returned rather than kept private so that a test can look at it: it is
        the second matrix ``ORDERING`` claims is structurally symmetric.
        """
        nz, nx, dx = self.nz, self.nx, self.dx
        n = nz * nx
        idx = np.arange(n).reshape(nz, nx)

        kv, kh, kw = self.link_conductivity()

        rows, cols, vals = [], [], []
        pairs = [(idx[:-1, :].ravel(), idx[1:, :].ravel(), kv.ravel()),
                 (idx[:, :-1].ravel(), idx[:, 1:].ravel(), kh.ravel())]
        if self.network.periodic_x:
            # The wrap link closes the section onto itself, so there are no
            # side walls at all. A no-flow wall forces the lateral flow to
            # vanish there, which with subhorizontal joints manufactures a
            # domain-scale circulation and a drainage divide down the middle.
            pairs.append((idx[:, -1], idx[:, 0], kw))
        for a_, b_, k in pairs:
            rows += [a_, a_, b_, b_]
            cols += [a_, b_, b_, a_]
            vals += [k, -k, k, -k]

        # Infiltration into the top row [m2/s per unit thickness].
        rhs = np.zeros(n)
        rhs[idx[0, :]] = self.infiltration * dx

        # Base: the drainage boundary, psi = 0, so H = -depth. Applied as a
        # conductance to an external fixed head rather than by overwriting the
        # row. Overwriting pins the head but destroys continuity IN that row,
        # and the transport step then treats those cells as ordinary ones --
        # which cost about half a percent in the solute balance while the water
        # balance stayed exact, because water is solved and solute is swept.
        self._k_base = np.where(self.network.cell[-1, :],
                                self.k_fracture, self.k_matrix_at_T)
        self._h_base = -(nz - 0.5) * dx - 0.5 * dx
        rows.append(idx[-1, :]); cols.append(idx[-1, :])
        vals.append(self._k_base)
        rhs[idx[-1, :]] += self._k_base * self._h_base

        A = sp.coo_matrix((np.concatenate(vals),
                           (np.concatenate(rows), np.concatenate(cols))),
                          shape=(n, n)).tocsc()
        return A, rhs

    def _solve_head(self, A, b):
        """
        Solve ``A H = b`` for the head, warm-started from the previous head.

        The flow matrix is a conductance Laplacian: **exactly symmetric with a
        positive diagonal**, checked in the tests, so conjugate gradients
        applies. And between two solves the rock has changed by at most
        ``flow_tolerance``, so the previous head is a very good starting guess
        and the previous factorisation is a very good preconditioner.

        That combination is what makes a converged ``flow_tolerance``
        affordable. Re-solving the head is a third of the run at 0.05 and half
        at 0.01, almost all of it factorisation; warm-started CG replaces most
        of those with a handful of back-substitutions. Measured over 2000 kyr
        at 5 cm, against the direct route at the same tolerance:

            flow_tolerance 0.01   13.94 s, 850 factorisations   direct
                                  10.41 s,  21 factorisations   this
                                  max|dM| between them 4.9e-4

        The factorisation is refreshed when CG fails to converge or needs more
        than ``max_head_iterations`` -- either means the preconditioner has
        drifted too far from the operator to be worth keeping. It is a
        preconditioner, so a stale one costs iterations rather than accuracy:
        the convergence test is on the residual, which does not care where the
        iteration started or what preconditioned it.

        NOTE the opposite result for the SOLUTE operator, where keeping a
        stale factorisation across a rebuild is 2.1x SLOWER -- see FRAME (f).
        The two are not symmetric: that one is rebuilt because its diagonal
        moves every step, this one because the conductance field moves rarely.
        """
        if self._flow_lu is None or self._H_prev is None:
            self._flow_lu = spl.splu(A, permc_spec=ORDERING)
            x = self._flow_lu.solve(b)
        else:
            n = {"i": 0}
            P = spl.LinearOperator(A.shape, matvec=self._flow_lu.solve)
            x, info = spl.cg(A, b, x0=self._H_prev, M=P,
                             callback=lambda xk: n.__setitem__("i", n["i"] + 1),
                             maxiter=self.max_head_iterations * 4,
                             # atol=0 so the RELATIVE tolerance governs alone.
                             # scipy's default absolute tolerance would end the
                             # iteration on the scale of b rather than on the
                             # accuracy asked for -- the same defect as
                             # np.allclose's atol, which made several tests in
                             # this suite unable to fail.
                             atol=0.0, **{_RTOL: self.head_tol})
            if info != 0 or n["i"] > self.max_head_iterations:
                self._flow_lu = spl.splu(A, permc_spec=ORDERING)
                x = self._flow_lu.solve(b)
        self._H_prev = x
        return x

    def solve_flow(self):
        """
        Steady Darcy head, and the link fluxes it implies.

        Finite volume on square cells, so the geometric factor is one: the flux
        along a link is ``K * (H_i - H_j)`` in m2/s per unit thickness. ``H`` is
        TOTAL head with elevation already in it -- adding a separate gravity
        term to the link flux double-counts it and manufactures water.

        The head is re-solved whenever the rock has changed by more than
        ``flow_tolerance``, which at a converged tolerance is most of the run's
        cost -- half of it at 0.01. See :meth:`_solve_head`.
        """
        nz, nx, dx = self.nz, self.nx, self.dx
        kv, kh, kw = self.link_conductivity()
        # The tortuosity follows the rock on the SAME trigger, and must: it
        # depends on M, and the transport operator is cached and refactorised.
        # Refreshing it every step would rebuild and refactorise the operator
        # every step. Refreshing it here ties it to the rock-change tolerance
        # that already governs when the conductivity is allowed to move, so
        # the two halves of "weathering opens the rock" advance together.
        self._tort = self.link_tortuosity()

        A, b = self.flow_operator()
        H = self._solve_head(A, b).reshape(nz, nx)
        self.H = H
        self.q_v = kv * (H[:-1, :] - H[1:, :])      # positive downward
        self.q_h = kh * (H[:, :-1] - H[:, 1:])      # positive rightward
        if self.network.periodic_x:
            self.q_wrap = kw * (H[:, -1] - H[:, 0])  # last column -> first
        else:
            self.q_wrap = np.zeros(nz)

        # Per-cell inflows, split by where they come from. Vertical flow is
        # strictly downward (checked in prototypes/probe_d_darcy.py), so the
        # only vertical inflow to a cell is from the row above it.
        self._in_above = np.zeros((nz, nx))
        self._in_above[0, :] = self.infiltration * dx
        self._in_above[1:, :] = np.maximum(self.q_v, 0.0)
        self._in_left = np.zeros((nz, nx))          # from the cell to the left
        self._in_left[:, 1:] = np.maximum(self.q_h, 0.0)
        self._in_right = np.zeros((nz, nx))         # from the cell to the right
        self._in_right[:, :-1] = np.maximum(-self.q_h, 0.0)
        # Across the seam, column nx-1 is "left of" column 0 and vice versa.
        self._in_left[:, 0] += np.maximum(self.q_wrap, 0.0)
        self._in_right[:, -1] += np.maximum(-self.q_wrap, 0.0)
        self.q = self._in_above + self._in_left + self._in_right
        # What leaves the domain through the base, per bottom-row cell.
        self.q_out_base = self._k_base * (H[-1, :] - self._h_base)
        self._M_flow = None if self.M is None else self.M.copy()
        self.flow_solves += 1
        return self.q

    # ---- lifecycle

    def initialize(self):
        """Set the rock fresh and route the flow. Run once, before update()."""
        self._T = None
        self._T_key = None
        self._A = None
        self._diag = None
        self._lu = None
        self._x = None
        self._dt = None
        self._c_held = None
        self._drift = None
        self.rejected_steps = 0
        self.factorisations = 0
        self._M_flow = None
        self.flow_solves = 0
        self.M = np.ones((self.nz, self.nx))
        self.c = np.zeros((self.nz, self.nx))
        self.t = 0.0
        self.solve_flow()
        return self

    def update(self, dt=None, dt_limit=None):
        """
        Advance the rock state by one step, and return the step actually taken.

        Solute transport is solved at steady state -- it is fast next to the
        rock changing -- so there is no storage term and the model carries no
        porosity, and therefore no residence time.

        What the rock loses is what the water gains:

            d(M/M0)/dt = - r (1 - c) / tau

        and ``r`` is proportional to ``M``, because the reactive surface area
        is. Over a step in which ``c`` is held the equation is therefore linear
        in ``M``, and it integrates EXACTLY:

            M(t + dt) = M(t) exp(-lambda dt), lambda = (r / M) (1 - c) / tau

        Forward Euler stood here before, taking the tangent to that
        exponential. The tangent always undershoots, which is the only reason
        the step ever needed clipping at zero, and it was six times less
        accurate at identical cost -- the same solves, one np.exp added.
        ``lambda`` is formed from :attr:`specific_reaction_coefficient` and so
        never divides by a mineral content approaching zero.

        WHAT LIMITS THE STEP. Holding ``c`` across it is the only
        approximation left, so that is what is bounded: after each solve the
        new field is compared with the one that was held, and the next step is
        scaled so the drift lands near ``c_drift_max``.

        ``dx_max`` used to do this job and did it badly. It bounds the change
        in ``M`` in the single fastest-dissolving cell, which is a proxy for
        the drift and not the drift; which cell that is jumps from step to
        step, so the sequence is erratic and the error is NOT MONOTONE in the
        limiter. Measured on the 3 m section at ``dt_max`` = 50 kyr:
        ``dx_max`` 0.05 gave 1.18e-4, 0.10 gave 2.66e-5, 0.20 gave 2.27e-4.
        Nothing can be chosen against a curve like that. Against
        ``c_drift_max`` the error is monotone on every case tried and close to
        linear, so it is a dial:

            c_drift_max   3 m app    45 deg   12 x 9 m   3 m at dx.02
                  0.003   8.1e-06   2.1e-03    3.3e-03        1.2e-05
                  0.01    3.1e-05   8.2e-03    1.1e-02        4.7e-05
                  0.03    1.1e-04   2.5e-02    3.1e-02        1.6e-04
                  0.10    5.4e-04   7.8e-02    9.3e-02        8.1e-04

        The default, 0.03, is a CHOSEN ERROR BUDGET and not a derived
        quantity: about 3 % of full scale on a field that lives in [0, 1],
        which is invisible in the demo and beats the old Euler answer on three
        of those four cases. It is one line to change, and the row above says
        what changing it costs.

        The budget is ENFORCED, not aimed at: a step is taken, the drift it
        actually produced is measured, and a step that overran is thrown away
        and retried shorter. Predicting the step from the previous one's drift
        is cheaper but leaves the FIRST step uncontrolled, and that one step is
        enough to put a floor under the whole run -- measured, an opening step
        of 7124 yr held the error at 1.1e-2 no matter how tight the budget,
        while the same run with the opening step controlled went to 6.3e-5.
        A rejection costs one solve, and only when the prediction was wrong.

        Passing ``dt`` explicitly overrides all of this. ``dt_limit`` instead
        caps the automatic choice without replacing it, which is what
        :meth:`run` uses to land exactly on the time asked for.
        """
        if self._c_held is None:
            self._c_held = self.solve_solute(self.reaction_coefficient)
        c_held = self._c_held
        lam = self.specific_reaction_coefficient * (1.0 - c_held) / self.tau

        def advance(step):
            # The clip is a guard, not a mechanism: lambda >= 0 because c <= 1,
            # so the exponential cannot leave (0, 1] on its own.
            return np.clip(self.M * np.exp(-lam * step), 0.0, 1.0)

        if dt is not None:
            M_new = advance(dt)
            self.M = M_new
            self._c_held = self.solve_solute(self.reaction_coefficient)
            self.c = self._c_held
            self.t += dt
            return dt

        rate = lam * self.M                        # d(M/M0)/dt [1/s]
        want = min(self.dt_max, self.dx_max / max(rate.max(), 1e-30))
        if self._dt is not None:
            # PREDICT from the drift the last step actually produced, rather
            # than reaching for dt_growth every time. Drift is first order in
            # the step, so this lands near the budget and the rejection below
            # becomes the safety net it is meant to be. Reaching for the growth
            # cap instead means overshooting and being rejected on nearly every
            # step: measured, 73 rejections for 75 accepted steps, which is two
            # solves per step and gave back the whole saving.
            predicted = self._dt * min(
                self.dt_growth,
                0.9 * self.c_drift_max / max(self._drift, 1e-300))
            want = min(want, predicted)
        step = want if dt_limit is None else min(want, dt_limit)

        while True:
            M_new = advance(step)
            saved_M, self.M = self.M, M_new
            c_new = self.solve_solute(self.reaction_coefficient)
            self.M = saved_M
            drift = np.abs(c_new - c_held).max()
            if drift <= self.c_drift_max or step <= self.dt_min:
                self._drift = drift
                break
            self.rejected_steps += 1
            step = max(self.dt_min,
                       0.9 * step * self.c_drift_max / max(drift, 1e-300))

        # Remember what the control would have allowed, undoing any cap the end
        # of the run imposed -- otherwise a short final step would be read as
        # evidence that short steps are needed, and the answer would depend on
        # how the run was chopped into calls.
        self._dt = step if dt_limit is None or step < dt_limit else want

        self.M = M_new
        self._c_held = c_new
        self.c = c_new
        self.t += step

        # Re-solve the head once the ROCK has changed enough to have moved the
        # conductance, not every so many steps. A step count would tie the
        # answer to the step size, and the drift control exists precisely so
        # that it does not: halving the budget would silently double how often
        # the flow was updated and change the result. This trigger is a
        # property of the rock, so refining the time step converges rather
        # than wanders.
        if self._M_flow is not None and \
                np.abs(self.M - self._M_flow).max() >= self.flow_tolerance:
            self.solve_flow()
            self._c_held = None           # the flow moved; c must be re-solved
        return step

    def run(self, years):
        """
        Advance to ``years`` of model time, initializing if needed.

        Lands ON the time asked for. It used to step past it by up to one step
        -- ``while t < target: update()`` and nothing trimming the last one --
        so ``run(years=30e3)`` returned a model at 32.5 kyr when the steps were
        long. That is invisible in any single run and poisonous in a
        comparison: two settings would be compared at two different times, and
        the difference read as the error of the coarser one. It put a floor of
        about 1e-2 under a convergence study that should have gone to zero, and
        it is the reason that study looked non-monotone.
        """
        if self.M is None:
            self.initialize()
        target = years * YEAR
        while self.t < target - 1e-9 * YEAR:
            self.update(dt_limit=target - self.t)
        return self

    def finalize(self):
        """Nothing to release; present so the lifecycle is complete."""
        return self

    # ---- output

    @property
    def dissolved_fraction(self):
        """Fraction of the soluble phase that has dissolved, in [0, 1]."""
        return 1.0 - self.M

    @property
    def is_grus(self):
        """
        Cells past the arbitrary ``x_grus`` cut-off.

        NOT a statement that the rock is grus. See ``x_grus``: the weathering
        grades are a matter of fabric and mineralogy, and this is one number
        thresholded.
        """
        return self.dissolved_fraction > self.x_grus

    @property
    def is_corestone(self):
        """
        Cells below the arbitrary ``x_core`` cut-off.

        NOT a statement that a cell is part of a corestone, which is a SHAPE
        -- a rounded block surrounded by weathered rock. Intact bedrock at
        depth passes this test and is not a corestone. See ``x_core``.
        """
        return self.dissolved_fraction < self.x_core

    @property
    def darcy_speed(self):
        """
        Speed of the water at each cell centre [m/s].

        The solver works on LINKS -- a flux through each face -- and this is
        the cell-centred vector reassembled from them and its magnitude taken:
        the mean of the two vertical faces and of the two horizontal faces,
        divided by the cell size to turn a flux per unit thickness into a
        specific discharge. At the surface the face flux is the infiltration
        and at the base it is the drainage, so no cell is missing a face.

        This is the CAUSE the rest of the model is the effect of, and unlike
        the affinity field it is not a restatement of the dissolved fraction:
        early on the matrix carries about a thousandth of the mean infiltration
        rate while the joints carry twenty times it, and as the rock opens the
        matrix comes up to carry nearly all of it.
        """
        nz, nx, dx = self.nz, self.nx, self.dx
        qz = np.zeros((nz, nx))
        qz[:-1, :] += self.q_v
        qz[1:, :] += self.q_v
        qz[0, :] += self.infiltration * dx
        qz[-1, :] += self.q_out_base
        qx = np.zeros((nz, nx))
        qx[:, :-1] += self.q_h
        qx[:, 1:] += self.q_h
        qx[:, -1] += self.q_wrap
        qx[:, 0] += self.q_wrap
        return np.hypot(0.5 * qx, 0.5 * qz) / dx

    @property
    def affinity(self):
        """The bracket, ``1 - C/C_eq``: how much capacity the water has left."""
        return 1.0 - self.c

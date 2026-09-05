"""
Probe I: what regime is the model in once oxygen drives it?

OPEN. Run before writing step 2 of design/08-BUILD.md, and it changed what
step 2 means. Nothing here is implemented; this file is the arithmetic that
says what implementing it would do.

Design 08 presents the switch from silica to oxygen as a re-pointing: same
operator, same assembly, right-hand side to zero, inlet to one. That is true
of the CODE. It is not true of the physics, and this probe is what found the
difference.

THE THREE FINDINGS

1.  The regime inverts. On silica the section-scale Damkohler is 6.56 --
    saturation-limited, water arrives at the base within exp(-6) of the
    ceiling, and that is what shelters a block interior. On oxygen it is
    0.004 to 0.11 across the whole parameter bracket, which is
    REACTION-LIMITED by the model's own classifier. Water crosses the section
    barely touched; oxygen is everywhere at once.

    So the module docstring's headline -- "Corestones are a
    saturation-limited phenomenon" -- is false under oxygen, and the
    sheltering mechanism is a different one. A corestone survives because O2
    cannot DIFFUSE into intact matrix, not because the water arrived spent.
    Both are "the water never got there"; they are not the same sentence and
    the exercise is built on the first.

2.  The front the exercise exists to show is at or below the grid scale. The
    reaction-diffusion boundary layer into fresh rock is sqrt(D_fresh / r) =
    2.4 to 13 cm, centre about 6. The demo offers 5, 2.5 and 2 cm cells, so
    the front is 0.5 to 2.7 cells wide at the coarse end. Cell size is a
    student-facing control, and the page promises it is the numerical grid
    and not the rock.

3.  The mechanical feedback is load-bearing, not decoration. Through fresh
    rock at tortuosity 1e4 the front advances at 0.07 to 0.41 m/Myr -- slower
    than the silica model's measured 0.81 and one to two orders below the
    4-7 m/Myr the field gives. The same expression at the weathered
    tortuosity of 10 gives 5.2 m/Myr, inside the field range. Cracking is the
    difference between those, so design 08 step 3 is not a garnish on step 2;
    it is what decides whether the answer is right.

    CAVEAT, and it is not small: at tortuosity 10 the penetration length is
    1.87 m against a 3 m section, so the sharp-front assumption behind the
    velocity expression has failed and 5.2 m/Myr is a scaling indication
    rather than a front rate. What survives is the RATIO -- cracking is worth
    up to about thirty-fold -- not the endpoint.

WHAT IS ASSUMED, MARKED
    k_ox      from the literature rate ~1e-13 mol m-2 s-1 at 0.25 mM O2 and
              25 C, good to a factor of three (design 08; White & Yee 1985
              remains unread). Bracketed 0.33x to 3x here.
    A         6 phi / d, the model's own geometric convention, on the BIOTITE
              volume fraction. phi_biotite IS A GAP: 3 to 10 % is the range
              for granite and nothing here narrows it. Bracketed.
    D_O2      2.1e-9 m2/s in free water. Not the model's silica value; O2 is
              about twice as mobile.
    tortuosity, grain size, infiltration, f_FeO, V_FeO: the model's own.

This is a one-dimensional scaling analysis and not the model. It says what to
expect and what to resolve; it does not say what the model will do.

Run:  PYTHONPATH=src python3 prototypes/probe_i_oxygen_regime.py
"""

import numpy as np

from corestone.weathering import oxygen_solubility, YEAR

#: Free-water diffusivity of dissolved O2 [m2/s] at 25 C. About twice the
#: silica value the model carries, because O2 is a small neutral molecule.
D_O2 = 2.1e-9

#: The literature rate and its factor-of-three bracket [mol m-2 s-1 at
#: 0.25 mol/m3], and the biotite volume fractions that span granite.
RATES = (0.33e-13, 1.0e-13, 3.0e-13)
PHI_BIOTITE = (0.03, 0.05, 0.10)

#: The model's own values, repeated here rather than imported, so that this
#: file states every number it uses.
Q = 0.30 / YEAR                  # infiltration [m/s]
GRAIN = 2.0e-3                   # grain diameter [m]
TORT_FRESH = 1.0e4
TORT_WEATHERED = 10.0
F_FEO = 0.011
V_FEO = 12.00e-6                 # [m3/mol Fe]
DEPTH = 3.0                      # the demo's section [m]
T_REF = 285.0                    # [K]
CELL_SIZES = (0.05, 0.025, 0.02)  # what the demo offers [m]


def reaction_coefficient(rate, phi):
    """``r = k_ox A`` [1/s], from a measured rate and a mineral abundance."""
    return (rate / 0.25) * (6.0 * phi / GRAIN)


def report():
    C_O2 = oxygen_solubility(T_REF)
    tau_ox = 0.25 * F_FEO / (V_FEO * C_O2)
    print("oxygen at %.2f K: %.4f mol/m3, tau_O2 = %.0f\n" % (T_REF, C_O2, tau_ox))

    print("           the two length scales, and what they cost")
    print("  k_ox [m/s]   A [m2/m3]   Da(3 m)   front [cm]   cells@5cm   "
          "v [m/Myr]")
    lo_pen, hi_pen, lo_v, hi_v, lo_da, hi_da = 9e9, 0.0, 9e9, 0.0, 9e9, 0.0
    for rate in RATES:
        for phi in PHI_BIOTITE:
            r = reaction_coefficient(rate, phi)
            D_fresh = D_O2 / TORT_FRESH
            pen = np.sqrt(D_fresh / r)          # reaction-diffusion layer [m]
            v = np.sqrt(D_fresh * r) / tau_ox   # sharp-front speed [m/s]
            da = DEPTH * r / Q                  # advective Damkohler [-]
            print("  %.2e     %6.0f    %7.4f   %8.2f   %9.1f   %9.3f"
                  % (rate / 0.25, 6.0 * phi / GRAIN, da, 100 * pen,
                     pen / CELL_SIZES[0], v * YEAR * 1e6))
            lo_pen, hi_pen = min(lo_pen, pen), max(hi_pen, pen)
            lo_v, hi_v = min(lo_v, v), max(hi_v, v)
            lo_da, hi_da = min(lo_da, da), max(hi_da, da)

    print("\n  Damkohler        %.4f to %.4f   -- silica is 6.56"
          % (lo_da, hi_da))
    print("  front width      %.1f to %.1f cm" % (100 * lo_pen, 100 * hi_pen))
    print("  front speed      %.3f to %.3f m/Myr  -- silica model 0.81, "
          "field 4-7" % (lo_v * YEAR * 1e6, hi_v * YEAR * 1e6))

    print("\n  how many cells wide is the front, at each size the demo offers?")
    for dx in CELL_SIZES:
        print("    dx = %4.1f cm   %.1f to %.1f cells"
              % (100 * dx, lo_pen / dx, hi_pen / dx))

    print("\n  and once it has cracked -- tortuosity 1e4 -> 10:")
    r = reaction_coefficient(RATES[1], PHI_BIOTITE[1])
    D_w = D_O2 / TORT_WEATHERED
    pen_w, v_w = np.sqrt(D_w / r), np.sqrt(D_w * r) / tau_ox
    print("    penetration %.2f m, speed %.1f m/Myr -- a factor of %.0f"
          % (pen_w, v_w * YEAR * 1e6,
             v_w / (np.sqrt(D_O2 / TORT_FRESH * r) / tau_ox)))
    print("    BUT the penetration is now %.2f m against a %.0f m section, so"
          % (pen_w, DEPTH))
    print("    the sharp-front assumption has failed and the speed is a")
    print("    scaling indication. The RATIO is the finding; the endpoint is not.")


if __name__ == "__main__":
    report()

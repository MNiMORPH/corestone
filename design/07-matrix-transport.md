# 07 -- Matrix transport: evolving tortuosity, and the scale of the dispersivity

## The question

Two terms carry solute through the matrix, and both were wrong in the same
direction -- they let fresh rock transport far too freely, which is exactly
where and when the weathering rind forms.

1. `tortuosity` was fixed at 10, a WEATHERED value: saprolite at ~30 % porosity
   gives `D_eff/D_0 ~ 0.1`. Intact granite measures 2e-14 to 1.3e-12 m2/s
   against a free-water 1e-9, so a tortuosity of 1e3 to 1e5, centre 1e4. The
   conductivity already interpolates with `M` and this did not.

2. `dispersivity` was 0.05 m, from Gelhar, Welty & Rehfeldt (1992), whose
   scaling is roughly a tenth of the transport distance. That is a
   MACROdispersivity: it stands in for metre-scale heterogeneity in the flow
   paths. This model *draws* that heterogeneity, as the joint network, so
   using the field-scale value in the matrix counts the joints twice.

## Options, each with its full cost

**(a) Leave both.** Free, and the model keeps a rind about a metre thick where
real rinds are centimetres. The known cheat stays documented.

**(b) Evolving tortuosity alone.** Rejected on the probe below: it hands the
fresh matrix to the dispersion term, 4:1, using a bulk-aquifer dispersivity at
1e-11 m/s. That moves the error rather than removing it.

**(c) Both.** Tortuosity interpolates with `M` as `k` does; the dispersivity
becomes pore-scale, tied to the grain size already used for the reactive
surface area. Molecular diffusion then governs the matrix at both ends and the
mechanism stays diffusive. Cost: the model slows again (t90 1191 -> 3755 kyr),
so the animation cadence and the Show cap both move, and every measured number
in the demo comments and on the exercise page is re-measured. Fourth such pass.

## The probe

`prototypes/probe_i_evolving_tortuosity.py`, and a second part run against it.

Which term governs the matrix, median matrix |v| = 5.82e-12 m/s:

    dispersivity  tortuosity       dispersive    molecular   mol/disp
    0.05 m        10 weathered      2.911e-13    6.923e-11     237.86
    0.05 m        1e4 fresh         2.911e-13    6.923e-14       0.24
    2 mm          10 weathered      1.164e-14    6.923e-11    5946.62
    2 mm          1e4 fresh         1.164e-14    6.923e-14       5.95

So (b) inverts which term rules the fresh matrix and (c) does not.

Sensitivity to the fresh tortuosity, which spans two orders in the literature:
t90 moves only 1.5x from 1e3 to 1e5 (2256 to 3329 kyr). That looks like
robustness and is not -- it is dispersion having taken over, so lowering
molecular diffusion further changes little. Under (c) the same sweep is the
one that matters and molecular diffusion is still in charge.

What it costs and buys:

                                     t90 kyr    front       rind
    as shipped (fixed 10, 0.05 m)       1191   2.52 m/Myr   30.6 %
    evolving tortuosity, 0.05 m         2925   1.03          5.8 %
    evolving tortuosity, 2 mm           3755   0.80          4.6 %

## Decision, and why

**(c).** The rind is what this model exists to show and it goes from about a
metre of gradational zone to fourteen centimetres, which is the scale real
rinds have.

The awkward part, recorded rather than buried: correcting the transport makes
the predicted weathering front WORSE against field data, from 2.52 m/Myr to
0.80 against 4-7 measured for temperate granite regoliths (White et al. 2001).
That is accepted and reported rather than tuned away. The gap lives in the
reactive surface area, which is geometric here at 900 m2/m3 while BET for
granite is 3e5-3e6; closing it needs a factor of five, still five hundred
times below BET. A discrepancy sitting inside a range the field itself has not
resolved (White & Brantley 2003) is a better thing to show a student than a
match obtained by choosing a surface area.

It also corrects a claim that was on the live page. At 30.6 % of the section
part-dissolved, the shipped model was not advancing a front at all -- dividing
3 m by t90 and calling it a front rate is only meaningful once the front is
sharp. The number that looked like agreement was partly an artefact of that.

**Unverified, and marked as such:** the pore-scale dispersivity being of order
the grain diameter is standard in the porous-media literature, and Perkins &
Johnston (1963) is the canonical review, but I could not read a source stating
it in the session where this was decided. Gelhar is verified for the
field-scale value it replaces.

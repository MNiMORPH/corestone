# 13 -- Why the model is slow, measured lever by lever

Written 2026-09-23, after design 12 restored the temperature correction and
made the discrepancy worse. **Nothing here is changed in the model.** It
records which levers were measured, what each delivers, and which of them is
defensible on evidence rather than on hitting a target.

Read design 11 (the ceiling) and design 12 (the reference temperature) first;
this is the measurement that sits on top of both.

## 1. The gap, restated honestly

Design 12 halved the model's speed twice over. At the defaults the model reacts
90 % of a 3 m section in about 15 000 kyr, against measured temperate fronts of
4 m/Myr (Davis Run) and 7 m/Myr (Panola), which would take a 3 m section
430-750 kyr. **The model is 20-35x slow.**

Two things first, because both were wrong in this file before today:

**`3 m / t90` is not a front rate.** At t90, 91 % of the section is
part-reacted, at 1.5, 3 and 6 m alike, and `D/t90` rises with the depth drawn
-- 0.333, 0.444, 0.706 m/Myr. The profile shows why: the top cell goes to
M = 0.02 and everything below sits flat at the joint fraction. The joint
network dissolves at every depth at once, and nothing propagates downward.

**The gap is NOT in the reactive surface area**, which this file used to assert.
Raising the area a thousandfold raises the rock consumed in 300 kyr by 1.72.
Water in a joint does not touch the block interior however reactive that
interior is, so area cannot buy what contacting does not deliver.

## 1a. The comparison was between two different quantities (2026-09-23)

Andy: *"Are you comparing directly field rates against weathering rates in
fractures and at the surface?"* **No -- and worse.** `3 m / t90` is BULK
dissolution of a section that is mostly inert block. A field weathering front
is a BOUNDARY advancing through rock that is reacting. Those are different
quantities, and section 1 above had already rejected the first one.

Measured properly -- how fast the reacted rind eats into a block from its
bounding joint, which IS a front:

| time | rind depth (dissolved = 0.5) | advance |
|---|---|---|
| 400 kyr | 0.050 m | 0.25 m/Myr |
| 1600 kyr | 0.100 m | 0.06 m/Myr |
| 3200 kyr | 0.150 m | 0.03 m/Myr |

**It decelerates, as sqrt(t).** That is not a calibration problem; it is the
signature of diffusion into a fixed geometry, and it is exactly the case
Fletcher et al. (2006) show cannot sustain a steady front. Their fractured case
advances at constant rate; their unfractured case goes parabolic and falls
behind. **This model is the unfractured case**, whatever joints it was seeded
with, because the joints never change.

**The mechanism is diagnosed but not wired.** `cracking_number`,
`fracture_energy` and `cracking_threshold` exist, are computed from Goodfellow's
criterion (design 10) and are printed in `thermo_report` -- and nothing writes
to `link_v` or `link_h` after seeding. The model knows when the rock should
crack and never cracks it.

**So "the model is 20-35x too slow" was the wrong frame.** The model is doing
something structurally different from the field sites: a fixed-geometry
diffusion problem against a self-fracturing front. Closing the gap by moving
`C_eq`, the surface area or the conductivity would be fitting a constant-rate
observation with a decelerating model. The levers in section 3 are still
correctly measured; what changes is that none of them addresses this.

## 2. The identity that organises everything

Section-integrated, the rate of rock consumption is

    v = q C_eq phi / N_0

where `phi` is the flux-weighted saturation of the water leaving the base -- the
**contacting efficiency**. Measured at the defaults: `q C_eq / N_0` = 3.41 m/Myr
is the supply ceiling, `phi` = 0.09 to 0.14, and the product is the model's
actual rate. **1/phi is most of the gap, and it is a measurement, not a
hypothesis.**

Every lever below acts on one of the three factors, and they are not
independent: raising `C_eq` lengthens `L = q C_eq/(k A)`, so the water takes
longer to fill and `phi` FALLS. That is why the ceiling alone buys so little.

## 3. What each lever delivers

Rock consumed in 300 kyr, a 3 m section at 1 m joint spacing, 12 C, relative to
the baseline; `phi` and the effluent silica beside it. **Two metrics on
purpose:** a lever that raises the rate while leaving the water clean has not
touched the defect.

| lever | rate | phi | effluent |
|---|---|---|---|
| baseline | 1.0x | 0.143 | 0.47 mg/L |
| ceiling x4.3 alone | 1.7x | 0.058 | 0.81 |
| ceiling x7.7 alone | 2.1x | 0.044 | 1.11 |
| area x5 alone | 1.4x | 0.115 | **0.38** (falls) |
| tortuosity_fresh 1e4 -> 10 | 2.1x | 0.397 | 1.29 |
| K_sat_intact x100 | 1.8x | 0.339 | 1.10 |
| K_sat_intact x1e4 | 4.4x | 0.866 | 2.82 |
| joint spacing 1.0 -> 0.25 m | 2.7x | 0.441 | 1.44 |
| ceiling x4.3 + area x4.3 (L held) | 2.5x | 0.096 | 1.34 |
| ...+ K_sat_intact x100 | 7.0x | 0.303 | 4.24 |
| ...+ K_sat_intact x1e4 | 13.3x | 0.197 | 2.77 |

**Caveat on the metric:** rock-consumed-in-300-kyr saturates as a run
approaches completion, so the 13.3x row (93 % consumed) understates its true
factor. It is a cheap proxy, not a front rate -- see section 1.

**Only the flow levers raise the effluent.** Area LOWERS it. That is the
signature of a contacting problem rather than a kinetic one.

## 4. The one conclusion that does not depend on the rate at all

    model C_eq at 12 C           3.26 mg/L SiO2
    granitic groundwaters       14-25 mg/L SiO2

**The model's water cannot reach field concentrations even at perfect
contact**, because its ceiling sits below them. This is a comparison against
measured water chemistry, not a fit to a rate, and it is the only lever here
that is already evidenced. It needs a factor of 4.3 to 7.7 -- which is design
11's range, arrived at from the water side instead of the thermodynamic side.

## 5. The trap in "make it match"

Ceiling x4.3 + area x4.3 + `K_sat_intact` x100 gives ~7x, which puts t90 near
550 kyr against a field 430-750 kyr. **It matches.** Its effluent is
4.24 mg/L, still far under 14-25.

So the model can be made to match the RATE without matching the WATER, by being
fast rather than by being well contacted. The effluent is the number that
discriminates between a model that is right and a model that has been aimed,
and it is the one to hold this model to. Choosing parameters to make the rate
come out is still the move that would make the number meaningless.

## 6. What needs evidence before anything moves

| lever | status |
|---|---|
| **ceiling x4.3-7.7** | **evidenced** -- field waters exceed the model's ceiling. Value still needs a source: what silica activity does plagioclase-kaolinite equilibrium actually fix? |
| `K_sat_intact` 5e-10 | **open.** That is LABORATORY-INTACT granite. Rock at a weathering front is microfractured by exfoliation, biotite oxidation and stress release. Is there a measured value for the transition zone? |
| area x4.3 | **open.** Geometric 900 m2/m3 against BET 3e5-3e6. Is there a principled intermediate, and what do published reactive-transport models of granite actually use? |
| the comparison target | **open.** What did White et al. (2001) measure at Panola and Davis Run, over what interval, and is a 3 m jointed section with rain entering at zero silica the right analogue for a front under metres of loaded regolith? |

**Nothing is decided.** A literature chase on the four is running; its answers
belong in this file, and only then is there a change to propose.

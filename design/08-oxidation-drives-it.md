# 08 -- Biotite oxidation sets the pace, and its expansion cracks the rock

## The question

This model dissolves plagioclase, limited by silica approaching quartz
saturation, and indexes every rock property on how much plagioclase is left.
The literature says that is the wrong driver.

Buss, Sak, Webb & Brantley (2008), abstract, on Rio Icacos:

> "Chemical, petrographic, and spectroscopic evidence demonstrates that
> **biotite oxidation is the most likely fracture-inducing reaction**. ... In
> the corestone, biotite oxidation induces spheroidal fracturing, facilitating
> the influx of fluids that react with other minerals, **dissolving
> plagioclase** and chlorite, creating additional porosity"

They imaged oxidised biotite 2.7 cm inside nominally fresh corestone and found
**no** plagioclase weathering there. Goodfellow et al. (2016) -- the paper this
model already takes k_matrix and k_weathered from -- says the same in its own
abstract, and adds the sentence that matters most here: "**Major changes in
rock properties can occur with only minor element leaching.**"

So: oxidation goes first, cracks the rock, and dissolution follows through the
cracks. We have the sequence backwards.

## Options, each with its full cost

**(a) Leave it.** The transport physics is right and the rind is the right
shape; only the label on the reaction is wrong. Cheap, and dishonest in a way
the page cannot admit without undermining itself.

**(b) Swap the driver: one solute, O2 instead of silica.** M becomes the
fraction of biotite Fe(II) still unoxidised. The existing k(M) and
tortuosity(M) interpolations then index on the right variable -- which is
Goodfellow's own independent variable, since their conductivity series is
against Fe(III), not against leaching. Cost: the solute flips from product to
reactant, so the transport right-hand side goes to zero and the inlet boundary
becomes c = 1; the rate becomes proportional to c M rather than (1 - c) M. The
operator, the sparse assembly and the exact exponential integrator all survive.
Every measured number in the demo and on the page is re-measured. Again.

**(c) (b) plus plagioclase as a second solute.** The honest two-stage picture:
oxidation paces and opens, dissolution removes the mass and makes the grus.
Doubles the transport solve, and the second stage adds nothing the pictures
show. Recorded as the physically complete option and not taken.

## The mechanical feedback, without fracture mechanics

The expansion has to crack the rock or none of this works, and the whole
fracture problem is out of scope for a teaching model. An **energy budget**
is enough, and it is Fletcher, Buss & Brantley's own criterion: fracture when
the elastic strain energy stored by the expansion reaches the energy needed to
make new surface.

    eps(x) = phi_FeO * beta * x        volumetric strain, x = fraction oxidised
    U(x)   = 0.5 * K * eps^2           elastic energy stored per unit volume
    U_c    = 2 gamma / d               energy to open one grain-scale crack
    N(x)   = U / U_c                   cracks the budget can pay for

Four lines, no crack geometry, no stress field, no failure criterion beyond an
energy balance. `d` is already a model parameter -- the grain size that sets
the reactive surface area and the dispersivity -- so the mechanism costs one
new length scale of nothing.

Inputs, with how well each is known:

| | value | source |
|---|---|---|
| `phi_FeO` | 0.05 | Fletcher et al. (2006) Table 1, volume fraction of the FeO component |
| `V_FeO` | 12 cm3/mol | Fletcher et al. (2006) Table 1 |
| `V_ferrihydrite` | 21.9 cm3/mol Fe | 168.70 g/mol per 2 Fe at 3.85 g/cm3 (mindat). **Needs a better source** |
| `beta` | 0.83 | derived from the two molar volumes above |
| `K` | 20-40 GPa | granite E = 30-60 GPa at nu = 0.25. **Needs a source** |
| `gamma` | ~1 J/m2 | theoretical estimate for silicate minerals, quartz 1.16, albite 0.93, orthoclase 0.89. **Citation needs pinning** -- retrieved from a Tromans & Meech (2002) PDF whose summary returned Lawn & Marshall (1979); one of those two is right and I have not established which |

**It behaves, and it predicts the observation rather than assuming it.**
Cracking begins at x_c = 0.005 to 0.013 across the whole plausible range of
gamma and K, and by x = 0.1 the budget pays for 57 to 379 grain-scale cracks
per grain volume -- pervasive cracking. That is Goodfellow's "~10 % rise in
Fe(III) gives one to three orders of magnitude in K", and their "major changes
with only minor element leaching", coming out of an energy balance rather than
being written in. It also replaces the smooth log-linear k(M) ramp with a
threshold, which is what the percolation literature reports (Navarre-Sitchler
et al. 2009: effective porosity stays at zero until total porosity reaches
about 9 %).

x_c varies by less than a factor of three over gamma = 0.9-3.0 J/m2 and
K = 20-40 GPa, so it is not a tuning knob in disguise.

## Why this is worth the fourth re-measurement

    tau, plagioclase / silica (today)   47744   ceiling on the front  6.3 m/Myr
    tau, biotite Fe(II) / O2             3086                       97.2 m/Myr
    field, temperate granite                                         4-7 m/Myr

Oxygen is fifteen times less limiting per unit rock, so the brake stops being
the solute budget and becomes what Fletcher's model says it is: diffusive
penetration of O2 into the corestone. The model currently reaches only 12 % of
its own stoichiometric ceiling, which is the real discrepancy and is not a
surface-area problem.

And the teaching line gets better rather than worse: **a corestone is not
tougher rock -- the oxygen never got there.**

## Decision

**(b), with the energy budget.** (c) is recorded as the complete picture and
declined: the second solute doubles the cost and changes no picture.

## Still to settle before code

- Does the colour bar become oxidation extent? It stops meaning "soluble phase
  dissolved", which changes what the exercise is about. Andy's call.
- E_a for biotite oxidation, and the rate constant. Not in Palandri & Kharaka,
  which is a dissolution compilation. Buss et al. (2008) measured 8.2e-14
  mol m-2 s-1 at Rio Icacos -- a rate, not an Arrhenius pair.
- Whether the page keeps a plagioclase sentence saying what the model does not
  track, since that is what actually makes the sand.

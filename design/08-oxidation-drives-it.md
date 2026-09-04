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
fracture problem is out of scope. An **energy budget** is enough, and it is
Fletcher, Buss & Brantley's own criterion: fracture when the elastic strain
energy stored by the expansion reaches the energy needed to make new surface.

    eps(x) = x * f_FeO * (dV/V) / 3     LINEAR strain; the 3 converts from
                                        volumetric, Fletcher eq. 12
    U(x)   = E eps^2 / (1 - nu)         Fletcher eq. 11, laterally confined
    U_c    = 2 Gamma / d                energy to open one grain-scale crack
    N(x)   = U / U_c                    cracks the budget can pay for

Four lines, no crack geometry, no stress field, no failure criterion beyond an
energy balance. `d` is already a parameter -- the grain size that sets the
reactive surface area and the dispersivity -- so this costs no new length
scale.

**The energy form is Fletcher's, not the obvious one.** The natural guess is
`0.5 K eps^2` with a volumetric strain, and it is wrong here: the oxidising
rind is laterally confined by the rock around it, so the strain parallel to
the front is zero and the modulus is the constrained combination `E/(1-nu)`,
not the bulk modulus. Fletcher states the constraint explicitly. Getting this
wrong changes which modulus the model needs.

| | value | source |
|---|---|---|
| `f_FeO` | **0.011** | granite. USGS reference granites G-1/G-2/G-3 give 0.0083-0.0116 and Goodfellow's granodiorite 0.0107. Fletcher's 0.05 implies 10.9 wt% FeO and is too high even for his own rock, which measures 4.6 |
| `V_FeO` | 12.00 cm3/mol | Robie & Hemingway (1995), USGS Bulletin 2131, p. 16 |
| `V_goethite` | 20.82 cm3/mol | same table, same page |
| `dV/V` | **0.735** | the two above. Both are single-Fe formula units, so no factor-of-2 ambiguity |
| `E`, `nu` | 76 GPa, 0.24 | SKB R-05-83 Tables 2-7/2-9: 52 intact granite-granodiorite cores, Forsmark. Near-surface cracked granite is softer -- 40 GPa is the more defensible choice here, and the model is about rock that is cracking |
| `x_c` | **0.10** | **calibration.** The oxidation extent at which cracking begins, from Goodfellow et al. (2016). Implies Gamma = 7.3 J/m2, against a floor of 1-2 (Brace & Walsh 1962) and a specimen-scale ceiling of 200-500 (Friedman et al. 1972) |

**The product is goethite, not ferrihydrite.** Fletcher's footnote quotes 0.7
for "wustite to ferrihydrite" citing Robie & Hemingway -- but that compilation
contains no ferrihydrite entry at all, and 0.7 is only reproducible as
FeO to goethite (0.735). Lebedeva & Brantley (2020), the same group, later
write the reaction to goethite explicitly with the identical volumes. Using a
ferrihydrite molar volume would be following a citation its own source does
not support.

**Which fracture energy: the bracket, not a number.** This was chased hard and
it does not close. Gamma is bounded below by the thermodynamic floor -- twice
the specific surface free energy, and Brace & Walsh (1962) measured quartz at
0.40-1.00 J/m2 by cleavage, muscovite at 0.375, and put common minerals at
0.1-2. It is bounded above by specimen-scale fracture energy, the 200 J/m2
Fletcher uses from Friedman, Handin & Alani (1972), which is overwhelmingly
process-zone dissipation rather than surface creation.

At granite f_FeO = 0.011 that bracket is enormous:

    Gamma   1 J/m2  (2 gamma, thermodynamic floor)   x_c =  3.7 %
    Gamma   2       (Goodfellow's high end)                 5.2 %
    Gamma  20                                              16.6 %
    Gamma 200       (specimen-scale G_c)                   52.5 %

**Why it cannot be narrowed, stated rather than fudged.** Goodfellow justify a
near-thermodynamic value by asserting that crack-tip yielding "occurs at a
scale smaller than the constituent grains". Zang, Wagner, Stanchits & Janssen
(2000), JGR 105:23651-23661, MEASURED the process zone in Aue granite at
**2 to 9 grain diameters** -- larger than a grain, not smaller. So the
justification for the low end is contradicted by direct measurement on
granite, while the high end is a specimen-scale quantity that a grain-scale
crack cannot possibly develop. The truth is in between and nothing locates it.

**So the model anchors on the measurement and reports the fracture energy as
a by-product.** Rather than choose a Gamma from a 200-fold range and predict a
threshold, take the threshold from the one direct observation there is --
Goodfellow's ~10 % rise in Fe(III) accompanying one to three orders of
magnitude in hydraulic conductivity -- and report what fracture energy it
implies:

    x_c = 10 %   =>   Gamma = 7.3 J/m2

That number is checkable, and it checks out:

    2 gamma, thermodynamic floor              1 - 2 J/m2   Brace & Walsh (1962)
    Fletcher's own single-crystal estimate    2 - 5        0.01 x their 200-500
    THIS MODEL                                7.3
    specimen-scale G_c                      200 - 500      Friedman et al. (1972)

One order above the floor and one to two below specimen scale, which is where
a crack spanning a process zone measured at 2 to 9 grain diameters (Zang et
al. 2000) should sit -- above pure surface creation, far below a mature
fracture's full dissipation. And within about 1.5x of Fletcher's own
single-crystal estimate, arrived at independently. That agreement is a
consistency check, not a fit: nothing in this chain was adjusted to produce
it.

**This is calibration, and it is labelled as calibration.** x_c is the one
number in the cracking criterion taken from an observation rather than derived,
and the page must say so. What it buys is that the uncertainty now sits in one
visible place instead of propagating through five inputs, and that the model
makes a statement about fracture energy which a reader can check against the
rock-mechanics literature.

An earlier draft of this document claimed the opposite direction -- that the
criterion PREDICTS 12 % against Goodfellow's 10 %. That rested on Fletcher's
f_FeO = 0.05, wrong for granite by 4.5x and too high even for his own rock,
and on a Gamma taken from one end of the bracket. Both halves were errors
pulling opposite ways, so the product looked right. Recorded because it was
the strongest sentence here and it was false.

Two documents would close this and neither is open access: Zang et al. (2000)
in full, and its companion Janssen, Wagner, Zang & Dresen (2001), Int. J.
Earth Sciences 90:46-59, which almost certainly underlies Zang's fracture
surface energy estimate. Worth an interlibrary request alongside White & Yee
(1985) for the oxidation rate.

**Rate law: first order in O2, not Fletcher's C^0.25.** Fletcher's exponent is
the stoichiometric quarter of eq. (2) adopted as a concentration exponent,
with no citation and no experiment behind it. Lebedeva & Brantley (2020) eq.
(15), same group, same reaction, is linear in dissolved O2 and keeps the
quarter where it belongs, in the stoichiometry.

**There is no activation energy, and that is a finding.** No measured Ea
exists for aqueous oxidation of structural Fe(II) in a silicate by dissolved
O2. Hogg & Meads (1975), a dedicated Mossbauer kinetics study of biotite
Fe(II) oxidation, contains no Arrhenius parameters anywhere -- verified
directly: 4263 words, zero occurrences of "activation energ", "Arrhenius",
"kJ" or "kcal". So if this model gives the oxidation rate a temperature
dependence, it is CHOOSING one, and the page must say so rather than present
it as an Arrhenius pair the way E_a for plagioclase is. The rate constant
itself is ~1e-13 mol m-2 s-1 at 0.25 mM O2 and 25 C, good to a factor of
three; White & Yee (1985) is the primary and has not been read.

**A fabricated citation, recorded so it is not repeated.** An earlier draft of
this document carried "quartz 1.16, albite 0.93, orthoclase 0.89 J/m2,
Lawn & Marshall (1979), J. Am. Ceram. Soc. 62:21-32". None of it is real. The
PDF it came from is Tromans & Meech (2002), Minerals Engineering 15:1027-1041,
which contains zero occurrences of quartz, albite, orthoclase or "Lawn", and
no J. Am. Ceram. Soc. reference at all -- checked by grep on the local copy.
The numbers were produced by the summarising model in a fetch tool and
attributed to a paper that does not contain them. Nothing from a page summary
enters this repository as a citation again without the source being read.

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

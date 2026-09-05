# 09 -- Oxygen in the water, and biotite weathering: what the sources say

Written 2026-09-05, from PDFs read directly rather than from search summaries,
because this repository has already had one fabricated citation reach a design
document (see design 08). Every number below is quoted from a table or a
sentence in the paper named. Where a source is calibrated rather than
measured, that is said.

Sources read for this: Fletcher, Buss & Brantley (2006) in full including
Table 1; Bazilevskaya et al., *Reviews in Mineralogy & Geochemistry*; Behrens
et al. (2015); Buss et al. (2010), *Chemical Geology* 269:52-61.

## 1. Oxidation is the FIRST and DEEPEST weathering reaction, and it is measured

Behrens et al. (2015), on the Hakgala profile in Sri Lanka:

> "As oxidation is the first weathering reaction, the supply of O2 is a
> rate-limiting factor for chemical weathering."

> "the dissolution front of plagioclase (consumption of protons, zone 3) lies
> higher than the O2 consumption front (zone 2)"

That is design 08's premise, measured in a profile: the oxidation front is
BELOW the plagioclase front. Not inferred, not modelled -- ordered by depth in
the field.

They also state the mechanism this model now runs on:

> "Fe(II) in bedrock minerals is oxidized after O2 transport from the soil
> surface into the rock"

> "The strength of the feedback depends on the relative weight of advective
> versus diffusive transport of O2 through the weathering profile."

## 2. WHY THIS MODEL HAS NO DEPTH PROFILE, AND WHY THAT IS RIGHT FOR GRANITE

The flip to oxidation removed the model's weathering profile with depth: at
t90 the extent runs 0.95, 0.90, 0.94, 0.89, 0.94, 0.89 down a 3 m section,
which is flat. That looked like a defect. It is not, and Bazilevskaya et al.
say so directly:

> "In rocks with low FeO concentrations, oxygen may not be consumed at shallow
> depths, and oxidation may therefore be the deepest weathering reaction."

> "the high content of ferrous iron may deplete weathering fluids in O2 at
> depth so that the deepest weathering reaction is CO2- rather than
> O2-promoted"

The controlling quantity is a ratio, and they name it:

> "In this ratio of FeO to base cation oxides, the numerator and denominator
> summarize the relative capacity of the rock to consume oxygen versus the
> capacity to consume carbon dioxide during weathering, respectively."

Granite is the low-FeO end -- `f_FeO = 0.011` here -- so oxygen is NOT consumed
over the top few metres and oxidation runs deep. A flat oxidation profile over
3 m of granite is the published expectation. A model of diabase would have to
behave differently.

**And there is a field measurement of exactly this.** Fletcher's Table 1 gives
dissolved O2 in pore fluid at the bedrock-saprolite interface at Rio Icacos,
sampled by lysimeter at 8.5 m depth:

    cR   = 2.3e-7 mol/cm3  =  0.230 mol/m3   at 8.5 m
    catm = 2.6e-7 mol/cm3  =  0.260 mol/m3   atmospheric equilibrium

**88 % of saturation at 8.5 m depth.** Oxygen really does reach the weathering
front nearly undepleted, which is what a Damkohler of 0.023 says and what the
flat profile shows.

Two caveats Fletcher states himself and which must travel with that number:

> "no attempts were made to keep atmospheric oxygen from interacting with the
> sample (White, A., pers. comm.). The value of cR, a major unknown in the
> model, could vary significantly toward very low values"

> "biotic reactions will typically use oxygen and create more reducing
> conditions at the bottom of weathering regolith as compared to surficial
> systems"

So 88 % is an upper bound. The sink this model does not have is BIOLOGICAL --
organic matter and root respiration -- and Buss et al. (2010) make it concrete:

> "microbial consumption of oxygen during respiration may affect the rate of
> spheroidal fracturing, which is thought to be controlled by the
> concentration of oxygen in the porewater"

> "The preferential mobilization of light Fe isotopes at the bedrock-saprolite
> interface suggests that microorganisms are likely influencing concentrations
> of Fe and oxygen in the porewater"

## 3. THE RATE CONSTANT: OUR PROVENANCE IS MUCH WEAKER THAN IT WAS WRITTEN

Design 08 records the oxidation rate as "~1e-13 mol m-2 s-1 at 0.25 mM O2 and
25 C, good to a factor of three", attributed to White & Yee (1985), GCA
49:1263-1275, unread. Two things are now known and neither is comfortable.

**Fletcher (2006) does not cite White & Yee.** His reference list has no such
entry; the only 1985 paper in it is Rice, Buol & Weed, *Soil Sci. Soc. Am. J.*
49:171-178, on saprolite profiles. So whatever secondary chain produced the
1e-13 attribution, it does not run through Fletcher.

**And Fletcher's own value is four thousand times larger.** Table 1:

| | Fletcher (2006) | this model | ratio |
|---|---|---|---|
| areal rate at c = 0.23 mol/m3 | 3.9e-10 mol m-2 s-1 | 9.2e-14 | **4239x** |
| reactive area per m3 of rock | 3.0e4 m2/m3 (S = 0.2 m2/g, BET-like) | 1.8e2 (geometric) | 167x |
| volumetric rate | 1.17e-5 mol m-3 s-1 | 1.66e-11 | **7.1e5x** |
| effective O2 diffusivity in corestone | 5.0e-10 m2/s | 1.45e-13 | 3455x |
| tortuosity of corestone | **3** | **1e4** | 3300x |

Seven hundred thousand fold apart on the volumetric rate. That is not a
disagreement about a parameter; it is a disagreement about what kind of rock
the oxygen is moving through.

**But Fletcher's numbers are a CALIBRATION, and he says so.** K and Theta are
back-calculated from two observations:

> "The two quantities, W and tcrack, are assumed to equal the observed value
> (Figure 2) and the value inferred from the postulate of steady-state
> denudation (R = W/tcrack = 1 cm/100 y) respectively."

W = 2.6 cm is the observed mean rindlet thickness; tcrack = 260 y follows from
assuming the fracturing front keeps pace with denudation. Equation (20) then
*yields* K = 7.8e-11 mol/(g s) and Theta = 5e-6 cm2/s. They are outputs of a
fit to a landscape, not inputs from a laboratory.

His transport parameter is strained even so, by his own admission:

> "A conventional estimate for diffusivity of aqueous solutes in water,
> D = 1e-5 cm2/s, requires an unrealistically high value of porosity of
> phi = 0.5 in the case that m = 1. ... The defensible choices of phi = 0.01
> and m = 1 require D = 5e-4 cm2/s. Given the uncertainties ... within 1.5
> orders of magnitude of the conventional value, is acceptable."

D = 5e-4 cm2/s is **23.8 times free-water O2 diffusivity**. A solute cannot
diffuse through rock faster than through water, so that parameter is absorbing
something else -- almost certainly the advection along the fractures his own
model creates.

**Where that leaves this model.** The honest statement is that the oxidation
rate constant is the WEAKEST parameter in the model, not a factor-of-three
quantity, and that the two candidate values differ by three to four orders of
magnitude with the larger one calibrated to a tropical denudation rate. The
model uses the smaller, laboratory-flavoured value paired with a GEOMETRIC
surface area, which is the same pairing as the plagioclase side and the
standard way to land near field behaviour (White & Brantley 2003). Consistency
between the two reactions is worth more here than agreement with Fletcher.

The consequence is visible and should be stated on the page: at this rate the
oxidation length is 132 m and nothing is consumed over 3 m. Were the rate a
hundred times faster it would be 1.3 m and a depth profile would appear.

## 4. The fracture energy: Fletcher chose it for fit, and says which way

Design 08 asserts that Fletcher's Gamma = 200 J/m2 comes from Friedman, Handin
& Alani (1972) and is specimen-scale. Verified in the source, with a detail
worth having:

> "the value of Gamma is obtained for tensile fracture in polycrystalline rock
> [34]; values from single crystals are smaller by a factor of ~0.01 and
> values obtained for shear fractures are larger by a factor of ~100. **The
> best model fit is obtained for the smallest value, and we adopt that here.**"

So the 200 J/m2 is the smallest of the 2-5e5 erg/cm2 range, adopted because it
fit best -- and Fletcher's own single-crystal figure is 0.01 of that, i.e.
2 J/m2. Design 08's bracket (floor 1-2, this model 7.3, specimen scale
200-500) survives contact with the source, and the "Fletcher's own
single-crystal estimate 2-5" line is confirmed as his 0.01 factor.

## 5. Rindlets: the length scales to match

Fletcher (2006), on the Rio Blanco quartz diorite:

> "The transition between the coarse-grained Rio Blanco quartz diorite and the
> saprolite is not a sharp interface, but consists rather of a **20-60 cm-thick
> zone** characterized by fracture-bound concentric shells, termed here
> rindlets"

> "results consistent with the observed thickness of rindlets in the Rio Icacos
> bedrock (**~2-3 cm**) and a time interval between fractures (**~200-300 a**)"

> "fractures occur every ~250 years, ferric oxide is fully depleted over a four
> rindlet set in ~1000 years, and saprolitization is completed in ~5000 years
> in the zone containing ~20 rindlets"

This model's O2 penetration into intact rock is 4.5 cm and its rind is 15 cm
and up. The penetration is within a factor of two of an individual rindlet;
the rind is a factor of 1.3 to 4 below the rindlet ZONE. Both are the right
order, and neither was tuned.

## 6. The fracture feedback decides the answer, and there is a natural experiment

Bazilevskaya et al., on paired sites of contrasting lithology:

> "both the regolith itself and the plagioclase weathering front are 20-times
> thicker on the metagranite ... than on a nearby diabase"

> "the greater thickness of the regolith and the reaction front on the
> metagranite compared to the diabase was due to differences in solute
> transport during early weathering: WIF allowed solute transport by advection
> that in turn led to thicker regolith. In contrast, on the nearby diabase ...
> microfracturing was not observed and solute transport was limited to
> diffusive processes."

Weathering-induced fracturing changing transport from diffusive to advective
is worth a factor of twenty in regolith thickness, in the field. Probe I put
the same feedback at a factor of about thirty in this model. That is an
agreement between an independent field comparison and a scaling estimate, and
it is the strongest argument yet that design 08 step 3 -- the cracking
criterion -- decides the answer rather than decorating it.

## What to do about it

1. **Correct the rate constant's provenance in the code.** It is written as
   "good to a factor of three". It is not.
2. **Say on the page that the oxidation rate is the weakest link**, and that a
   hundredfold faster rate would put a depth profile back.
3. **Do not adopt Fletcher's parameters to close the gap.** They are a fit to a
   tropical denudation rate, and his diffusivity is 24x free water.
4. Still worth an interlibrary request: White & Yee (1985), to find out whether
   the 1e-13 attribution has a source at all.

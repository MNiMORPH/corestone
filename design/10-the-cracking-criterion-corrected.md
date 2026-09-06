# 10 -- The cracking criterion: Goodfellow's, not Fletcher's

Written 2026-09-06, before implementing design 08 step 3, in answer to one
question: *is this strain energy, or more?*

**It is strain energy. It is also, and equivalently, a stress criterion. And
it is more than both, in one way that is measured and that design 08 misses.**

Read this instead of design 08's "The mechanical feedback, without fracture
mechanics" section, which is superseded in four places.

## 1. Goodfellow et al. (2016) already did this calculation, for this rock

Design 08 proposed building Fletcher, Buss & Brantley's (2006) energy budget.
Goodfellow et al. (2016) -- the paper this model already takes `k_matrix`,
`k_weathered` and the oxidation mechanism from -- ran the same criterion on
the same kind of rock, and reported the answer. Their section 4.2:

> "In unweathered rock, and using a high Gamma (2 x 10^3 ergs cm^-2), L is at
> the meter scale ... after about 20% of the original Fe(II) in biotite has
> been oxidized ... L has decreased to the centimeter scale and reaches the
> millimeter scale once about 50-60% of the original Fe(II) has oxidized."

> "matrix fracturing is initiated over a range of **20-65% oxidation of
> biotite Fe** in our samples. This process will occur sooner after the onset
> of weathering in more Fe-rich granitoids, and later in Fe-poor granitoids."

Their criterion is a LENGTH: `L*`, "the distance over which the rock supports
fracturing", from their equation 13, and matrix fracturing occurs when `L*` is
at or below the crystal length. That is algebraically design 08's
`U >= 2 Gamma / d` rearranged -- `L* = 2 Gamma / U` -- so the two documents
propose the same criterion, and only one of them has run it.

## 2. FOUR NUMBERS IN DESIGN 08 ARE WRONG, AND ALL FOUR COME FROM GOODFELLOW

| quantity | design 08 | Goodfellow et al. (2016), measured or chosen |
|---|---|---|
| `x_c`, oxidation at cracking | **0.10**, "from Goodfellow's ~10 % Fe(III) rise" | **0.20 to 0.65**, their own stated result |
| `Gamma` | 7.3 J/m2 implied | **0.2 to 2 J/m2** (2e2 to 2e3 erg/cm2), the range they use |
| `dV/V` | 0.735, FeO -> goethite, MINERAL level, applied to `f_FeO` | **0.04 (ferrihydrite) or 0.05 (altered biotite)**, GRAIN level, applied to the biotite fraction |
| product | goethite | **both**: altered biotite AND ferrihydrite, run as separate scenarios |

The `dV/V` row is the one to be careful about, because the two numbers are not
the same quantity and neither is wrong on its own terms. Fletcher's 0.7 is the
volume ratio of the FeO component to its oxidation product, so it multiplies
`f_FeO`. Goodfellow's 0.05 is the expansion of the biotite CRYSTAL -- the
d-spacing going from 10 A to 10.5 A -- so it multiplies `phi_biotite`. Bulk
volumetric strain at full oxidation:

    design 08   f_FeO * 0.735    = 0.011 * 0.735 = 0.0081
    Goodfellow  phi_bt * 0.05    = 0.06  * 0.05  = 0.0030

Use Goodfellow's. It is measured on the rock rather than computed from a
thermodynamic compilation, it is the one their own 20-65 % rests on, and it
comes with its second pathway (ferrihydrite, 0.04) already quantified.

`x_c = 0.10` should not have been called a calibration "from Goodfellow". They
say 20-65 %. The 10 % appears to be their Fe(III) rise accompanying the
conductivity jump, which is a different observation from the fracture
threshold.

## 3. STRAIN ENERGY AND STRESS ARE THE SAME CRITERION, AND THAT IS FREE CHECK

    stress:  sigma_t     = E eps / (1 - nu)
    energy:  2 Gamma / d = E eps^2 / (1 - nu)

Eliminate `eps` and the two coincide exactly when

    Gamma = sigma_t^2 d (1 - nu) / (2 E)

which is Griffith with the GRAIN as the flaw. So there is no choice to make
between them; there is one criterion with two dials, and requiring them to
agree pins the pair.

**And they do agree, on Goodfellow's own numbers, independently.** They CHOSE
Gamma in 0.2 to 2 J/m2 on a crack-tip argument. They separately MEASURED the
tensile strength of their least weathered granodiorite at 6.3 MPa. Through the
relation above, at `d` = 2 mm with SKB's E = 76 GPa and nu = 0.24:

    sigma_t 6.3 MPa  (Goodfellow, least weathered)   ->  Gamma = 0.40 J/m2
    sigma_t 3.0 MPa  (Goodfellow, cobble group 3)    ->  Gamma = 0.09
    sigma_t 1.0 MPa  (Goodfellow, ITRT)              ->  Gamma = 0.01
    sigma_t 13.5 MPa (SKB R-05-83, Forsmark granite) ->  Gamma = 1.82

0.40 sits inside the 0.2-2 they chose, and 1.82 does too. Two independent
routes to the same order, and it also shows why Fletcher's 200 J/m2 cannot be
right at this scale: it would require a tensile strength of 140 MPa.

**Report both.** The model should print the energy ratio AND the equivalent
stress, so a reader sees they are one statement, and so a wrong parameter
shows up as the two disagreeing.

## 4. WHERE IT IS MORE THAN STRAIN ENERGY: THE STRENGTH IS NOT A CONSTANT

This is the part neither design 08 nor a fixed threshold can represent, and
Goodfellow measured it (their section 4.2.2, Figure 7, Table S8):

    least weathered cobbles (groups 1, 2)   6.3 MPa
    cobble group 3                          3.0
    cobble group 4                          2.8
    ITRB                                    2.3
    cobble group 5                          1.6
    ITRT                                    1.0
    saprolite                               0.1-0.2

> "Tensile strength therefore appears **unaffected by slight
> weathering-related oxidation** of biotite Fe."

> "The total decrease in tensile strength from unweathered granodiorite to
> saprolite exceeds 6 MPa, with the **largest decline (exceeding 50%)
> occurring early** in weathering, from the parent rock to cobble group 3."

A factor of sixty, flat at first and then collapsing. So the criterion is not
"reach a fixed threshold and crack": the threshold FALLS as the rock cracks,
which is a runaway, and it is the mechanical half of the same positive
feedback that `k(M)` and `tortuosity(M)` already carry on the transport side.

That the strength is flat through *slight* oxidation and then drops sharply is
also the cleanest evidence there is for a threshold existing at all.

## 5. What it still leaves out, stated rather than hidden

- **Time.** Subcritical crack growth means rock under stress below the
  instantaneous threshold still cracks, given long enough, and this model's
  steps are kiloyears. An energy balance has no clock.
- **The process zone.** `U_c = 2 Gamma / d` counts only new surface. Zang et
  al. (2000) measured the process zone in Aue granite at 2-9 grain diameters,
  so real dissipation exceeds surface creation, and the effective `Gamma` at
  crack scale is above the thermodynamic floor. This is why the bracket
  cannot be closed from below.
- **Porosity relief.** Some of the expansion goes into pore space rather than
  straining the matrix. Granite porosity is ~1 % and Goodfellow measure it
  RISING with weathering, so this relief grows exactly as the strain does.
  Ignoring it overestimates the strain, and by an amount nobody has measured.

## 6. What to build

1. `bulk_volumetric_strain(x)` from `phi_biotite` and a grain expansion that
   defaults to Goodfellow's 0.05, with 0.04 available as the ferrihydrite
   scenario.
2. `tensile_strength(x)`, interpolating Goodfellow's measured sequence, so the
   strength falls as the rock oxidises.
3. `cracking_number(x)` -- report BOTH forms, and assert they agree through
   the Griffith relation.
4. Report `x_c` as a BRACKET and check it against Goodfellow's 20-65 %, rather
   than calibrating to a single number.
5. Do NOT calibrate. Design 08 calibrated `x_c` because Gamma spanned 200-fold;
   with Goodfellow's own Gamma and their own tensile strength the spread is a
   factor of ten, and a prediction inside their stated 20-65 % is worth more
   than a fit to it.

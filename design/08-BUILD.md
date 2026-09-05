# Build plan for design 08 -- self-contained

Written 2026-09-05. **Assume the reader has no memory of the session that
produced it.** Everything needed to execute is here or in
`design/08-oxidation-drives-it.md`; nothing depends on a conversation.

Read `design/08-oxidation-drives-it.md` first for *why*. This file is *how*,
in order, with what verifies each step.

---

## 0. Not blocked. What was recorded here as a blocker, and why it is not one

This section used to say "BLOCKED until Andy answers two questions" and to
forbid step 3. Both questions were raised in the design document by its own
author and neither came from Andy; on 2026-09-05 he read them and the block
was lifted. They are still open, and they are page decisions rather than code
decisions:

- **Does the colour bar become oxidation extent?** The demo's right panel is
  labelled "soluble phase dissolved". After the flip the same field means how
  much of the iron has rusted, not how much mineral has left the rock. A
  caption and a paragraph, best decided with a rendered figure in hand rather
  than in the abstract -- which is why it should never have stopped the code.
- **Does the page keep a sentence on plagioclase dissolution**, the mass loss
  this model will not track? A note already exists on the page ("One mineral,
  where a rock has many") and may only need rewording.

**The decision that IS load-bearing was found later and is not either of
these.** Probe I (`prototypes/probe_i_oxygen_regime.py`) measured the
section-scale Damkohler at 0.0227 on oxygen against 6.56 on silica, so the
model changes limit: from saturation-limited to reaction-limited. The page's
central claim -- a corestone survives because the water arrived spent --
becomes false, and is replaced by a corestone surviving because O2 cannot
diffuse in over 4.5 cm. Both are "the water never got there"; they are not the
same lesson. Read probe I before writing step 2.

---

## 1. The parameters, with sources. Do not re-derive these.

| symbol | value | where it came from |
|---|---|---|
| `f_FeO` | **0.011** | granite, volume fraction of FeO component (all Fe as FeO). USGS reference granites G-1/G-2/G-3 give 0.0083-0.0116; Goodfellow's granodiorite 0.0107. **Fletcher's 0.05 is wrong** -- it implies 10.9 wt% FeO and his own rock measures 4.6 |
| `V_FeO` | 12.00 cm3/mol | Robie & Hemingway (1995), USGS Bulletin 2131, p. 16 |
| `V_goethite` | 20.82 cm3/mol Fe | same table. **Goethite, not ferrihydrite** -- that compilation has no ferrihydrite entry, and Fletcher's quoted 0.7 is only reproducible as FeO -> goethite |
| `dV/V` | **0.735** | the two above; both single-Fe formula units |
| `E`, `nu` | 76 GPa, 0.24 | SKB R-05-83 Tables 2-7/2-9, 52 intact granite cores. Use 40 GPa if you want near-surface cracked rock; it moves `x_c` by 1.6x |
| `d` | 2 mm | already `self.grain_size`. Sets reactive surface area AND dispersivity AND this |
| `x_c` | **0.10** | **CALIBRATION**, from Goodfellow's ~10 % Fe(III) rise. Implies Gamma = 7.3 J/m2, which is checkable and lands between the 1-2 floor and the 200-500 specimen-scale ceiling |
| `C_O2` | 0.338 mol/m3 | dissolved O2 at 12 C (10.8 mg/L). Tabulated; 0.456 at 0 C, 0.237 at 30 C |
| `tau_O2` | ~679 | `0.25 * f_FeO / (V_FeO * C_O2)`. Compare the current 47744 for silica: oxygen is 70x less limiting, so the brake becomes diffusion, not budget |
| rate | ~1e-13 mol m-2 s-1 | at 0.25 mM O2, 25 C, **good to a factor of 3 only**. White & Yee (1985) is the unread primary; two secondary renderings disagree 1.5-2.4x and reverse a rank order |
| `E_a` | **does not exist** | no measured activation energy for aqueous Fe(II)-silicate oxidation by O2. Verified: Hogg & Meads (1975), a dedicated Mossbauer kinetics study, has zero occurrences of "activation energ", "Arrhenius", "kJ" or "kcal" in 4263 words. **If you give the rate a temperature dependence you are CHOOSING one -- say so on the page, do not present it as an Arrhenius pair** |
| O2 order | **1, not 0.25** | Fletcher's C^0.25 is the stoichiometric quarter borrowed as a concentration exponent, with no experiment behind it. Lebedeva & Brantley (2020) eq. 15, same group and reaction, is first order |

## 2. The structural change -- a re-pointing, not a rewrite

The solute flips from **product** to **reactant**. Concretely, in
`solve_solute` and `_transport_operator`:

    now:  div(q c) - div(D grad c) + r c = r        inlet c = 0, c -> 1
    then: div(q c) - div(D grad c) + r M c = 0      inlet c = 1, c -> 0

Same operator, same sparse assembly, same cached factorisation. Only the
right-hand side (goes to zero) and the inlet boundary (becomes 1) change.

The rock update keeps its exact exponential integrator: the rate becomes
proportional to `c M` instead of `(1 - c) M`, so `M(t+dt) = M exp(-lam dt)`
still holds with `c` held, only `lam ∝ c` rather than `∝ (1 - c)`.

`k(M)` and `link_tortuosity(M)` need **no change at all**. That is the point:
`M` becomes unoxidised Fe(II), which is Goodfellow's own independent variable,
so the interpolations finally index on the thing their data was measured
against.

## 3. The cracking criterion

    eps(x) = x * f_FeO * (dV/V) / 3        LINEAR strain, x = 1 - M
    U(x)   = E eps^2 / (1 - nu)            Fletcher eq. 11, laterally confined
    U_c    = 2 Gamma / d
    N(x)   = U / U_c

**Not `0.5 K eps^2`.** The oxidising rind is confined by the rock around it,
so the strain parallel to the front is zero and the modulus is `E/(1-nu)`, not
the bulk modulus. The `/3` converts volumetric to linear strain.

Implement by calibration: set `x_c = 0.10` and derive `Gamma` from it, rather
than setting `Gamma` and predicting `x_c`. Gamma spans 200-fold in the
literature and will not close; `x_c` is one observation. Report the implied
Gamma as a checkable output -- 7.3 J/m2.

## 4. Order of work, with the check for each

1. **DONE 2026-09-05.** Parameters and the O2 solubility function. Check: `thermo_report()`
   prints `tau_O2` near 679 at 12 C and the oxygen ceiling on front advance
   near 442 m/Myr.
2. **Flip the solute.** Check: with `M = 1` everywhere, `c` is 1 at the
   surface and falls with depth; the cell balance test in
   `test_stated_equations.py` must be re-transcribed to the new RHS and must
   still pass per cell.
3. **Cracking criterion + `x_c`.** Check: a test that `N(x_c) == 1` and that
   the implied Gamma is 7.3 J/m2 to two significant figures.
4. **Register every new displayed equation** in `tests/test_equation_coverage.py`
   against a test that verifies it, and bump `EXPECTED_BLOCKS` /
   `EXPECTED_DISTINCT`. The suite fails otherwise, by design.
5. **Prove one test bites** by reverting the change it guards, per the repo
   rule. The obvious one: pin `c` to its old sense and watch the depth
   profile invert.
6. **Re-measure everything.** `t50/t90/t99` for the six settings, frame costs
   at three cell sizes and two temperatures, watch times, the front rate.
   Scripts to copy: the ones described in FRAME (g). Then update, in this
   order: the demo comments in `interactive_demo/corestone_panel.py`, the
   module docstring's validation section, then the exercise page.
7. **Rebuild with `--strip-vendored`.** Without it the shared panel wheel goes
   9.4 MB -> 30.3 and every reader pays 21 MB. There is now a warning, but do
   not rely on reading it.

## 5. Do not

- **Do not use Fletcher's `f_FeO = 0.05`, `C^0.25`, or `Gamma = 200 J/m2`.**
  All three are wrong for this model and each is wrong for a different reason;
  see the table and design 08.
- **Do not claim the criterion predicts Goodfellow's 10 %.** An earlier draft
  did. It was a coincidence of two compensating errors and is withdrawn.
- **Do not add a temperature dependence to the oxidation rate silently.**
  There is no measured activation energy. Choosing one is allowed; hiding the
  choice is not.
- **Do not re-walk the two dead ends in FRAME (f)** -- reusing the solute LU
  across a rebuild is 2.1x slower, and raising Show's step cap moves the
  weathering front by 0.22.

## 6. What "done" looks like

The model runs on O2 and biotite, the rind still sharpens, the demo's numbers
and the page agree with the model, the equation ledger is complete, and the
page says plainly that the oxidation rate's temperature dependence is a choice
and that `x_c` is a calibration. The front rate will move; report it, do not
tune to it.

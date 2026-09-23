# 12 -- T_ref is doing two jobs, and the constants belong to the other one

Written 2026-09-23. **This is a defect, not a decision.** It is recorded here
before any fix because fixing it moves every timescale in the model and on the
teaching page, and because which of the two repairs to make is a judgement.

## 1. The finding

`T_ref = 285.0` K (11.85 C) is declared as *"reference temperature"*. Both
temperature corrections are taken against it:

    rate_factor       = k(T) / k(T_ref)          exp(-(E_a    /R)(1/T - 1/T_ref))
    solubility_factor = C_eq(T) / C_eq(T_ref)    exp(-(dH_r   /R)(1/T - 1/T_ref))

so at `T = T_ref` both are 1, and `L_ref` and `pore_volumes_ref` are *by
definition* the values at `T_ref`.

**They are not.** Both were computed from 25 C literature:

    L_ref = 0.457 m  =  q C_eq / (k A)  with
        C_eq  quartz saturation, llnl.dat AT 25 C
        k     10^-11.84 mol m-2 s-1, Palandri & Kharaka AT 25 C
    pore_volumes_ref = 47744 = N_0 / C_eq, same C_eq AT 25 C

Verified by recomputation: `q C_eq/(kA) = 0.4567 m` from the stated inputs,
matching the 0.457 in the file. The commit that introduced it
(`cf7ec6e`, 2026-09-04, "Derive the saturation length from grain size") does
not mention temperature anywhere; `T_ref = 285` predates it.

So the model applies NO correction at its default temperature, and therefore
runs 25 C chemistry while reporting 12 C.

## 2. Why this reads as an oversight rather than a convention

The same file gets this right twice, and in both cases by naming the
temperature *at which the data were measured*:

    T_K_ref = 293.15    "temperature of the measured conductivities"
    T_D_ref = 298.15    "silica [m2/s] AT this temperature"
    T_ref   = 285.0     "reference temperature"          <- the odd one out

And `self.T = self.T_ref` at line 904. **That is the root of it: one variable
is serving as both the thermodynamic reference and the default ambient
temperature.** Those are different numbers. The data are at 25 C; a temperate
default is 12 C. Conflating them silences the correction exactly at the
default.

## 3. What it costs

At 285 K relative to 298.15 K: `k` is 0.2727, `C_eq` is 0.5421.

| quantity | as shipped (really 25 C) | true 12 C value |
|---|---|---|
| `L_ref` | 0.457 m | **0.908 m** |
| `pore_volumes_ref` | 47 744 | **88 073** |
| supply ceiling `q C_eq / N_0` | 6.284 m/Myr | **3.406 m/Myr** |

So the correction makes the model **slower**, by about 1.8x on the supply-limited
arm. The model is already too slow against field fronts, so this **widens** the
discrepancy rather than closing it. See design 11 section 2a: the gap was
4.9-8.6x on the page's own number, and my own t90 measurement makes it larger
still.

## 4. The two repairs

**A. Move the reference to the data, and split the default out.**
`T_ref = 298.15`, keep `L_ref = 0.457` and `pore_volumes_ref = 47744`, and give
the default operating temperature its own name so `self.T` no longer inherits
the reference. Conceptually right -- a reference state should sit where the
measurements are, which is how `T_K_ref` and `T_D_ref` already work -- and it
leaves the two sourced constants untouched and still checkable against P&K and
llnl.dat. Costs a second attribute and one line at `self.T = ...`.

**B. Keep `T_ref = 285` and convert the constants.**
`L_ref = 0.908`, `pore_volumes_ref = 88073`. No structural change, but the two
headline constants no longer match the sources they cite, so every future
reader must redo this conversion to check them against the literature. It also
leaves the one-variable-two-jobs conflation in place to be tripped over again.

**Recommended: A.** It fixes the cause rather than the symptom, and it matches
what the file already does twice.

## 5. Blast radius, whichever is chosen

- Every timescale in the module docstring, `pacing.txt`, and the teaching page.
- The demo's pacing table must be re-measured.
- Design 11's comparison table, which prices the ceiling against a baseline
  that is currently ~1.8x too fast.
- The page's validation section, which is already in question for a separate
  reason (see design 13 if written: `D/t90` is not a front rate).

**Not decided. Andy's call.**

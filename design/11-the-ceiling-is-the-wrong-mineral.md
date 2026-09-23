# 11 -- The ceiling is the wrong mineral's

Written 2026-09-23, before any code change, in answer to Andy's remark that
*"quartz/silica should always be saturated in water, I think"*.

**He is right, and the consequence is worse than the sloppy sentence it started
from. `C_eq` is quartz saturation. Natural waters run 2-4x ABOVE quartz
saturation. So `1 - C/C_eq` is not approximately zero at the far end of a flow
path -- it is NEGATIVE, everywhere, in the real system the model is about.**

Nothing is changed yet. `C_eq` is calibrated and every timescale hangs off it,
so the decision is Andy's and it is recorded here rather than acted on.

## 1. What the model does now

`weathering.py` sets

    self.delta_H_r = 32.9e3    # quartz, llnl.dat at 25 C
    C_eq: quartz saturation, 1.0e-4 mol/kg, i.e. 0.10 mol Si/m3

and drives dissolution with `1 - C/C_eq`, which stops the reaction when the
water reaches quartz saturation.

## 2. The arithmetic

Recomputed in session on 2026-09-23 from the concentrations below, not relayed:

| water | SiO2 | as Si | Omega vs quartz | `1 - Omega` |
|---|---|---|---|---|
| Hem median, surface water | 14 mg/L | 233 uM | 2.33 | **-1.33** |
| Hem median, groundwater | 17 mg/L | 283 uM | 2.83 | **-1.83** |
| granitic groundwater, typical | 25 mg/L | 416 uM | 4.16 | **-3.16** |
| amorphous silica (upper limit) | 116 mg/L | 1931 uM | 19.3 | -18.3 |

Quartz solubility at 25 C is 6.0 mg/L SiO2 = 100 uM, which is the 1.0e-4
mol/kg already in the file. The driving force does not approach zero; it goes
through it and out the other side. In most implementations that clamps
weathering to zero or flips its sign.

## 3. Why plagioclase weathers anyway

This is the part that resolves the apparent paradox, and it is the same
distinction as design 10's `delta_H_r` argument arriving from the other side.

**Plagioclase equilibrium is not a silica-only condition.** Its ion activity
product carries Na+, Ca2+, Al3+ and H+ as well as SiO2(aq). Al3+ is pinned near
1e-7 M by kaolinite/gibbsite precipitation, and H+ is resupplied by soil CO2, so
Q stays orders of magnitude below K no matter what silica does. Water can sit at
four times quartz saturation and still be wildly undersaturated with respect to
the feldspar that is dissolving in it. White & Brantley (2003) put it directly:
*"column effluents remained far from thermodynamic saturation with respect to
plagioclase."*

**What actually caps silica** is secondary-clay (kaolinite) buffering, with
amorphous silica as the absolute upper limit, modulated by flushing. Not quartz.
Hem (1985, p. 70): *"The direct precipitation of quartz is unlikely to control
solubility of silica in most natural waters at Earth surface temperatures."*

**The standard rate law writes the affinity with respect to the DISSOLVING
mineral**, never with respect to quartz -- Palandri & Kharaka (2004) Eq. 4-7,
`g(dG_r) = (1 - Omega^p)^q` with `Omega = Q/K` for that mineral.

## 4. Why quartz survives weathering, correctly weighted

Both kinetics and thermodynamics, but **thermodynamics is the larger term**, and
I had this backwards when I first said it out loud. From Palandri & Kharaka
(2004) Table 13, quartz sits only ~1.5-2.4 orders of magnitude below
albite/oligoclase -- not the five or six usually implied. The bigger reason
quartz grains are what grus is made of is that natural water has no
thermodynamic drive to dissolve them at all.

## 5. The options, with their costs

| option | what it costs |
|---|---|
| **Leave it.** Treat `C_eq` as a teaching abstraction: "the ceiling", never identified with a real mineral. | Free. But the page names quartz, and 32.9 kJ/mol IS quartz's, so the identification is already made in print. |
| **Move `C_eq` to amorphous silica** (1931 uM). | `pore_volumes` 47 740 -> 2 473, so **19.3x faster weathering**. Every timescale on the page and in `pacing.txt` is invalidated. The demo's pacing table must be re-measured. |
| **Move `C_eq` to kaolinite-buffered silica.** | Most defensible physically. Needs an activity diagram and a stated pH/Al assumption -- i.e. a thermodynamics where the model currently has a normalisation. Same timescale blast radius. |
| **Keep quartz, and say plainly that it is a proxy** for whatever stops the reaction. | Free, honest, and consistent with design 10 section on what `delta_H_r` belongs to. Does not fix the sign problem, only stops claiming it is not there. |

**Not decided. Andy's call.** The blast radius on options 2 and 3 is the reason
it was not done on the spot: they are not parameter edits, they are a
recalibration.

## 6. What was actually done on 2026-09-23

Subtraction only, and Andy confirmed that was the right instinct. The page had
claimed:

> *"quartz does not dissolve in this model for a reason rather than by decree:
> the water is already at its saturation, so its driving force `(1 - C/C_eq)` is
> zero"*

That fails where the model is most active -- rain enters the joints at `C = 0`,
maximally undersaturated with respect to quartz, so on that argument quartz would
dissolve fastest exactly there. Cut to "quartz itself is held inert here", which
is what the code always said. **No new claim was put in its place**, because the
claim worth making is that the ceiling is the wrong mineral's, and that is a
model change nobody has made.

## 7. Provenance, including the gaps

**Read directly this session:** Hem (1985, USGS WSP 2254); Palandri & Kharaka
(2004, USGS OFR 2004-1068, Tables 4 and 13); White & Brantley (2003); Maher et
al. (2009, Table 1); LLNL `llnl.dat` and USGS `wateq4f.dat`; Jacks (1978, KBS TR
88); Hynek et al. (2022, Table 2).

**NOT read, and any claim attributed to them here is secondary:** Drever;
Langmuir; Stumm & Morgan; Appelo & Postma; Garrels & Christ (1964); Feth et al.
(1964); Davis (1964); Morey et al. (1962, 1964).

**A live units trap, unresolved.** Rimstidt (1997) is quoted elsewhere as "about
11 mg/L" in **mg Si/L**, against Hem's 6 mg/L as **SiO2** (= 2.8 mg/L as Si). One
of the two carries a units error and neither primary source was reached. **That
number must not go near the teaching page until someone reads the original.**

**A process note worth keeping.** The first automated read of the Palandri &
Kharaka PDF returned albite as `-10.16 / -12.56 / -13.71` and oligoclase as
`-9.74 / -12.26`. Both were wrong; Table 13 actually reads albite
`-10.16 / -12.56 / -15.60` and oligoclase `-9.67 / -11.84`. It was caught only by
extracting the PDF text and reading the table. Summarised numbers from a PDF are
not measurements.

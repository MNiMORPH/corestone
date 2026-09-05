# Frame -- read this first

The read-first frame for this repo, per `~/.claude/COMPACTION_PLAYBOOK.md`.
After a compaction, read this **before** acting, and verify every structural
claim below against git and disk before trusting it. Current to `f3bd9c4`.

## (a) Origin -- why this model exists

To teach one idea: weathering is a race between how fast water delivers fresh
fluid and how fast rock dissolves into it. The payload is the misconception it
corrects -- **a corestone is not tougher rock**. Same granite, same minerals,
same temperature; it survives because the water never reached it, or reached it
already saturated.

Delivery is an interactive browser demo built with
[artesian](https://github.com/MNiMORPH/artesian), which constrains the model:
it must run under Pyodide, so numpy and scipy only, nothing compiled, and fast
enough to press Run. **Verified, not assumed:** `artesian check numpy scipy
matplotlib` reports all three bundled by Pyodide.

## (b) Plan and trajectory -- as the next action

0. **Two decisions wait on Andy.** Does the colour bar become oxidation extent
   once design 08 lands -- that changes what the exercise is *about*. And the
   "What to do" section on both exercise pages is still `*(To be written.)*`,
   which is the only thing between this and a usable assignment.

1. **Next: implement design 08. `design/08-BUILD.md` is the ordered,
   self-contained plan -- start there.** Rationale in
   `design/08-oxidation-drives-it.md`. Designed in full,
   corrected against sources, not one line written. Biotite Fe(II) oxidation
   by dissolved O2 replaces plagioclase dissolution as the driver; `M` becomes
   the unoxidised fraction, so the existing `k(M)` and `tortuosity(M)`
   interpolations index on Goodfellow's own variable. The solute flips from
   product to reactant: same operator, right-hand side goes to zero, inlet
   boundary c = 1, rate proportional to `c M` instead of `(1-c) M`. The
   exponential integrator survives unchanged. Cracking is an elastic-energy
   budget, four lines, with `x_c = 0.10` calibrated on Goodfellow and the
   implied fracture energy of 7.3 J/m2 reported as the checkable by-product.

2. **Re-measure designs 02-06.** Still open and now much worse: 2026-09-04
   changed tau, the saturation length, both matrix transport terms, the joint
   conductivity, the diffusivity and the pace. Nearly every number in those
   documents is stale.

3. Done 2026-09-05: the head solve is warm-started (`e141329`) and
   `flow_tolerance` is converged to 0.02 (`cc4098d`). Flow factorisations fell
   from 850 to 21 over 2000 kyr, which is what made the tolerance affordable.
   Error against a 0.005 reference: 0.381 at the old 0.05, 0.141 now, 0.045 at
   0.01 -- which was rejected because it puts the warm end at 61 ms against a
   33 ms frame budget. **0.02 leaves the warm end ON the budget** (32 and
   35 ms in two measurements), so 30 C drops the occasional frame. Tolerable
   only because artesian's animator now yields every frame, so a dropped frame
   stretches the run instead of freezing the controls.

4. `P21` carries a systematic +8 % bias from counting wall joints at full
   length; `test_a_full_orthogonal_grid_has_the_analytic_intensity` hides it
   behind `rel=0.15`. Untouched.

5. **Three interlibrary requests would close real uncertainties.** White & Yee
   (1985), GCA 49:1263-1275 -- the primary behind the only biotite oxidation
   rate, whose two secondary renderings disagree by 1.5-2.4x and reverse a
   rank order. Zang et al. (2000), JGR 105:23651-23661, and its companion
   Janssen et al. (2001), Int. J. Earth Sci. 90:46-59 -- these would turn the
   fracture-energy bracket into a number.

6. Done and struck from this list: `f_inert` removed (`43895c9`); `E_a` and
   `delta_H_r` sourced (`bdaab2f`); the demo built and deployed.

## (c) Key current data and objects

Branch `master`, **everything committed and pushed** as of 2026-09-05 across all
three repositories (corestone, artesian, GeomorphOnline.github.io). Working
trees clean.

Live: <https://geomorphonline.github.io/exercises/corestone-weathering/>

Note that **other Claude sessions push to the site repository while this one
works**. A push was rejected on 2026-09-04 and rebased cleanly; check
`origin/master` before assuming a push landed.

## (d) Guardrails and irreversibility state

`github.com/MNiMORPH/corestone` is **public**; pushing publishes. 17 commits
are local only. Tags, releases, version bumps and closing issues each need
explicit authorisation in the current message.

Two pushes exist in the reflog: `0a57fcd` (mine, authorised) and `18cb016`
(not from my session -- Andy's).

## (e) Results, each with the method that verified it

**AS OF 2026-09-04 THE MODEL IS FULLY PARAMETERISED FROM SOURCES.** Nothing in
the chemistry or the flow is fitted: rate constant and activation energy from
Palandri & Kharaka for oligoclase; solute ceiling and enthalpy from quartz
saturation; matrix conductivities from Goodfellow et al. (2016); joint
conductivity from a 100 um hydraulic aperture through the cubic law
(Witherspoon et al. 1980); silica diffusivity from Rebreanu et al. (2008)
scaled by Stokes-Einstein; tau and the saturation length from the mineralogy
and a 2 mm grain size. Which makes the weathering timescale a PREDICTION:
**0.79 m/Myr against 4-7 measured** for temperate granite regoliths, so five
to nine times slow, reported rather than tuned. The gap sits in the reactive
surface area, geometric at 900 m2/m3 against a BET 3e5-3e6.

The rind is 4.6 % of the section part-dissolved, about 14 cm, against 30.6 %
before the matrix transport was corrected. Real rindlet zones run 20-60 cm.

Watching to 90 % dissolved: 342 s at 0 C, 126 s at 12, 35 s at 30, at 1 kyr a
frame. **That is a slower demo than it was, and better science.** Whether the
trade is right for a teaching page is a live question, not a settled one.

**WARNING -- designs 02 to 06 contain numbers measured before three corrections
and are not to be trusted without re-measurement.** In order: `a505892` made the
saturation length scale with flux; `c0d7749` added diffusion and the `C_eq`
temperature term and replaced the whole transport operator; `7cbd0a7` made
non-axis-aligned joints conduct at all. Any figure or number in a design
document predating those is stale. The prose reasoning survives; the numbers do
not.

Current and verified:

- Solute is conserved to 1e-15 at every domain size and resolution checked
  (`test_what_the_rock_loses_is_what_the_water_carries_out_of_the_base`).
- Water is conserved to ~1e-10 per cell, checked per cell rather than globally,
  because a global balance telescopes and is blind to interior error.
- Diffusion takes the fraction of the domain that is undersaturated from 8.7 %
  to 100 %, and corners sit 41 % further from saturation than faces at equal
  distance -- the rounding, with no oxidation and no fracture mechanics.
- Temperature acts through **solubility**, not the rate constant: 275->315 K
  moves mean dissolved fraction 189 %, monotone. The model is transport-limited
  almost everywhere, and there the amount dissolved scales with `C_eq`.
- Snapped orientations leave the seam within 13 % of the interior.

**Speed, 2026-09-02.** Measured by alternating the old and new trees as
subprocesses, min of three, so machine load hits both equally -- the first
attempt at this ran the variants sequentially and produced a 3.2x that was
load drift and not the solver. Mean `M` agrees to twelve decimals throughout;
none of it changes an answer.

    case      before    after     ratio
    3 m app    0.415 s   0.294 s   x1.41
    3 m dx.02 16.318 s   6.833 s   x2.39
    45 deg     0.121 s   0.074 s   x1.64
    12 x 9 m   2.726 s   1.434 s   x1.90

Three changes, in order of what they bought: minimum-degree ordering on
`A + A.T` instead of SuperLU's COLAMD default (fill halved, `565e0e4`); the
warm start (`8395d97`); the base boundary and the reaction term both assembled
without a format round trip (`ed30aca`, `739754b`).

**The integrator and the step control, 2026-09-02.** Same method as above;
`c_drift_max` at its default 0.03.

    case      before    after     ratio   steps
    3 m app    0.511 s   0.215 s   x2.38   107 -> 79
    45 deg     0.117 s   0.060 s   x1.95    33 -> 26
    12 x 9 m   3.477 s   1.276 s   x2.72    57 -> 30
    3 m dx.02 21.829 s   6.309 s   x3.46   108 -> 82

- **`M(t+dt) = M(t) exp(-lambda dt)`** replaces forward Euler (`78eb750`). `r`
  is proportional to `M`, so with `c` held the step integrates exactly. 1.1x
  to 8.4x more accurate at identical step counts, for one `np.exp`.
- **`c_drift_max` replaces `dx_max` as the control** (`ebf8f65`). It bounds the
  model's one time-step approximation -- `c` held while the rock moves --
  rather than a proxy for it, and the error is monotone in it and nearly first
  order. It is **a chosen error budget**, ~3 % of full scale, one line to
  change; the table is in `update`.
- The budget is **enforced by rejection**, and the step is **predicted** from
  the last drift (`f3bd9c4`). Predicting matters as much as enforcing: reaching
  for the growth cap instead was rejected on 73 of 75 steps.

**Still awaiting a decision: the refactorisation cap.**
`max_krylov_iterations = 15` has never fired; lowering it to 5 so that it does
is worth ~15 %. Left alone because the optimum moves between 4 and 5 case to
case. Table in the `solve_solute` docstring.

## (f) Negative results, and why not

- **No existing fracture-network generator was usable.** Eight examined
  (design 03). fractopo analyses but does not generate; FracSim2D is Python 2
  with Windows-only extensions; ADFNE is MATLAB; HatchFrac C++; dfnWorks 3D.
  Pyodide rules out compiled code regardless.
- **Free segments from a power-law length distribution do not make a network.**
- **Connectivity is the wrong measure for corestones.** A joint is a conduit;
  rock is continuous across it. The measure is distance from the network.
- **Full saturation is the wrong setting.** Grus and corestones are
  vadose-zone products; stagnant saturated water equilibrates.
- **No-flow side walls are not neutral.** They manufacture a domain-scale
  circulation and an 89 % block-to-block spread. The section is periodic in x.
- **A continuous rotation slider cannot preserve exact tiling** -- only lattice
  angles do, hence the snapping.
- **Two mechanisms for the horizontal-joint effect were measured and
  falsified** before the third (lateral bypass around the matrix, which is
  where the dissolving happens) was supported.
- **Defect correction does not beat BiCGSTAB here.** `x <- x + LU^-1(b - Ax)`
  cannot break down and needs one back-substitution per iteration where
  BiCGSTAB needs two, so it looked certain. Alternated in one process it is a
  wash: x1.06, x0.97, x1.17. Written, measured, and left on the stash.
- **The BiCGSTAB breakdowns are not staleness.** All 14 refactorisations over
  108 steps are `info = -10` after two iterations, which happens because the
  preconditioner is very good. The iteration cap plays no part: raising it to
  a million changes nothing.
- **A convergence study that stalls is probably not measuring convergence.**
  `run(years=X)` stepped PAST X by up to one step, so two settings were
  compared at two different model times and the gap was read as the coarser
  one's error. It floored a study at 1e-2 and made it look non-monotone
  (`f2645cd`).
- **An uncontrolled first step floors the whole run.** One opening step of
  7124 yr held the error at 1.1e-2 however tight the budget; controlled, the
  same run reached 6.3e-5. This is why the step control rejects rather than
  merely predicts.
- **Reusing the LU across a transport-operator rebuild is 2.1x SLOWER.** It is
  only a preconditioner, so keeping it looked free: 146 factorisations avoided
  out of 385, answer identical to 3e-10. It cost 10.51 s against 4.92, because
  the stale preconditioner costs more in extra BiCGSTAB iterations than the
  factorisation saves. Do not retry this.
- **Raising the step cap for Show trades the front for the speed.** Show does
  not animate, so stepping at the animation frame length looks wasteful:
  raising `dt_max` from 1 kyr to 50 is 2.3x faster and the bulk answer barely
  moves (0.9760 to 0.9746 dissolved). But `max|dM| = 0.22` -- the weathering
  front lands somewhere else, which is the thing a reader is looking at.
  Tightening the drift budget to compensate gives back the speed and not the
  accuracy.
- **Count solves, not steps.** The step control cut steps by 30 % and the run
  took exactly as long, because it was paying two solves per step in
  rejections.
- **`spilu` is not worth it as the preconditioner.** Half the build cost of a
  full `splu`, but the factorisation is reused across many steps, so build
  cost is not what dominates.

## (g) Reproduction

```sh
pip install -e ".[test]" && pytest                       # 60 tests, ~100 s
PYTHONPATH=src python3 examples/figure_three_panel.py    # the 20 m figure
PYTHONPATH=src python3 examples/figure_three_panel.py \
    --width 5 --depth 5 --rotation 45 --xperiod 2.0 --kyr 25 \
    --out examples/figure_45deg.png
PYTHONPATH=src panel serve examples/app.py --show        # the demo
artesian check numpy scipy matplotlib                    # browser feasibility
# probe_c needs fractopo, deliberately not a dependency:
PYTHONPATH=src <venv-with-fractopo>/bin/python prototypes/probe_c_topology.py
```

## The standing lesson

*(the original entry follows the additions of 2026-09-04)*

**A default tolerance can make a test unable to fail.** `np.allclose` carries
`atol=1e-8`. Every conductance and diffusivity in this model is smaller than
that, so several transcription tests -- the ones whose whole purpose is to
catch a docstring drifting from its code -- could not fail, and had not been
able to for their entire existence. Removing the default exposed two real
errors immediately. Every `np.allclose` in the suite now passes `atol=0.0`.

**Nothing from a page-summarising tool becomes a citation without reading the
source.** A fetch tool returned "quartz 1.16, albite 0.93, orthoclase 0.89
J/m2, Lawn & Marshall (1979)". The PDF it summarised contains none of those
minerals and no such reference. The numbers and the citation were both
invented by the summarising model, and they reached a design document before
a grep caught them.

**An unanchored substitution hits the neighbouring row.** Three times in one
day: a regex on a provenance SHA rewrote GRLP's and artesian's rows, a rename
of `m.k_weathered` missed `m.k_matrix *` and left a test comparing corrected
against uncorrected values, and a blanket rename of "transport-limited"
rewrote the sentence that defined the term. Anchor them, then diff.

**Send an agent to attack a result, not to confirm it.** Asked to confirm the
cracking criterion, an agent instead found that `f_FeO` was wrong by 4.5x and
that the fracture energy had been taken from one end of a 200-fold range --
two errors pulling opposite ways, whose product had looked like a prediction
landing on a measurement. Confirmation would have found supporting citations
for a wrong number.

Every real defect this model has had was caught by looking at a picture or by
reading an equation against its code -- never by an aggregate statistic. The
`L_eq` bug moved joint-versus-matrix dissolution by a factor of three thousand
while the mean dissolved fraction changed in the third decimal. Rebuild and
send the figure after every substantive change.

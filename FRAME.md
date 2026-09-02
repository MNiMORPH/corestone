# Frame -- read this first

The read-first frame for this repo, per `~/.claude/COMPACTION_PLAYBOOK.md`.
After a compaction, read this **before** acting, and verify every structural
claim below against git and disk before trusting it. Current to `0c76150`.

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

0. **Two decisions are waiting on Andy**, both measured and neither
   implemented -- see (e). The exponential integrator is the larger
   acceleration available to this model, and it is a parameter change, so it
   is his call and not a detail.
1. **Next: `artesian build` the demo.** `examples/app.py` exists and runs
   locally; it has never been compiled. That is the real test of whether a
   scipy sparse solve behaves under Pyodide.
   `HANDOFF-geomorphonline-demo.md` (untracked, written 2026-09-02 by the
   session that shipped the GRLP exercise) is the guide for this, and Andy
   has flagged it as the thing to take up next.
2. **Re-measure the stale design documents** (see the warning in (e) -- this is
   not optional bookkeeping, they contain numbers that are now wrong).
3. `f_inert` is set, claimed in the README and `design/02` as the second solid
   phase, and **never used**. Implement it or delete the claim.
4. `P21` carries a systematic +8 % bias from counting wall joints at full
   length; `test_a_full_orthogonal_grid_has_the_analytic_intensity` hides it
   behind `rel=0.15`.
5. Read `E_a` and `delta_H_r` out of Palandri & Kharaka rather than assuming.

## (c) Key current data and objects

Branch `master`, HEAD `0c76150`, **25 commits unpushed**, working tree clean
apart from the untracked `HANDOFF-geomorphonline-demo.md`.

- `src/corestone/fractures.py` -- the joint network. Seeded via `seed()`, or
  supplied wholesale via `from_masks()`. `tiling_angles()` / `tiling_spacing()`
  give the orientations and spacings that tile a periodic width.
- `src/corestone/weathering.py` -- **the model**. Steady Darcy flow solved once;
  steady advection-diffusion-reaction for the solute, one sparse solve per step
  with the LU cached and reused as a preconditioner.
- `examples/app.py` -- the four-slider demo. `panel serve examples/app.py`.
- `examples/figure_three_panel.py` -- the figure, parameterised:
  `--width --depth --dx --spacing --kyr --rotation --xperiod --out`.
- `prototypes/probe_[a-e]_*.py` -- the evidence behind each design decision.
- `tests/test_stated_equations.py` + `test_equation_coverage.py` -- the
  transcription tests and the ledger that stops an equation entering a
  docstring without a check. **These exist because a docstring drifted from its
  code for six revisions.**
- `tests/test_solver.py` -- the properties the fast paths rely on. Structural
  symmetry of both matrices (which is what licenses the ordering), the
  diagonal being present in the pattern (which is what licenses writing the
  reaction term in place), and the warm-started field matching a cold direct
  solve (which is what makes the reused guess a guess and not an
  approximation).

74 tests, ~14 s. **The "~100 s" this file used to claim was wrong** -- measured
at 14.4 s on the pre-speed-up tree, so it was never true, not merely stale.

## (d) Guardrails and irreversibility state

`github.com/MNiMORPH/corestone` is **public**; pushing publishes. 17 commits
are local only. Tags, releases, version bumps and closing issues each need
explicit authorisation in the current message.

Two pushes exist in the reflog: `0a57fcd` (mine, authorised) and `18cb016`
(not from my session -- Andy's).

## (e) Results, each with the method that verified it

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

**Two measured proposals, not implemented, awaiting a decision:**

- **The exponential integrator** (`prototypes/probe_f_integrator.py`,
  `0c76150`). `r` is proportional to `M`, so with `c` held over a step the
  equation is `dM/dt = -lambda M` and its solution is an exponential. Euler
  takes the tangent. At equal cost the exponential form is ~6x more accurate;
  at matched error over 200 kyr it is **18 steps against 107**, the largest
  single acceleration left. It needs `dt_max` and `dx_max` relaxed, which are
  parameters, and the error is not monotone in `dx_max`.
- **The refactorisation cap.** `max_krylov_iterations = 15` has never fired;
  lowering it to 5 so that it does is worth ~15 %. Left alone because the
  optimum moves between 4 and 5 case to case. Table in the `solve_solute`
  docstring.

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

Every real defect this model has had was caught by looking at a picture or by
reading an equation against its code -- never by an aggregate statistic. The
`L_eq` bug moved joint-versus-matrix dissolution by a factor of three thousand
while the mean dissolved fraction changed in the third decimal. Rebuild and
send the figure after every substantive change.

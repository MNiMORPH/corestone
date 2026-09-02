# Frame -- read this first

The read-first frame for this repo, per `~/.claude/COMPACTION_PLAYBOOK.md`.
Keep it current to HEAD. After a compaction, read this **before** acting, and
verify every structural claim below against git and disk before trusting it.

## (a) Origin -- why this model exists

To teach one idea: weathering is a race between how fast water delivers fresh
fluid and how fast rock dissolves into it. The payload is the misconception it
corrects -- **a corestone is not tougher rock**. Same granite, same minerals,
same temperature; it survives because the water never reached it, or reached it
already saturated.

Intended delivery is an interactive browser demo built with
[artesian](https://github.com/MNiMORPH/artesian), which constrains the model
hard: it must run under Pyodide, so numpy and scipy only, nothing compiled, and
fast enough to move a slider.

## (b) Plan and trajectory -- as the next action

1. **Next**: the artesian front end: sliders for temperature, rainfall, joint
   spacing and elapsed time, plus reseeding the network. Two panels -- the
   affinity field `(1 - c)` beside the rock -- so the equation is seen rather
   than narrated.
2. Read `E_A` out of Palandri & Kharaka (2004) rather than assuming 60 kJ/mol.

Open decision left with the author: the abutting-set `density`. 1.0 gives clean
well-bounded blocks and is the current default; low density gives realistic
Y-dominated outcrop topology with less well-bounded blocks. See design 03.

## (c) Key current data and objects

- `src/corestone/fractures.py` -- seeds the joint network.
- `src/corestone/weathering.py` -- **the model**. Verified to reproduce the
  prototype bit-for-bit (max|diff| = 0 on flux, concentration and dissolved
  fraction) when it was promoted out of `prototypes/`.
- `prototypes/probe_a_fracture_seeder.py` -- distance-to-fracture measurement.
- `prototypes/probe_b_weathering.py` -- the timing and temperature tables;
  now imports the package rather than carrying its own copy.
- `prototypes/probe_c_topology.py` -- fractopo validation. Needs an environment
  with fractopo; it is not a dependency.
- `examples/figure_three_panel.py` / `.png` -- the three-panel figure.
- `design/01`, `02`, `03` -- fracture seeding; teaching scope; throughgoing
  joints.

## (d) Guardrails and irreversibility state

Public repository: `github.com/MNiMORPH/corestone`. Pushing publishes.
Tags, releases, version bumps and closing issues each need their own explicit
authorisation.

## (e) Results, each with the method that verified it

- Corestones emerge from the affinity term alone. At 100 kyr the median
  distance-to-joint of corestone cells is 0.40 m against 0.00 m for grus cells
  (`probe_b`).
- Raising temperature 30 K multiplies the rate constant more than tenfold, yet
  grus changes by 0.8 percentage points while corestone survival *rises* by 15
  -- hotter water saturates sooner and the work concentrates at the joints
  (`probe_b`).
- The generator reproduces the analytic intensity of an orthogonal grid,
  P21 = 2/S = 1.33 against 1.28 measured (test).
- The network has zero I nodes and 2.00 connections per branch by fractopo's
  Sanderson & Nixon classification (`probe_c`) -- though `spans=1` guarantees
  I=0 by construction, so this confirms intent rather than discovering it.
- Results are grid-independent: grus 28.4 / 29.2 / 29.6 % at dx = 0.20 / 0.40 /
  0.50 m.

## (f) Negative results, and why not

- **No existing fracture-network generator was usable.** Eight examined and
  tabulated in design 03. fractopo analyses but does not generate; FracSim2D is
  Python 2 with Windows-only compiled extensions; ADFNE is MATLAB; HatchFrac is
  C++; dfnWorks is 3D. Pyodide rules out compiled code regardless.
- **Free segments from a power-law length distribution do not make a network.**
  The joints came out shorter than the gaps between them and nothing linked up.
- **Connectivity is the wrong measure for corestones.** Removing fractured links
  can never disconnect anything, because a joint is a conduit and rock is
  continuous across it. The measure is distance from the network.
- **Full saturation is the wrong setting.** Grus and corestones are vadose-zone
  products; stagnant saturated water equilibrates and stops weathering. The
  saturated zone sets *where* the front sits, not what does the weathering.

## (g) Reproduction

```sh
pip install -e ".[test]" && pytest              # 23 tests
PYTHONPATH=src python3 prototypes/probe_b_weathering.py
PYTHONPATH=src python3 examples/figure_three_panel.py
PYTHONPATH=src python3 examples/seed_a_joint_network.py
# probe_c needs fractopo, which is not a dependency:
PYTHONPATH=src <venv-with-fractopo>/bin/python prototypes/probe_c_topology.py
```

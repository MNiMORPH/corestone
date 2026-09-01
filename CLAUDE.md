# corestone

Project-scoped notes. Global standing instructions live in `~/.claude/CLAUDE.md`;
this file holds only what is specific to *this* model.

## Commands

```sh
pip install -e ".[test]"    # from a clone
pytest                      # the whole suite
pyright                     # type check (config is in pyproject.toml)
```

## Conventions

These are assumed silently by every routine, so an unstated one becomes a bug
that looks like a result. Change one only deliberately, and grep for what
depends on it.

- **Units**: SI throughout. Length m, time s, mass kg, temperature K.
  Aqueous concentration mol m^-3 *of solution*. Mineral content mol m^-3 *of
  bulk rock*, not of solid -- the two differ by porosity and the difference is
  a standard place to lose a factor.
- **Vertical coordinate is DEPTH, positive downward**, zero at the ground
  surface. A weathering profile is naturally depth-referenced, and every
  boundary condition here is stated relative to the surface. `z = -depth` is
  computed only for plotting. (This inverts the usual elevation-positive-up
  convention of surface-process work; it is deliberate, and it is the one
  convention most likely to bite.)
- **Array indexing is `[iz, ix]`**: `iz` increases downward (row 0 is the
  ground surface), `ix` increases to the right. Shape is `(nz, nx)`.
- **Cells hold the matrix; links hold the fractures.** Cell-centred arrays are
  `(nz, nx)`: porosity, mineral content, saturation, temperature, cohesion.
  Fracture state lives on the links *between* cell centres, in two arrays:
  `(nz - 1, nx)` for vertical links and `(nz, nx - 1)` for horizontal links.
  Anything indexed by link must say which of the two it is.
- **Flux sign**: positive downward on vertical links, positive rightward on
  horizontal links. Infiltration is therefore a positive flux at the top.
- **Saturation** is the fraction of pore volume occupied by water, in [0, 1],
  defined separately for matrix pores and fracture apertures.
- **Time**: seconds internally. Years appear only at the input and output
  edges, and the conversion is named where it happens.
- **Boundary conditions**: infiltration flux prescribed at the surface
  (row 0); the sides are no-flow; the base is the drainage boundary. Stated
  here rather than in the solver so that changing it is a visible act.

## Layout

- `src/corestone/corestone.py` -- the model: state, parameters, time evolution.
- `tests/` -- one test per claim; the test name states the claim.
- `examples/` -- runnable scripts, meant to be read.
- `design/` -- a design document per decision, written *before* the code.
- `prototypes/` -- throwaway executable probes that settle a design question.

## Working rules for this repo

- A design decision gets a `design/*.md` and a runnable probe in `prototypes/`
  **before** the implementation. If `prototypes/` is empty while `src/` is
  growing, that step is being skipped.
- Any threshold, cut-off, filter, or default not asked for is a proposal. Say
  it in one sentence; do not implement it silently.
- A regression test must be shown to fail without its fix.

## Provenance instrumentation (not yet wired)

The `trust/` kit (provenance banner + stop-gate hooks) is under evaluation in
`~/projects/lidar-diff-icp` and is deliberately *not* vendored here yet. If it
proves out, copy `trust/` in and add its hook stanza to
`.claude/settings.local.json`.

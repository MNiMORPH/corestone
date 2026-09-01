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

FILL THESE IN BEFORE WRITING PHYSICS. They are assumed silently by every
routine, so an unstated one becomes a bug that looks like a result.

- **Units**: SI throughout unless a docstring says otherwise.
- **Sign of z**: (positive up? positive down?)
- **Grid axes**: what the first and second index mean, and which way each runs.
- **Time**: seconds internally; converted only at the input and output edges.
- **Boundary conditions**: which are held, which are free, and where.

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

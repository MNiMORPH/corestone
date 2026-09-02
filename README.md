# corestone

**Fracture-controlled chemical weathering of granite, in a 2D vertical cross
section.**

Named for what the model leaves behind. A corestone is not tougher rock — it is
the same granite, the same minerals, the same temperature as the sand around
it. It survives because the water never reached it, or reached it already
saturated. The rock the water *did* reach falls apart into grus.

The physics is one equation:

```
R = k(T) · A · (1 − C/C_eq)
```

Dissolution runs at an Arrhenius rate constant multiplied by how far the pore
water is from equilibrium. Water that has equilibrated stops weathering rock,
however soluble the rock and however warm the water. Fresh water descends the
joints, so the joints decide where weathering happens.

The quantity that decides everything is the **equilibration length**

```
L_eq = q · C_eq / (k(T) · A)
```

— how far water travels before it is saturated. Rock further than `L_eq` from a
joint never sees undersaturated water. That is a corestone.

## Status

**Early and incomplete.** This is a *teaching* model, built to be run as an
interactive browser demo. What exists:

| | |
| --- | --- |
| `corestone.FractureNetwork` | implemented and tested — seeds a conjugate joint network |
| the weathering step | working, but still in `prototypes/`, not yet a module |
| the browser front end | not started |

**Every physical parameter in this model is a placeholder.** None is measured.
They are tabulated as such in `design/02-teaching-scope.md`, and the figure says
so on its face. Do not take a number out of this repository and use it.

## What it deliberately does not do

Simplifications made on purpose, and recorded with their costs in `design/`:

- **Flow is steady, gravity-driven descent**, one sweep from surface to base —
  no pressure solve, and no unsaturated (Richards) flow. Under partial
  saturation, wide joints can act as capillary *barriers* rather than conduits,
  which this model cannot represent. That is a real loss and it is deliberate.
- **Two solid phases**, one soluble and one inert, with no aqueous speciation,
  no secondary minerals, no oxygen, and no biotite oxidation.
- **No aperture evolution**, so no dissolution-driven channelization.
- **Disaggregation is a threshold**, not mechanics.

## Sources

The method, rather than the parameters, rests on:

- **Palandri & Kharaka (2004)**, *A compilation of rate parameters of
  water–mineral interaction kinetics*, USGS Open-File Report 2004-1068 — the
  rate-law form and the constants this model still needs.
- **Sanderson & Nixon (2015)**, the X/Y/I node classification used to check
  that the generated joint network is topologically plausible.
- **Fadakar-Alghalandis (2017)**, ADFNE, and **FracSim2D** (TU Delft) for the
  discrete fracture network construction. Neither was usable here — see
  `design/03-throughgoing-joints.md` — so the method is reimplemented and the
  network is validated against
  [fractopo](https://github.com/nialov/fractopo).
- **Rempe & Dietrich (2014)**, *PNAS*, and the vadose-zone weathering
  literature, for why the model sits above the water table rather than below
  it. Summarised in `design/02-teaching-scope.md`.

There is no paper to cite yet. `CITATION.cff` describes the software.

## Installation

Not on PyPI. From a clone:

```sh
git clone https://github.com/MNiMORPH/corestone.git
cd corestone
pip install -e ".[test]"
pytest
```

If your computer shields the core Python install from external packages, either
pass `--break-system-packages` (fine in my experience, but packages can clash)
or build a separate environment.

## Using it

```python
import numpy as np
from corestone import FractureNetwork, conjugate_sets

# A 20 x 15 m section at 20 cm resolution.
net = FractureNetwork(nz=75, nx=100, dx=0.2).seed(
    sets=conjugate_sets(dip_primary=90.0, dip_secondary=0.0, spacing=1.5),
    rng=np.random.default_rng(12345))

d = net.distance_to_fracture()       # metres from each cell to the nearest joint
print(net.p21, np.median(d))         # fracture intensity, median distance
```

`examples/seed_a_joint_network.py` runs that and plots it.

## Repository layout

- `src/corestone/` — the model.
- `design/` — a design document per decision, written *before* the code, with
  the probe that settled it and the parameters it introduced.
- `prototypes/` — the runnable probes themselves. Ugly on purpose.
- `tests/` — one test per claim; the test name states the claim.

## Conventions

Units, signs and grid orientation are stated in [CLAUDE.md](CLAUDE.md) and
assumed everywhere in the code — in particular that **the vertical coordinate
is depth, positive downward**, which inverts the usual convention of
surface-process work. Read that before changing anything numerical.

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).

# 01 -- Seeding the joint network

## The question

Granite weathers where water reaches it. Water reaches it along joints. So the
seeded fracture network decides, before any chemistry exists, which rock can
weather and which cannot. Does a plausible joint network leave rock at
distances a weathering front can actually cross -- and what network parameters
control that?

## Two corrections this probe forced

**A fracture is a conduit, not a barrier.** The first version of Probe A
measured connected components of cells with fractured links removed, i.e. it
treated joints as walls partitioning the rock into blocks. Rock is continuous
across a joint, mechanically and hydraulically, so that measurement can never
disconnect anything: a single fracture spanning the entire domain "leaked" at
every orientation tested. The test was wrong, not the seeder.

**A corestone is defined by distance, not by enclosure.** It is rock the water
never reaches -- rock whose distance from the network exceeds how far weathering
penetrates from a fracture wall in the time available. The right measurement is
therefore the distance transform from the network, and the right summary is its
distribution, not its maximum. The maximum is dominated by the domain edge.

## Options, with costs

**Placement.** Poisson-random positions are one line of code but cluster, and
clustered joints do not make blocks. Spacing-controlled placement -- offsets
drawn along each set normal -- costs a loop and reproduces the blocky structure
granite actually has. *Chosen: spacing-controlled.*

**Edges.** Seeding only inside the domain leaves edge cells artificially far
from any fracture, because no fracture can be centred just outside. Seeding
over a padded box and clipping costs one parameter (`PAD`) and removes the
artifact. *Chosen: padded.*

**Number of sets.** Two conjugate near-vertical sets are the textbook granite
case and were the plan. The probe says they are not enough (below).

## The probe

`prototypes/probe_a_fracture_seeder.py`. Three sets on 20 x 15 m at dx = 0.05 m
(300 x 400 cells), then the distance transform from the network -- benchmarked
against *the same measurement* on a perfectly regular, fully persistent network
at the same nominal spacing, so the comparison is like for like rather than
against an analytic guess.

```
fractures           56
P21 intensity       0.685 m/m2   (trace length per area)
distance to nearest fracture [m]:
  median 0.43   p90 1.12   max 2.31
  within 0.05 m of a fracture:  10.1 % of rock
  within 0.25 m of a fracture:  33.2 % of rock
  within 1.00 m of a fracture:  85.8 % of rock

REFERENCE -- same three sets, perfectly regular, fully persistent, spacing 1.5 m:
  median 0.15   p90 0.40   max 0.72
```

### Findings

**1. Persistence controls connectivity, and the first parameters were wrong.**
With trace lengths from a power law with median ~1 m against a 1.5 m spacing,
joints were *shorter than the gaps between them*: 29 fractures, P21 = 0.17,
median distance 1.12 m, max 5.98 m. Sweeping persistence L/S:

```
  Lmin  medL/S  P21    med d   p90 d
   0.5     0.7  0.17    1.12    2.55
   1.0     1.3  0.29    0.85    2.10
   2.0     2.7  0.42    0.63    1.57
   4.0     5.3  0.58    0.46    1.27
  16.0    21.3  1.05    0.29    0.85
```

Returns flatten past L/S ~ 5. Real joint sets are persistent, so this is
physically the right end of the range anyway.

**2. Two near-vertical sets cannot bound a block vertically.** No amount of
persistence closes it -- at L/S = 21 the median stalls at 0.29 m. Adding a
sub-horizontal (sheeting) set is what moves it (median 0.29 -> 0.20, p90 0.85
-> 0.58 at matched persistence). Sheeting joints were scheduled as optional and
later; the measurement says they are structural.

**3. A realistic network is about three times coarser than its nominal
spacing.** Median 0.43 m against the regular network's 0.15 m at the same
nominal 1.5 m. Lognormal spacing variability and finite trace length do that,
and it is real -- it is what makes the occasional large corestone. But it means
"spacing = 1.5 m" does not deliver 1.5 m blocks, and the nominal spacing must
be chosen against the measured distance distribution, not assumed.

## Decision

Spacing-controlled placement, padded seeding, three sets (two conjugate
near-vertical plus one sub-horizontal), von Mises orientation scatter, power-law
trace lengths at L/S ~ 5, lognormal spacing. Sets are a list, so a fourth is
configuration rather than code.

Reported per run: P21, and the distance-to-fracture distribution with the
regular-network reference beside it.

## Parameters chosen here

**Every value below is proposed, not measured.** They are placeholders until
field joint spacings and orientations replace them.

| parameter | value | note |
| --- | --- | --- |
| domain | 20 x 15 m | holds tens of blocks; from the solve-cost probe |
| `dx` | 0.05 m | resolves a cm-scale rind; 0.76 s per pressure solve |
| set dips | +75, -75, +5 deg | conjugate pair plus sheeting |
| `kappa` | 20, 20, 40 | von Mises orientation scatter |
| spacing | 1.5 m, lognormal sigma = 0.35 | **effective spacing is ~3x this** |
| trace length | power law, exponent 2, min 4 m | L/S ~ 5 |
| `PAD` | 3.0 m | seed outside, clip |

## Open

- Nominal spacing must be set from field data, against finding 3.
- Aperture distribution is not yet seeded; it belongs with the flow solver,
  since aperture only means something once it becomes conductance.
- Whether the sheeting set should have depth-dependent spacing -- sheeting
  joints tighten downward in real granite -- is deferred to the flow design.

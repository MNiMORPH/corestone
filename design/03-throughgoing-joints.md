# 03 -- Throughgoing joints, and an abutting conjugate set

## The question

Design 01 seeded both joint sets as free segments drawn from a power-law length
distribution. Real joint sets are not like that: they are *persistent*, and a
younger set *terminates against* an older one. Can an existing generator be
used instead, and does the change matter?

## Existing code: looked, and could not use any

| tool | why not |
| --- | --- |
| **fractopo** (Python, pip, JOSS) | analysis only -- validates and characterises trace maps, does not generate |
| **FracSim2D** (TU Delft, GPL-3) | Python 2, Windows-only `.pyd` extensions plus CGAL/Boost DLLs, needs the proprietary `digifract` core |
| **ADFNE** (Fadakar-Alghalandis 2017) | MATLAB |
| **HatchFrac** | C++, Mendeley data deposit rather than a maintained repository |
| **dfnWorks** (LANL) | 3D, requires LaGriT and PFLOTRAN |
| **GeoMop/dfn** | Gmsh dependency, undocumented, not installable |
| **pySimFrac** (LANL) | one rough fracture *surface*, not a network |
| **PorePy** | not on PyPI; takes fractures as input rather than generating them |

FracSim2D was the closest fit and was downloaded and unpacked before it was
ruled out. The `artesian` target settles it independently: the demo runs under
Pyodide, which supplies numpy and scipy and nothing compiled, so a C++ or
MATLAB dependency is not merely inconvenient but impossible.

**Adopted instead:** the published *method*, reimplemented, and validated
against an established tool. The construction is standard DFN practice; the
topological description is Sanderson & Nixon's X/Y/I node classification.

## The construction

- a **primary** set of *throughgoing* joints, each clipped from an infinite
  line to the domain, so both tips land on the boundary and never in intact
  rock;
- a **secondary** set that **abuts** the primary one: each trace is cut back to
  run from one primary joint to the next and terminate there.

`conjugate_sets(dip_primary=90, dip_secondary=0)` is the default -- vertical
joints cut by horizontal ones, 90 degrees apart, which is the simplest geometry
that bounds a block on all four sides. A symmetric shear pair is
`conjugate_sets(45, -45)`: still 90 degrees apart, differently oriented
relative to the surface. **Reading the request as the conjugate ANGLE being 90
degrees; the orthogonal default is one line to change if the intent was a
+/-45 pair.**

## What it bought

Same domain, same seed, same 20 x 15 m section at dx = 0.40 m:

```
                              traces  P21    median d   p90 d   max d
free segments, 3 sets (old)     56    0.69     0.40      1.08    2.24 m
throughgoing + abutting (new)   24    0.69     0.40      0.80    2.00 m
```

**Less than half the traces at the same intensity, and a tighter distribution.**
Persistent joints cover the rock more evenly than a scatter of short ones, which
is the point.

## Validation, by an outside tool

`prototypes/probe_c_topology.py`, measured with `fractopo`:

```
                          traces  P21     X    Y    I   Y/(Y+I)  CpB
3 sets, all throughgoing     41   1.86   303    0    0    0.00   2.00
conjugate 90/0, abutting     24   0.69     8   19    0    1.00   2.00
conjugate +/-45, abutting    33   0.70     7   15    0    1.00   2.00
conjugate 90/0, S = 0.8 m    43   1.21    27   35    0    1.00   2.00
```

Zero I nodes: no joint ends in intact rock. Connections per branch 2.00. By the
standard measure the network is fully connected.

**Two honest caveats.** First, `spans = 1` *guarantees* I = 0 by construction,
so fractopo is confirming that the implementation matches the intent rather
than discovering connectivity. Second, that makes the network **idealised**:
real outcrops do carry some I nodes, joints that die out in intact rock. Adding
a fraction of early-terminating traces would be more realistic; for teaching, a
cleanly connected network is the clearer object.

The all-throughgoing row is the contrast: with no abutting rule every joint
crosses every other, giving an X-dominated network at nearly three times the
intensity asked for. Outcrop networks are Y-dominated. The abutting rule is
what produces that.

## What the tests pin

Two tests failed on first writing, and both times the test was wrong rather
than the code -- worth recording, because the reason is a real property of the
geometry. With orientation scatter a "throughgoing" near-vertical joint seeded
near a side exits through *that side*, so its vertical extent is not the domain
depth. Throughgoing means **boundary to boundary**, not top to bottom. The
second failure came from classifying primary versus secondary traces by their
vertical extent, which fails for the same reason; traces now record which set
they came from (`segment_set`, `segments_of`), which is also what the figure
needs to colour them.

## Open

- Some fraction of I-node terminations, if realism outweighs clarity.
- Whether the sheeting set should have depth-dependent spacing; deferred with
  the rest of design 01's research scope.

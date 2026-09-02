# 05 -- The side walls were not neutral

## What was seen

Andy, looking at the figure: *"There is still less flow in the center than on
the edges."*

## What was measured

Two false starts worth recording, because both looked like answers.

**First**, quarter-integrating the downward flux gave edges 5.26e-8 against
centre 4.03e-8 -- a 30 % excess, apparently confirming it. It was a binning
artifact: the outer quarters happened to contain four joint columns and the
inner ones three.

**Second**, sampling the flux at joint columns derived from the trace
coordinates gave a factor of a thousand between the left and right halves.
That was floating-point drift in `round(x / dx)` landing off the joint, not
in the flux.

Taking the joint columns from the link mask instead:

```
flux down each joint / interior mean
0.987 0.990 0.994 0.999 1.003 1.006 1.007 1.007 1.006 1.003 0.999 0.994 0.990 0.987
```

**The flow field is uniform to about 1 %.** The observation was right but the
quantity was wrong: what varies is the weathering.

```
dissolved fraction per block, relative to the mean
1.03  1.20  1.23  1.17  1.00  0.70  0.34  0.70  1.00  1.17  1.23  1.20  1.03
                              ^ centre
```

A factor of 3.6 between the centre block and the blocks two in from the walls,
and perfectly symmetric.

## The cause

```
both sets,  19.55 m : 1.03 1.20 1.23 1.17 1.00 0.70 0.34 0.70 1.00 1.17 1.23 1.20 1.03
VERTICAL ONLY       : 1.30 1.08 0.97 0.92 0.90 0.89 0.89 0.89 0.90 0.92 0.97 1.08 1.30
both sets,  ~31 m   : 1.01 1.21 1.28 1.29 1.25 1.17 1.04 0.88 0.70 0.51 0.31 0.51 ...
```

The dip needs the subhorizontal joints, and it **scales with the width of the
section** rather than staying near the edges. So it is not a boundary layer: a
no-flow wall forces the lateral flow to vanish at both ends, and in a section
whose horizontal joints carry lateral flow that manufactures a domain-scale
circulation with a drainage divide down the middle. With vertical joints alone
there is only a mild two-block edge effect.

**A no-flow wall is a physical statement, not a neutral one.** It says the
section is bounded by impermeable rock. A section cut out of a jointed mass is
not.

## The fix

Periodic in x: the section wraps onto itself, so there are no walls.

- The Darcy assembly gains a wrap link joining column `nx - 1` to column `0`,
  conductive wherever a joint reaches both edges.
- The joint pattern must tile the period: `periodic_grid_shape()` returns
  `n * spacing / dx` columns, with no joint repeated at the seam, and **raises**
  if the spacing is not a whole number of cells. 1.5 m at 20 cm is 7.5 cells;
  silently rounding it would tile at 1.6 m while the caller believed 1.5 m.
- The abutting traces cross the seam, which appears as two pieces because the
  section is drawn cut open.
- The solute row balance becomes **cyclic** tridiagonal -- cell 0 draws on cell
  `nx - 1` and vice versa, putting entries in the far corners. Sherman-Morrison
  folds them back into two banded solves.

## Result

```
             block-to-block spread   mass-balance error
no-flow                     89.3 %             1.1e-09
periodic                     0.0 %             1.1e-09
```

Every block now weathers identically. Total weathering falls from 4.90 % grus
to 3.65 %, because the no-flow case was concentrating weathering into the
near-wall blocks rather than distributing it.

## Open

- Vertical remains non-periodic, correctly: the surface and the drainage base
  are real boundaries, not artifacts.
- The no-flow option is kept rather than deleted. It is the right model for a
  section genuinely bounded by impermeable rock, and it is a good illustration
  of a boundary condition quietly shaping an interior result.

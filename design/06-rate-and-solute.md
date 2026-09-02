# 06 -- The rate equation did not match the code

> **NUMBERS IN THIS DOCUMENT MAY BE STALE.** They were measured before one or
> more of: `a505892` (the saturation length made to scale with flux),
> `c0d7749` (diffusion, the C_eq temperature term, and the replacement of the
> whole transport operator), and `7cbd0a7` (non-axis-aligned joints made to
> conduct at all). The reasoning survives; re-measure before quoting a figure.
> See `FRAME.md` section (e).

## What was found

Andy: *"Search your code for the equation setting dissolution rate alongside
how you track solute concentrations in water."*

The module docstring states

```
L_eq = q * C_eq / (k(T) * A)
```

and the code computed

```python
L_eq = self.equilibration_length / np.maximum(self.M, 1e-6)   # no q anywhere
```

`L_eq` was a scalar per timestep. Since the per-cell dissolution is
`Q * beta * (1 - c)` with `beta = expm1(dx / L_eq)`, a uniform `L_eq` makes
**dissolution proportional to the local flux**: a joint carrying thirty times
the flux dissolved thirty times faster per unit volume at the same
undersaturation.

That is wrong. The rate per unit volume is `k A (1 - c)` -- a property of the
rock, not of how fast water moves past it. With `L_eq` proportional to `Q` the
`Q` cancels and the rate is flux-independent, which is the correct limit.

## The size of the error

Per-cell dissolution capacity `Q * beta`, relative to the mean infiltration:

```
                       joints    matrix     ratio
uniform L_eq (wrong)    1.066     0.003   348 : 1
L_eq ~ Q (correct)      0.100     0.871   0.115 : 1
```

A factor of three thousand in the relative weighting. The old code put
essentially all the dissolving in the joint cells; the correct physics puts the
*capacity* in the matrix, where it goes largely unused because the matrix water
is saturated.

**The totals barely moved** -- mean dissolved fraction 0.04598 to 0.04678 --
because what is dissolved is set by the undersaturation, not by the capacity.
So this was invisible in every aggregate number reported so far, and only shows
in where the dissolving happens.

## Open defect: solute is not conserved to better than a few percent

No solute enters at the surface, so everything produced by dissolution must
leave through the base. Measured over one step:

```
produced by dissolution   1.296e-07
leaving the base          1.265e-07
mismatch                       2.5 %
```

against a water balance exact to 1e-9. The residual of the per-cell solute
balance is concentrated entirely at the two boundary rows:

```
row 0            6.6e-02
rows 1-29        3.8e-04
interior         1.6e-04
last row         1.8e-01
```

The base boundary was changed from an overwritten Dirichlet row to a
conductance to an external fixed head, so that continuity holds in every cell
rather than being destroyed in the bottom row. That was right on its own terms
-- the water balance stayed exact and every cell now conserves -- but it did
not fix the solute residual.

Likely related: the concentration is not monotonic with depth, and the bottom
row runs lower than the rows above it, which is not explained.

**This is unresolved and is the next thing to chase.** Every aggregate result
in designs 04 and 05 carries a few percent of unexplained solute
non-conservation, which is small but is not nothing, and the boundary rows are
where the model is least trustworthy.

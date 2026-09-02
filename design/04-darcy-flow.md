# 04 -- Solve for the head, do not route the water

> **NUMBERS IN THIS DOCUMENT MAY BE STALE.** They were measured before one or
> more of: `a505892` (the saturation length made to scale with flux),
> `c0d7749` (diffusion, the C_eq temperature term, and the replacement of the
> whole transport operator), and `7cbd0a7` (non-axis-aligned joints made to
> conduct at all). The reasoning survives; re-measure before quoting a figure.
> See `FRAME.md` section (e).

## The question

Design 02 chose a gravity cascade for the flow: each cell hands its water to
the three cells below, split by their conductance. It was chosen to avoid a
pressure solve. Does that choice cost anything?

## It cost the entire horizontal joint set

A cascade moves water down and diagonally down. It can *enter* a subhorizontal
joint but never travel *along* one, so the horizontal set did nothing:

```
vertical + horizontal   133 traces   grus 24.65 %   corestone 61.89 %
vertical only            13 traces   grus 24.43 %   corestone 61.77 %
```

Deleting 120 of 133 traces moved the answer by 0.2 percentage points. On a
horizontal joint row the flux was 20 % above matrix flux, against a 1000:1
conductance ratio.

## The wrong fix, and why it was wrong

The first proposal was a rule: water arriving on a subhorizontal joint
*perches* and spreads along it before descending, mixing across each connected
run of fractured horizontal links. It has a real story behind it -- that is what
a capillary barrier does in the vadose zone -- and it was cheap and
parameter-free.

It was still a rule bolted on to preserve the cascade. Andy rejected it:
*"implement it as a general outcome of a water-flow equation and not as a
cheap-trick rule."* He was right, and the tell was already in the record: the
cascade had been justified by the cost of a pressure solve, and that cost had
been measured early -- 0.76 s at 120,000 cells, and **one-time**, because the
conductance never evolves. The justification had already been falsified.

## The flow equation

Steady Darcy flow for the hydraulic head:

```
div( K grad H ) = 0,        H = psi - d        (d is depth, positive down)
```

Infiltration prescribed at the surface, fixed head at the base, no-flow sides.
Conductance lives on the **links**, which is where the fracture network already
lives: a fractured link conducts, an intact one barely does. Lateral flow along
a joint is then not a special case -- it is what the head field does when a
low-resistance path exists.

One bug caught on the way: `H` is total head with elevation already in it, so
adding an elevation term to the link flux as well double-counts gravity. That
manufactured water -- 3500 % more leaving the base than entered the top.

## The transport that follows

A flow driven by the gradient of a potential cannot circulate, so the flux field
is acyclic -- measured, not assumed: **no vertical link carries upward flow.**
Rows can therefore still be swept in order. Within a row, water moving sideways
along a joint couples neighbouring cells, so each row is a tridiagonal system:

```
Q_i (1 + beta_i) c_i  -  sum(lateral inflow) c_j  =  Q_i beta_i  +  solute from above
```

with `beta = expm1(dx / L_eq)`. That choice is not an approximation: substituted
into the balance it reproduces the exact exponential
`c_out = 1 + (c_in - 1) exp(-dx/L_eq)` in one dimension, while keeping the
equation linear in two.

## What changed

```
                       grus %   corestone %   isolated corestones
cascade                 24.65        61.89              0
Darcy, no horizontal    10.31        83.02              0
Darcy, full network      9.50        78.91             33
```

**The corestones came back, in the perfect grid.** Removing the horizontal set
now changes the mean dissolved fraction by 1.1 %, and the horizontal joints
carry 51 % of all lateral flux. Water runs along them to reach the next vertical
joint, so weathering wraps around a block and closes it off from beneath --
which a purely vertical stripe never did.

**This corrects an earlier diagnosis.** With the cascade, only a network with
orientation scatter produced isolated corestones, and that was read as the
perfect verticality being at fault. It was not: it was the flow model. Testing
one deficiency while another is still present gives a confident and wrong
answer.

## It also corrected the temperature demonstration

Design 02 built the teaching demonstration on this: raise the temperature and
*more* rock survives as corestone, because hotter water saturates sooner and
concentrates the weathering at the joints. Under the cascade that was a clean
15-point rise, 52.1 to 67.5 %, across 275-305 K.

Under Darcy flow it nearly vanishes:

```
 T [K]  L_eq [m]     mean X   grus %  corestone %
   275     1.256    0.09510     8.70        84.42
   285     0.500    0.10220     9.08        82.47
   295     0.212    0.10387     9.26        84.09
   305     0.095    0.10424     9.70        86.06
   315     0.045    0.10427    10.21        87.49
```

A **28-fold** shortening of `L_eq` moves the total weathering by under a tenth,
and the last factor of four moves it by 0.3 %. The corestone fraction is not
even monotonic -- it dips between 275 and 285 K, and at coarse resolution the
275-to-305 trend reverses outright, which is how the test caught this.

The reason is that under a real flow equation the water is *already*
concentrated at the joints, so making `L_eq` smaller has little left to do. The
system is transport-limited: what dissolves is set by how much water arrives
and how much solute it can carry, not by how fast the rock dissolves.

That is a **sharper** lesson than the one it replaces, and a truer one -- the
high-Damkohler limit, where chemistry stops mattering entirely. But it is not
the lesson design 02 was written around, and the earlier dramatic version was
partly an artifact of the cascade. The slider that now does something is
rainfall, not temperature.

## Cost

```
head solve, 120,000 cells    ~3 s, ONCE (conductance does not evolve)
100 kyr of weathering        ~5 s   (was 3.7 s with the cascade)
figure, end to end           ~9.6 s
```

Mass balance is exact to 2e-16 at every row.

## Parameters this introduced

| parameter | value | note |
| --- | --- | --- |
| `k_fracture` | 1e-5 m/s | jointed granite. **Placeholder** |
| `k_matrix` | 1e-8 m/s | intact granite. **Placeholder** |
| base boundary | `psi = 0` | a drainage boundary at the base of the section |

The conductivity *ratio* is what sets the flow contrast; the absolute values
matter only against the infiltration rate, which they must exceed for the
matrix to accept water at all.

## Open

- The head-solve cost is now the floor on the browser demo. 10 cm resolution is
  equally converged and roughly 4x cheaper if the front end needs it.
- Unsaturated flow remains deferred: this is saturated Darcy, so joints cannot
  act as capillary barriers.

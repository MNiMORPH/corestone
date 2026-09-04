# 02 -- Teaching scope

> **NUMBERS IN THIS DOCUMENT MAY BE STALE.** They were measured before one or
> more of: `a505892` (the saturation length made to scale with flux),
> `c0d7749` (diffusion, the C_eq temperature term, and the replacement of the
> whole transport operator), and `7cbd0a7` (non-axis-aligned joints made to
> conduct at all). The reasoning survives; re-measure before quoting a figure.
> See `FRAME.md` section (e).

## The pivot

Design 01 and the literature review behind it were aimed at a research model:
variably saturated Richards flow with fracture-matrix capillary exchange,
multicomponent speciation, oxygen transport, cohesion mechanics. This model is
for **teaching, delivered as an `artesian` browser demo for visualisation and
equation learning**. Most of that machinery is not simplified away here -- it is
the wrong model for the purpose, and it is recorded in design 01 and below so
that it stays visible as a direction rather than quietly vanishing.

## The one idea

Weathering is a race between how fast water delivers fresh fluid and how fast
rock dissolves into it. The payload -- the misconception the model exists to
correct -- is that **a corestone is not tougher rock**. Same granite, same
minerals, same temperature. It survives because the water never reached it, or
reached it already saturated.

## The model, entire

Normalised concentration `c = C/C_eq` removes the need to assert a solubility.
Integrating the rate over a cell of height `dx` is then exact:

```
dc/dz = (1 - c) / L_eq        ->   c_out = 1 + (c_in - 1) * exp(-dx / L_eq)
```

with the **equilibration length**

```
L_eq = q * C_eq / (k(T) * A)          k(T) = k_ref * exp[-(Ea/R)(1/T - 1/T_ref)]
```

the distance water travels before it is saturated. This is the whole model in
one number, and it is the quantity to put on screen: corestones are the rock
further from a joint than `L_eq`.

Two solid phases. One dissolves; quartz does not. **The inert fraction is a
constant, not a state variable**, so the second phase costs one scalar: when the
soluble phase is gone, 30% of the original solid remains as loose grains. That
is the difference between grus and a cavity, and it is why the model can claim
grus at all -- a single soluble phase can only make a smooth dissolution front.

> **THIS WAS NEVER IMPLEMENTED, and the framing was wrong twice over.**
> `f_inert = 0.30` reached the parameter block and was read by nothing, so the
> second phase cost one scalar and bought nothing; it is now removed
> (`43895c9`). And quartz does not need declaring inert: `C_eq` is quartz
> saturation, so quartz sits at the ceiling and its driving force
> `(1 - C/C_eq)` is zero by construction. 0.30 was also the wrong number for
> the framing that replaced it -- it is the fraction for *everything but
> quartz*, both feldspars and biotite, whereas the kinetics now cited are
> plagioclase's alone.
>
> The paragraph's own argument still stands and is the reason to keep it
> visible: without grains left behind, this model cannot tell grus from a
> cavity. That is now recorded in the README as a limitation.

Flow was originally a gravity cascade -- each cell handing its water to the
three cells below -- chosen to avoid a pressure solve. **That was wrong, and it
is replaced by steady Darcy flow in design 04.** A cascade cannot move water
sideways, so the entire subhorizontal joint set was inert. See design 04 for
the replacement and what it cost.

## Probe B

`prototypes/probe_b_weathering.py`.

```
Weathering through time at T = 285 K:
   kyr   grus %  corestone %   mean X
    20      5.4         86.7    0.060
   100     29.3         60.8    0.295
   500     64.3         16.3    0.657
  1000     71.8         10.2    0.746
```

**The claim holds.** Corestones emerge from the affinity term alone: at 100 kyr
the median distance-to-joint of corestone cells is 0.40 m against 0.00 m for
grus cells. Nothing but `(1 - c)` and the plumbing produced that.

**The timescale is right, and it is geological.** A visible profile needs 10^5
to 10^6 years; the front advances at roughly 15 m/Myr, which is the observed
order for granite. A demo that runs 20 kyr shows a skin and teaches nothing.

### The demonstration (SUPERSEDED -- see design 04)

**The numbers below were produced by the gravity cascade and do not survive
the switch to Darcy flow.** Under a real flow equation a 28-fold change in
L_eq moves the total weathering by under a tenth, because the system is
transport-limited. Design 04 has the corrected sweep.


```
 T [K]  L_eq [m]   grus %  corestone %
   275     1.256     28.4         52.1
   285     0.500     29.3         60.8
   295     0.212     29.2         65.5
   305     0.095     29.2         67.5
```

A 30 K rise multiplies the rate constant by more than ten. **Grus changes by
0.8 percentage points, and corestone survival goes UP by 15.** Hotter water
saturates sooner, `L_eq` shrinks, and the work concentrates at the joints
instead of spreading. Warm does not mean weathered. That single slider move is
the lesson, and it inoculates against the commonest misconception in the field.

Then raising the rainfall slider *does* speed it up, and the mechanism clicks:
the limit was water supply, not chemistry.

## Numerics

Grid-independent to within the width of the result:

```
   dx   cells  grus %   mean X       dt cap    wall s
 0.20    7500    28.4    0.291         50 yr     1.49
 0.40    1900    29.2    0.295        200 yr     0.32
 0.50    1200    29.6    0.299        500 yr     0.28
```

**The 50-year step ceiling was mine, chosen arbitrarily, and it was the entire
performance bottleneck** -- the rate-based limiter never bound. Raising it to
500 years is 5.5x faster and changes the answer in the fourth significant
figure (0.2946 -> 0.2950). Settled at `dx = 0.40 m`, `DT_MAX_YR = 500`:
**0.28 s per 100 kyr**, interactive under Pyodide.

One diagnostic limitation, not a model limitation: the "fraction within L_eq of
a joint" column is meaningless once `L_eq < dx`, because the distance transform
cannot resolve it -- which is why 295 K and 305 K both report 30.5%. The model
itself is fine there; the per-cell integration handles `L_eq << dx` exactly
(`c_out -> 1`, the water saturates within the cell).

## Parameters

**None of these are measured.** They are placeholders producing the right orders
of magnitude, and the table exists so that no number in this model can be
mistaken for a result.

| parameter | value | note |
| --- | --- | --- |
| domain, `dx` | 20 x 15 m, 0.40 m | grid-independence checked |
| infiltration | 0.30 m/yr | plausible recharge |
| `L_EQ_REF` | 0.50 m at 285 K | **calibration choice**, not a measurement |
| `E_A` | 69.8 kJ/mol | oligoclase, neutral mechanism, Palandri & Kharaka (2004) Table 13 (was 60, unverified) |
| `TAU` = M0/C_eq | 6700 | saturated water volumes per rock volume |
| `X_GRUS`, `X_CORE` | 0.50, 0.05 | classification thresholds, chosen here |
| ~~`F_INERT`~~ | -- | removed in `43895c9`; never read, and quartz needs no declaring |
| K fracture : matrix | 1000 : 1 | routing conductance contrast |

## Deliberately deferred

Unsaturated flow and pore-water saturation (the largest remaining cut -- and
the largest scientific loss, since fracture-versus-matrix partitioning under
partial saturation is a real and unresolved control); aperture evolution and the
channelization feedback; multi-species transport, oxygen, and biotite
oxidation; cohesion mechanics as a solver rather than a threshold; permeability
feedback from porosity change.

## Open

- `E_A` should be read from Palandri & Kharaka rather than assumed.
- The `artesian` front end: two panels side by side -- the affinity field
  `(1 - c)` and the rock -- so the equation is seen rather than narrated.
- Students both move sliders **and** reseed the joint network (decided
  2026-09-02). Reseeding is itself a lesson: the same physics on a different
  network gives a different corestone field, which is the point -- the pattern
  is plumbing, not mineralogy.
- Sliders: temperature, rainfall, joint spacing, elapsed time. Temperature is
  the counterintuitive one and should be reachable first.

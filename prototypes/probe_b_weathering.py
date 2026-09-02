"""
Probe B: does the affinity term alone produce corestones and grus?

The whole model rests on one claim. Water enters the top fresh, descends
preferentially down the joints, and dissolves rock at

    R = k(T) * A * (1 - C/C_eq)

Where water moves fast it stays undersaturated and keeps working. Where it
creeps through the matrix it reaches C_eq, the bracket goes to zero, and
weathering stops -- with the rock just as soluble and the water just as warm.
The claim is that this, and nothing else, leaves corestones.

Working in normalised concentration c = C/C_eq removes the need to assert a
solubility. Integrating the rate over a cell of height dx is then exact:

    dc/dz = (1 - c) / L_eq        ->   c_out = 1 + (c_in - 1) * exp(-dx/L_eq)

with the EQUILIBRATION LENGTH

    L_eq = q * C_eq / (k * A)

the distance water travels before it is saturated. That is the whole model in
one number, and it is the quantity to put on screen: corestones should be the
rock further from a joint than L_eq.

Run:  PYTHONPATH=src python3 prototypes/probe_b_weathering.py
"""
import numpy as np

from corestone import (FractureNetwork, Weathering, orthogonal_grid,
                       periodic_grid_shape)

DX, SPACING = 0.05, 1.5
NZ, NX = periodic_grid_shape(20.0, 15.0, DX, SPACING)


def model(T=285.0):
    net = FractureNetwork(NZ, NX, DX, periodic_x=True).seed(sets=orthogonal_grid(SPACING),
                                           rng=np.random.default_rng(12345))
    m = Weathering(net)
    m.set_temperature(T)
    return m


net = FractureNetwork(NZ, NX, DX, periodic_x=True).seed(sets=orthogonal_grid(SPACING),
                                       rng=np.random.default_rng(12345))
dist = net.distance_to_fracture()
print("domain %d x %d at dx = %s m (%.0f x %.0f m), %d traces"
      % (NZ, NX, DX, net.lz, net.lx, len(net.segments)))
print("distance to joint: median %.2f m, p90 %.2f m"
      % (np.median(dist), np.percentile(dist, 90)))
print()
print("Weathering through time at T = 285 K:")
print("%6s %8s %12s %8s" % ("kyr", "grus %", "corestone %", "mean X"))
for kyr in (20, 100, 500, 1000):
    m = model().run(years=kyr * 1e3)
    print("%6d %8.1f %12.1f %8.3f" % (kyr, m.is_grus.mean() * 100,
                                      m.is_corestone.mean() * 100,
                                      m.dissolved_fraction.mean()))

print()
print("Turn up the temperature and the weathering does NOT keep pace: hotter")
print("water saturates sooner, so L_eq shrinks and the work concentrates at")
print("the joints. More rock survives as corestone, not less.")
print("%6s %9s %8s %12s" % ("T [K]", "L_eq [m]", "grus %", "corestone %"))
for T in (275.0, 285.0, 295.0, 305.0):
    m = model(T).run(years=100e3)
    print("%6.0f %9.3f %8.1f %12.1f" % (T, m.equilibration_length,
                                        m.is_grus.mean() * 100,
                                        m.is_corestone.mean() * 100))

"""
Probe F: is forward Euler the right integrator for the rock?

    d(M/M0)/dt = - r (1 - c) / tau,      r = r_ref * M * k(T)/C_eq(T)

r is PROPORTIONAL TO M -- the reactive surface area falls as the mineral is
consumed -- so with c held over the step the equation is dM/dt = -lambda M,
whose solution is an exponential and not a straight line. Forward Euler takes
the tangent, which always undershoots M and is why the step needs a dx_max
limiter and a clip at zero at all.

Substituting exp(-lambda dt) for (1 - lambda dt) is the same equation, more
exactly integrated, and costs one np.exp per step. The question is only how
much longer a step it buys at equal accuracy, since a step is one sparse solve.
"""
import time, numpy as np
from corestone import FractureNetwork, Weathering, orthogonal_grid, YEAR

def build(dx=0.05, nx=60, nz=61, spacing=0.5):
    net = FractureNetwork(nz, nx, dx, periodic_x=True).seed(
        sets=orthogonal_grid(spacing), rng=np.random.default_rng(12345))
    m = Weathering(net); m.set_infiltration(0.30/YEAR)
    return m

def run(dt_years, kyr, exponential, dx_max=0.05):
    m = build(); m.initialize()
    m.dt_max = dt_years * YEAR; m.dx_max = dx_max
    n = 0
    t0 = time.perf_counter()
    target = kyr * 1e3 * YEAR
    while m.t < target:
        r = m.reaction_coefficient
        c = m.solve_solute(r)
        if exponential:
            # lambda = r/M * (1-c)/tau, finite as M -> 0 because r is prop. to M
            lam = (m.reaction_coefficient / np.maximum(m.M, 1e-300)) \
                  * (1.0 - c) / m.tau
            lam = np.where(m.M > 0, lam, 0.0)
            step = min(m.dt_max, m.dx_max / max((lam * m.M).max(), 1e-30))
            m.M = np.clip(m.M * np.exp(-lam * step), 0.0, 1.0)
        else:
            rate = r * (1.0 - c) / m.tau
            step = min(m.dt_max, m.dx_max / max(rate.max(), 1e-30))
            m.M = np.clip(m.M - rate * step, 0.0, 1.0)
        m.c = c; m.t += step; n += 1
    return time.perf_counter()-t0, n, m.M

for KYR in (50, 200):
    print("=== %d kyr.  reference: Euler, dt_max = 50 yr, dx_max = 0.002" % KYR)
    _, nref, REF = run(50, KYR, False, dx_max=0.002)
    print("    %d reference steps" % nref)
    print("    %-12s %-9s %-8s %7s %7s %11s"
          % ("integrator", "dt_max", "dx_max", "time", "steps", "max|dM|"))
    rows = [("euler", 2000, 0.05)]          # as committed
    for dt in (5000, 20000):
        for dxm in (0.05, 0.10, 0.20, 0.50, 1.00):
            rows.append(("exponential", dt, dxm))
    for kind, dt, dxm in rows:
        t, n, M = run(dt, KYR, kind == "exponential", dx_max=dxm)
        print("    %-12s %-9d %-8.2f %6.2fs %7d %11.2e"
              % (kind, dt, dxm, t, n, np.abs(M-REF).max()))
    print()

# ---------------------------------------------------------------- result
#
# 200 kyr, 3 x 3 m at dx = 0.05, against a 4,071-step reference:
#
#     integrator   dt_max   dx_max   steps   max|dM|
#     euler          2000     0.05     107   1.30e-04     <- as committed
#     exponential    5000     0.05      64   6.28e-05
#     exponential   20000     0.10      30   3.80e-04
#     exponential   20000     0.20      18   1.51e-04     <- same error, 6x
#     exponential   20000     0.50      12   8.59e-04
#
# Two separate findings, and they should not be run together:
#
# 1. AT EQUAL COST the exponential form is about six times more accurate
#    (5.87e-3 against 3.52e-2 at 50 kyr, both at 28 steps). That is free: the
#    same solves, one np.exp added. It also cannot drive M below zero, so the
#    clip becomes a guard rather than a mechanism.
#
# 2. TURNING THAT INTO SPEED means relaxing dt_max and dx_max, and those are
#    parameters, not consequences. 18 steps against 107 at matched error is
#    the largest single acceleration available to this model -- larger than
#    everything in the solver -- but it is Andy's call, not a detail, because
#    it changes what the limiters mean and the error is not monotone in
#    dx_max (0.10 is worse than 0.20 at 50 kyr, since the limiter decides
#    which steps land where).
#
# Nothing here is implemented. This file is the measurement behind the
# proposal.

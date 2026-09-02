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

from corestone import FractureNetwork

YR = 365.25 * 24 * 3600.0

# ---- domain, teaching resolution ----------------------------------------------
NX, NZ, DX = 50, 38, 0.40             # 20 m wide x 15 m deep
DT_MAX_YR = 500.0                     # step ceiling; see design/02

# ---- PLACEHOLDER PARAMETERS -- none of these are measured ---------------------
INFILTRATION = 0.30 / YR              # m/s  (0.3 m/yr recharge)
K_FRACTURE = 1000.0                   # routing conductance, fracture cell
K_MATRIX = 1.0                        # routing conductance, intact rock
L_EQ_REF = 0.50                       # m, equilibration length at T_REF
T_REF = 285.0                         # K
E_A = 60.0e3                          # J/mol, feldspar-ish -- UNVERIFIED,
R_GAS = 8.314                         #        needs Palandri & Kharaka (2004)
TAU = 6700.0                          # M0/C_eq: saturated water volumes per
                                      #          rock volume to dissolve it
X_GRUS = 0.50                         # soluble fraction lost -> disaggregates
X_CORE = 0.05                         # below this, effectively unaltered
F_INERT = 0.30                        # quartz: never dissolves, becomes sand


def equilibration_length(T):
    """L_eq shrinks as temperature rises: hotter water saturates sooner."""
    return L_EQ_REF * np.exp((E_A / R_GAS) * (1.0 / T - 1.0 / T_REF))


def route_flow(fn):
    """
    Steady gravity-driven descent. One sweep, top to bottom: each cell hands
    its water to the three cells below, split by their conductance. No pressure
    solve -- gravity makes the grid a DAG ordered by depth.
    """
    K = np.where(fn.cell, K_FRACTURE, K_MATRIX)
    q = np.zeros((NZ, NX))
    q[0, :] = INFILTRATION * DX                      # m2/s per unit thickness
    for iz in range(NZ - 1):
        below = K[iz + 1, :]
        wl = np.concatenate([[0.0], below[:-1]])     # down-left receiver
        wr = np.concatenate([below[1:], [0.0]])      # down-right receiver
        tot = wl + below + wr
        f_l, f_c, f_r = wl / tot, below / tot, wr / tot
        send = q[iz, :]
        q[iz + 1, :] += send * f_c
        q[iz + 1, :-1] += (send * f_l)[1:]
        q[iz + 1, 1:] += (send * f_r)[:-1]
    return q, K


def weather(fn, years, T=T_REF, report=True):
    """Advance the rock state; return dissolved fraction X of the soluble phase."""
    q, K = route_flow(fn)
    L_eq0 = equilibration_length(T)
    M = np.ones((NZ, NX))                            # soluble mineral, M/M0
    t, dt = 0.0, DT_MAX_YR * YR

    while t < years * YR:
        # Surface area falls with the mineral that is left, so L_eq grows.
        L_eq = L_eq0 / np.maximum(M, 1e-6)
        c_in = np.zeros((NZ, NX))
        dissolved = np.zeros((NZ, NX))
        carry_q = np.zeros(NX)
        carry_qc = np.zeros(NX)
        for iz in range(NZ):
            qi = q[iz, :]
            ci = np.where(qi > 0, carry_qc / np.maximum(carry_q, 1e-300), 0.0) \
                 if iz > 0 else np.zeros(NX)
            c_out = 1.0 + (ci - 1.0) * np.exp(-DX / L_eq[iz, :])
            dissolved[iz, :] = qi * (c_out - ci) / DX     # mol-equivalents/m3/s
            c_in[iz, :] = ci
            # hand (q, q*c) down with the same splitting as the flow
            below = K[iz + 1, :] if iz < NZ - 1 else None
            if below is None:
                break
            wl = np.concatenate([[0.0], below[:-1]])
            wr = np.concatenate([below[1:], [0.0]])
            tot = wl + below + wr
            carry_q = np.zeros(NX); carry_qc = np.zeros(NX)
            for w, sl in ((wl / tot, -1), (below / tot, 0), (wr / tot, +1)):
                add_q, add_qc = qi * w, qi * w * c_out
                if sl == 0:
                    carry_q += add_q; carry_qc += add_qc
                elif sl == -1:
                    carry_q[:-1] += add_q[1:]; carry_qc[:-1] += add_qc[1:]
                else:
                    carry_q[1:] += add_q[:-1]; carry_qc[1:] += add_qc[:-1]

        rate = dissolved / (TAU * DX)                  # d(M/M0)/dt  [1/s]
        step = min(dt, 0.02 / max(rate.max(), 1e-30))
        M = np.clip(M - rate * step, 0.0, 1.0)
        t += step

    return 1.0 - M, q, c_in


fn = FractureNetwork(NZ, NX, DX).seed(rng=np.random.default_rng(12345))
dist = fn.distance_to_fracture()

print(f"domain {NZ} x {NX} at dx = {DX} m  ({NZ*DX} x {NX*DX} m), "
      f"{len(fn.segments)} fractures")
print(f"distance to joint: median {np.median(dist):.2f} m, "
      f"p90 {np.percentile(dist,90):.2f} m")
print()
print("Weathering through time at T = 285 K:")
print(f"{'kyr':>6} {'grus %':>8} {'corestone %':>12} {'mean X':>8}")
for kyr in (20, 100, 500, 1000):
    X, _, _ = weather(fn, years=kyr * 1e3)
    print(f"{kyr:6d} {(X > X_GRUS).mean()*100:8.1f} "
          f"{(X < X_CORE).mean()*100:12.1f} {X.mean():8.3f}")

print()
print("Turn up the temperature and the weathering does NOT keep pace: hotter")
print("water saturates sooner, so L_eq shrinks and the work concentrates at")
print("the joints. More rock survives as corestone, not less.")
print(f"{'T [K]':>6} {'L_eq [m]':>9} {'grus %':>8} {'corestone %':>12} "
      f"{'within L_eq %':>14}")
for T in (275.0, 285.0, 295.0, 305.0):
    X, q, _ = weather(fn, years=100000.0, T=T)
    L = equilibration_length(T)
    grus = (X > X_GRUS).mean() * 100
    core = (X < X_CORE).mean() * 100
    near = (dist <= L).mean() * 100
    print(f"{T:6.0f} {L:9.3f} {grus:8.1f} {core:12.1f} {near:14.1f}")

print()
X, q, c = weather(fn, years=100000.0, T=T_REF)
print(f"At T = {T_REF:.0f} K after 100 kyr:")
print(f"  mean dissolved fraction of the soluble phase: {X.mean():.3f}")
print(f"  grus (X > {X_GRUS}) is {(X > X_GRUS).mean()*100:.1f} % of the domain; "
      f"it retains {F_INERT:.0%} of its solid as loose quartz")
print(f"  corestones (X < {X_CORE}) are {(X < X_CORE).mean()*100:.1f} %")
print(f"  median distance-to-joint of corestone cells: "
      f"{np.median(dist[X < X_CORE]):.2f} m")
print(f"  median distance-to-joint of grus cells:      "
      f"{np.median(dist[X > X_GRUS]):.2f} m")

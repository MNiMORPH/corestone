"""
Probe D: does a real flow equation make the horizontal joints conduct?

The gravity cascade -- each cell handing water to the three cells below --
cannot move water sideways, so a subhorizontal joint could be entered but never
travelled along. Deleting the entire horizontal set changed the weathering by
0.2 percentage points. The horizontal joints were decorative.

The fix has to come out of a flow equation rather than a redistribution rule.
Solve steady Darcy flow for the hydraulic head,

    div( K grad H ) = 0,        H = psi - d      (d is depth, positive down)

with infiltration prescribed at the surface, a fixed head at the base, and
no-flow sides. Conductance lives on the LINKS between cells, which is where the
fracture network already lives: a fractured link is conductive, an intact one
is not. Lateral flow along a joint is then not a special case -- it is what the
head field does when a low-resistance path exists.

Two things this probe has to establish before any of it is worth building on:

  1. lateral flux appears on the horizontal joints, and is large;
  2. the flux field is acyclic, so solute can still be swept in one pass
     instead of solved implicitly every step. A flow driven by the gradient of
     a potential cannot circulate, so this should hold -- but it is the
     assumption the whole transport step rests on, so measure it.

Run:  PYTHONPATH=src python3 prototypes/probe_d_darcy.py
"""
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

from corestone import FractureNetwork, JointSet, orthogonal_grid

YR = 365.25 * 24 * 3600.0

NZ, NX, DX = 300, 401, 0.05
INFILTRATION = 0.30 / YR              # m/s
K_MATRIX = 1.0e-8                     # intact granite [m/s]   PLACEHOLDER
K_FRACTURE = 1.0e-5                   # jointed rock  [m/s]    PLACEHOLDER


def link_conductance(net):
    """K on every link: fractured links conduct, intact ones barely do."""
    kv = np.where(net.link_v, K_FRACTURE, K_MATRIX)      # (nz-1, nx)
    kh = np.where(net.link_h, K_FRACTURE, K_MATRIX)      # (nz, nx-1)
    return kv, kh


def solve_head(net, kv, kh):
    """
    Steady Darcy head on the cell centres.

    Finite volume on square cells, so the geometric factor is 1: the flux along
    a link is ``K * (H_i - H_j)`` in m2/s per unit thickness. Gravity enters
    through the base boundary condition and the elevation part of H.
    """
    nz, nx, dx = net.nz, net.nx, net.dx
    n = nz * nx
    idx = np.arange(n).reshape(nz, nx)

    rows, cols, vals = [], [], []

    def add(a, b, k):
        rows.extend([a, a]); cols.extend([a, b]); vals.extend([k, -k])

    a_v, b_v = idx[:-1, :].ravel(), idx[1:, :].ravel()
    a_h, b_h = idx[:, :-1].ravel(), idx[:, 1:].ravel()
    for a, b, k in ((a_v, b_v, kv.ravel()), (a_h, b_h, kh.ravel())):
        rows.append(a); cols.append(a); vals.append(k)
        rows.append(a); cols.append(b); vals.append(-k)
        rows.append(b); cols.append(b); vals.append(k)
        rows.append(b); cols.append(a); vals.append(-k)
    rows = np.concatenate(rows); cols = np.concatenate(cols)
    vals = np.concatenate(vals)
    A = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tolil()

    # Source: infiltration into the top row [m2/s per unit thickness].
    b = np.zeros(n)
    b[idx[0, :]] = INFILTRATION * dx

    # Base: fixed head, psi = 0 at the drainage boundary, so H = -depth.
    base = idx[-1, :]
    for i in base:
        A.rows[i] = [i]
        A.data[i] = [1.0]
    b[base] = -(nz - 0.5) * dx

    return spl.spsolve(A.tocsc(), b).reshape(nz, nx)


def link_fluxes(H, kv, kh, dx):
    """
    Flux along each link, from the gradient of the TOTAL head.

    ``H`` is total head, elevation included, so gravity is already in it --
    adding an elevation term here as well double-counts it and manufactures
    water (measured: 3500 % more leaving the base than entered the top).
    Positive is downward on vertical links and rightward on horizontal ones.
    """
    q_v = kv * (H[:-1, :] - H[1:, :])
    q_h = kh * (H[:, :-1] - H[:, 1:])
    return q_v, q_h


net = FractureNetwork(NZ, NX, DX).seed(sets=orthogonal_grid(1.5),
                                       rng=np.random.default_rng(12345))
kv, kh = link_conductance(net)

t0 = time.perf_counter()
H = solve_head(net, kv, kh)
t_solve = time.perf_counter() - t0
q_v, q_h = link_fluxes(H, kv, kh, DX)

print("head solve: %.2f s for %d cells" % (t_solve, NZ * NX))
print()
print("LATERAL flux on horizontal links [m2/s]:")
print("  on fractured links : mean |q| %.3e   max |q| %.3e"
      % (np.abs(q_h[net.link_h]).mean(), np.abs(q_h[net.link_h]).max()))
print("  on intact links    : mean |q| %.3e   max |q| %.3e"
      % (np.abs(q_h[~net.link_h]).mean(), np.abs(q_h[~net.link_h]).max()))
print("  ratio of the means : %.1fx"
      % (np.abs(q_h[net.link_h]).mean() / np.abs(q_h[~net.link_h]).mean()))
print()
print("VERTICAL flux on vertical links [m2/s]:")
print("  on fractured links : mean %.3e" % q_v[net.link_v].mean())
print("  on intact links    : mean %.3e" % q_v[~net.link_v].mean())
print("  ratio              : %.1fx"
      % (q_v[net.link_v].mean() / q_v[~net.link_v].mean()))
print()
lateral = np.abs(q_h).sum()
vertical = np.abs(q_v).sum()
print("total |lateral| / total |vertical| = %.3f" % (lateral / vertical))
print("fraction of vertical links flowing UPWARD: %.4f %%"
      % ((q_v < 0).mean() * 100))

# Mass balance: what goes in must come out of the base.
q_in = INFILTRATION * DX * NX
q_out = q_v[-1, :].sum()
print()
print("mass balance: in %.4e   out of base %.4e   (%.3f %%)"
      % (q_in, q_out, 100 * (q_out - q_in) / q_in))

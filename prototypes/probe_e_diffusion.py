"""
Probe E: does adding diffusion let blocks weather inward?

The model has pure advection. Water that is not on a flow path never carries
solute away, so a block interior saturates and then sits at c = 1 for ever: the
picture is binary, joints entirely dissolved and blocks entirely untouched,
with nothing in between. Real rock does not do that. The concentration gradient
between a saturated block interior and a fracture that is being flushed drives
**diffusive export of solute**, which keeps the interior slightly
undersaturated and lets it dissolve inward. That is how a weathering rind forms
at all, and it is also the geometric route to spheroidal rounding: a corner is
exposed on three faces, an edge on two, a face on one, so corners lose solute
fastest and round off.

Two transport coefficients, not one:

  D = D_molecular / tortuosity  +  alpha_L * |v|

The first is molecular diffusion. The second is **hydrodynamic (mechanical)
dispersion** from velocity variation within the pore space. It is often loosely
called turbulent diffusion, and the operator is identical, but at these fluxes
the flow is nowhere near turbulent -- this probe prints the Reynolds number to
make that concrete. The name matters only because it decides what sets
alpha_L: pore and aperture geometry, not eddies.

Adding diffusion breaks the one-pass sweep. Advection alone is acyclic, so the
solute could be swept row by row; diffusion is not directional, so the solute
field becomes a single coupled steady advection-diffusion-reaction solve over
the whole domain. This probe measures both whether it produces rinds and what
it costs.

A cleaner statement of the reaction while we are here. With L_eq proportional
to the flux, the reaction coefficient

    r = q / L_eq = q_ref / (dx * L_ref)          [1/s]

is a CONSTANT -- a property of the rock, with the flux cancelled out
explicitly rather than implicitly. Production per cell is r * dx^2 * (1 - c).

Run:  PYTHONPATH=src python3 prototypes/probe_e_diffusion.py
"""
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

from corestone import (FractureNetwork, Weathering, orthogonal_grid,
                       periodic_grid_shape)

YR = 365.25 * 24 * 3600.0

DX, SPACING = 0.05, 1.5
NZ, NX = periodic_grid_shape(20.0, 15.0, DX, SPACING)

# ---- transport coefficients. PLACEHOLDERS ------------------------------------
D_MOLECULAR = 1.0e-9        # m2/s, a typical aqueous ion
TORTUOSITY = 10.0           # matrix path tortuosity, dimensionless
ALPHA_L = 0.05              # longitudinal dispersivity [m]; of order the cell


def transport_coefficients(m, net):
    """D on every link: molecular plus velocity-dependent dispersion."""
    dx = m.dx
    # Darcy flux magnitude on each link, m/s.
    v_v = np.abs(m.q_v) / dx
    v_h = np.abs(m.q_h) / dx
    # Molecular diffusion is reduced by tortuosity in the matrix, not in a joint.
    dm_v = np.where(net.link_v, D_MOLECULAR, D_MOLECULAR / TORTUOSITY)
    dm_h = np.where(net.link_h, D_MOLECULAR, D_MOLECULAR / TORTUOSITY)
    # As a conductance for solute, in the same units as the fluxes Q (m2/s),
    # a link of length dx and face dx contributes D * (dx / dx) = D.
    return dm_v + ALPHA_L * v_v, dm_h + ALPHA_L * v_h


def solve_solute(m, net, D_v, D_h, r):
    """
    Steady advection-diffusion-reaction for the normalised concentration.

        sum_out f c_i - sum_in f c_j + sum_links D (c_i - c_j) + r dx^2 c_i
            = r dx^2

    Linear in c, so one sparse solve. Water entering the top carries c = 0, so
    it contributes nothing to the inflow sum -- the undersaturation enters the
    problem entirely through the reaction term on the right.
    """
    nz, nx, dx = m.nz, m.nx, m.dx
    n = nz * nx
    idx = np.arange(n).reshape(nz, nx)
    rows, cols, vals = [], [], []

    def add(i, j, v):
        rows.append(np.asarray(i).ravel())
        cols.append(np.asarray(j).ravel())
        vals.append(np.asarray(v).ravel())

    # --- advection, upwind, on every link
    for a, b, f in ((idx[:-1, :], idx[1:, :], m.q_v),
                    (idx[:, :-1], idx[:, 1:], m.q_h)):
        fwd = np.maximum(f, 0.0)          # a -> b
        rev = np.maximum(-f, 0.0)         # b -> a
        add(a, a, fwd);  add(b, a, -fwd)
        add(b, b, rev);  add(a, b, -rev)
    if net.periodic_x:
        fwd = np.maximum(m.q_wrap, 0.0)
        rev = np.maximum(-m.q_wrap, 0.0)
        add(idx[:, -1], idx[:, -1], fwd);  add(idx[:, 0], idx[:, -1], -fwd)
        add(idx[:, 0], idx[:, 0], rev);    add(idx[:, -1], idx[:, 0], -rev)
    # outflow through the base
    add(idx[-1, :], idx[-1, :], m.q_out_base)

    # --- diffusion, symmetric, on every link
    for a, b, D in ((idx[:-1, :], idx[1:, :], D_v),
                    (idx[:, :-1], idx[:, 1:], D_h)):
        add(a, a, D);  add(a, b, -D)
        add(b, b, D);  add(b, a, -D)
    if net.periodic_x:
        Dw = D_h[:, 0] * 0.0 + np.where(net.link_wrap, D_MOLECULAR,
                                        D_MOLECULAR / TORTUOSITY)
        add(idx[:, -1], idx[:, -1], Dw);  add(idx[:, -1], idx[:, 0], -Dw)
        add(idx[:, 0], idx[:, 0], Dw);    add(idx[:, 0], idx[:, -1], -Dw)

    # --- reaction
    src = r * dx * dx
    add(idx, idx, np.full((nz, nx), src))

    A = sp.coo_matrix((np.concatenate(vals),
                       (np.concatenate(rows), np.concatenate(cols))),
                      shape=(n, n)).tocsc()
    b = np.full(n, src)
    return spl.spsolve(A, b).reshape(nz, nx)


net = FractureNetwork(NZ, NX, DX, periodic_x=True).seed(
    sets=orthogonal_grid(SPACING), rng=np.random.default_rng(12345))
m = Weathering(net).initialize()

q_ref = m.infiltration * m.dx
r = q_ref / (m.dx * m.L_eq_ref)          # constant reaction coefficient [1/s]

D_v, D_h = transport_coefficients(m, net)

# How turbulent is this, really?
v_joint = np.median(np.abs(m.q_v)[net.link_v]) / DX
Re = 1000.0 * v_joint * 1.0e-4 / 1.0e-3      # rho v d / mu, d = 100 um aperture
print("flow regime")
print("  joint Darcy flux      %.3e m/s = %.2f m/yr" % (v_joint, v_joint * YR))
print("  Reynolds number       %.2e   -- laminar by five orders of magnitude," % Re)
print("                        so this is mechanical dispersion, not turbulence")
print()
print("transport coefficients [m2/s]")
print("  matrix links: molecular %.2e, dispersion %.2e"
      % (D_MOLECULAR / TORTUOSITY, ALPHA_L * np.median(np.abs(m.q_v)[~net.link_v]) / DX))
print("  joint links : molecular %.2e, dispersion %.2e"
      % (D_MOLECULAR, ALPHA_L * v_joint))
print()

lam = np.sqrt(D_MOLECULAR / TORTUOSITY / r)
print("predicted reaction-diffusion length  sqrt(D/r) = %.3f m = %.1f cells"
      % (lam, lam / DX))
print("  (real rindlets are about 2.5 cm)")
print()

t0 = time.perf_counter()
c_diff = solve_solute(m, net, D_v, D_h, r)
t_diff = time.perf_counter() - t0
c_adv = solve_solute(m, net, D_v * 0.0, D_h * 0.0, r)
print("solve cost: %.2f s for %d cells (once per step, matrix changes with M)"
      % (t_diff, NZ * NX))
print()

# Does the block interior stay saturated?
jc = np.nonzero(net.link_v[NZ // 2, :])[0]
mid_block = (jc[0] + jc[1]) // 2
row = NZ // 2
print("undersaturation (1-c) across one block at mid-depth:")
print("  distance from joint [cells]:", list(range(0, 16, 3)))
print("  advection only :", " ".join("%.2e" % (1 - c_adv[row, jc[0] + k])
                                     for k in range(0, 16, 3)))
print("  with diffusion :", " ".join("%.2e" % (1 - c_diff[row, jc[0] + k])
                                     for k in range(0, 16, 3)))
print()
print("fraction of the domain with (1-c) > 1e-6:")
print("  advection only %.1f %%   with diffusion %.1f %%"
      % (100 * ((1 - c_adv) > 1e-6).mean(), 100 * ((1 - c_diff) > 1e-6).mean()))

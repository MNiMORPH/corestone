#! /usr/bin/python3

"""
Dissolve granite along its joints, and see what is left.

The model is one equation. Dissolution runs at an Arrhenius rate constant
multiplied by how far the pore water is from equilibrium,

    R = k(T) * A * (1 - C / C_eq)

so water that has equilibrated stops weathering rock, however soluble the rock
and however warm the water. Fresh water descends the joints; the joints
therefore decide where weathering happens, and rock the water never reaches --
or reaches already saturated -- survives as a corestone.

Working in normalised concentration ``c = C / C_eq`` removes the need to assert
a solubility, and makes the integral over a cell of height ``dx`` exact:

    dc/dz = (1 - c) / L_eq       ->    c_out = 1 + (c_in - 1) * exp(-dx / L_eq)

with the **equilibration length**

    L_eq = q * C_eq / (k(T) * A)

the distance water travels before it is saturated. That single number carries
the model: rock further from a joint than ``L_eq`` never sees undersaturated
water. Raising the temperature *shrinks* ``L_eq``, so hotter water saturates
sooner and the weathering concentrates at the joints rather than spreading --
which is why warm does not mean weathered.

Flow is steady Darcy flow, solved for the hydraulic head,

    div( K grad H ) = 0,        H = psi - d        (d is depth, positive down)

with infiltration prescribed at the surface, a fixed head at the base and
no-flow sides. Conductance lives on the **links** between cells, which is where
the fracture network lives: a fractured link conducts, an intact one barely
does. Lateral flow along a subhorizontal joint is then not a special case -- it
is what the head field does when a low-resistance path exists.

An earlier version routed water down a gravity cascade instead, each cell
handing its water to the three cells below. That could enter a horizontal joint
but never travel along one, so the whole horizontal set was inert: deleting it
changed the weathering by 0.2 percentage points. The cascade was cheaper, but
the thing it could not represent was exactly the thing the joints are for.

Because the conductance does not evolve, the head is solved **once**. And
because a flow driven by the gradient of a potential cannot circulate, the flux
field is acyclic: solute is swept row by row, with the lateral exchange inside
each row solved as a tridiagonal system rather than iterated.

**Every parameter here is a placeholder.** None is measured. They are tabulated
in the design document, and no number from this module should be used as a
result.
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from scipy.linalg import solve_banded

#: Seconds in a Julian year. Time is seconds internally; years appear only at
#: the input and output edges, and this is where the conversion is named.
YEAR = 365.25 * 24 * 3600.0


class Weathering(object):
    """
    Weathering of a jointed rock section, driven by infiltrating rainwater.

    Takes a seeded :class:`~corestone.fractures.FractureNetwork` and evolves
    the fraction of the soluble mineral phase that has dissolved in each cell.
    """

    def __init__(self, network):
        self.network = network            # the seeded joint network
        self.dx = network.dx              # cell size [m]
        self.nz = network.nz              # rows, increasing downward
        self.nx = network.nx              # columns, increasing rightward

        # ---- parameters. ALL PLACEHOLDERS; see design/02-teaching-scope.md
        self.infiltration = 0.30 / YEAR   # recharge at the surface [m/s]
        self.k_fracture = 1.0e-5          # hydraulic conductivity, jointed
                                          # rock [m/s]        PLACEHOLDER
        self.k_matrix = 1.0e-8            # hydraulic conductivity, intact
                                          # granite [m/s]     PLACEHOLDER
        self.L_eq_ref = 0.50              # equilibration length at T_ref [m]
        self.T_ref = 285.0                # reference temperature [K]
        self.E_a = 60.0e3                 # activation energy [J/mol]
                                          # feldspar-ish, UNVERIFIED: this
                                          # needs Palandri & Kharaka (2004)
        self.R_gas = 8.314                # gas constant [J/mol/K]
        self.tau = 6700.0                 # M0/C_eq: volumes of saturated water
                                          # per volume of rock to dissolve it
        self.x_grus = 0.50                # soluble fraction lost -> grus
        self.x_core = 0.05                # below this, effectively unaltered
        self.f_inert = 0.30               # quartz: never dissolves, stays sand
        self.dt_max = 500.0 * YEAR        # step ceiling; see design/02
        self.dx_max = 0.02                # largest change in M per step

        # ---- state
        self.T = self.T_ref               # temperature [K]
        self.t = 0.0                      # model time [s]
        self.M = None                     # soluble mineral remaining, M/M0
        self.c = None                     # normalised concentration C/C_eq
                                          # LEAVING each cell, not entering
        self.H = None                     # hydraulic head [m]
        self.q = None                     # through-flux per cell [m2/s]
        self.q_v = None                   # flux on vertical links, down [m2/s]
        self.q_out_base = None            # flux leaving through the base
        self.q_h = None                   # flux on horizontal links, right
        self._in_above = None             # inflow from the row above
        self._in_left = None              # inflow from the left neighbour
        self._in_right = None             # inflow from the right neighbour

    # ---- parameter setters (one per parameter; units in the docstring)

    def set_temperature(self, value):
        """Temperature [K]. Higher shrinks the equilibration length."""
        self.T = value

    def set_infiltration(self, value):
        """Recharge at the ground surface [m/s]."""
        self.infiltration = value

    def set_equilibration_length(self, value):
        """Equilibration length at the reference temperature [m]."""
        self.L_eq_ref = value

    # ---- the physics

    @property
    def equilibration_length(self):
        """
        How far water travels before it is saturated [m], for water moving at
        the mean infiltration rate and through fresh rock.

        ``L_eq`` goes as ``1 / k(T)``, so it *shrinks* as the rock gets more
        reactive: a hotter system does less weathering per metre of flow path,
        because the water runs out of capacity sooner.

        This is the REFERENCE value. The local one scales with the local flux
        -- see :meth:`local_equilibration_length`.
        """
        return self.L_eq_ref * np.exp(
            (self.E_a / self.R_gas) * (1.0 / self.T - 1.0 / self.T_ref))

    def local_equilibration_length(self):
        """
        The equilibration length cell by cell [m].

            L_eq = q * C_eq / (k(T) * A)

        It is **proportional to the local flux**. Fast water in a joint travels
        far before it saturates; slow water in the matrix saturates almost at
        once. It also grows as the soluble mineral is consumed, because the
        reactive surface area falls with it.

        Getting this wrong is not a detail. Holding ``L_eq`` uniform makes the
        per-cell dissolution ``Q * beta * (1 - c)`` scale with ``Q``, so a joint
        carrying thirty times the flux dissolved thirty times faster per unit
        volume at the same undersaturation. The rate per unit volume is
        ``k A (1 - c)`` -- a property of the rock, not of how fast water moves
        past it. With ``L_eq`` proportional to ``Q`` the ``Q`` cancels and the
        rate is flux-independent, as it must be.
        """
        q_ref = self.infiltration * self.dx        # mean through-flux per cell
        return (self.equilibration_length * np.maximum(self.q, 1e-300) / q_ref
                / np.maximum(self.M, 1e-6))

    def solve_flow(self):
        """
        Steady Darcy head, and the link fluxes it implies.

        Finite volume on square cells, so the geometric factor is one: the flux
        along a link is ``K * (H_i - H_j)`` in m2/s per unit thickness. ``H`` is
        TOTAL head with elevation already in it -- adding a separate gravity
        term to the link flux double-counts it and manufactures water.

        Conductance is static, so this runs once.
        """
        nz, nx, dx = self.nz, self.nx, self.dx
        n = nz * nx
        idx = np.arange(n).reshape(nz, nx)

        kv = np.where(self.network.link_v, self.k_fracture, self.k_matrix)
        kh = np.where(self.network.link_h, self.k_fracture, self.k_matrix)

        rows, cols, vals = [], [], []
        pairs = [(idx[:-1, :].ravel(), idx[1:, :].ravel(), kv.ravel()),
                 (idx[:, :-1].ravel(), idx[:, 1:].ravel(), kh.ravel())]
        if self.network.periodic_x:
            # The wrap link closes the section onto itself, so there are no
            # side walls at all. A no-flow wall forces the lateral flow to
            # vanish there, which with subhorizontal joints manufactures a
            # domain-scale circulation and a drainage divide down the middle.
            kw = np.where(self.network.link_wrap, self.k_fracture,
                          self.k_matrix)
            pairs.append((idx[:, -1], idx[:, 0], kw))
        for a, b, k in pairs:
            rows += [a, a, b, b]
            cols += [a, b, b, a]
            vals += [k, -k, k, -k]
        A = sp.coo_matrix((np.concatenate(vals),
                           (np.concatenate(rows), np.concatenate(cols))),
                          shape=(n, n)).tocsr()

        # Infiltration into the top row [m2/s per unit thickness].
        b = np.zeros(n)
        b[idx[0, :]] = self.infiltration * dx

        # Base: the drainage boundary, psi = 0, so H = -depth. Applied as a
        # conductance to an external fixed head rather than by overwriting the
        # row. Overwriting pins the head but destroys continuity IN that row,
        # and the transport step then treats those cells as ordinary ones --
        # which cost about half a percent in the solute balance while the water
        # balance stayed exact, because water is solved and solute is swept.
        self._k_base = np.where(self.network.cell[-1, :],
                                self.k_fracture, self.k_matrix)
        self._h_base = -(nz - 0.5) * dx - 0.5 * dx
        A = A.tolil()
        A[idx[-1, :], idx[-1, :]] = A[idx[-1, :], idx[-1, :]] + self._k_base
        b[idx[-1, :]] += self._k_base * self._h_base

        H = spl.spsolve(A.tocsc(), b).reshape(nz, nx)
        self.H = H
        self.q_v = kv * (H[:-1, :] - H[1:, :])      # positive downward
        self.q_h = kh * (H[:, :-1] - H[:, 1:])      # positive rightward
        if self.network.periodic_x:
            kw = np.where(self.network.link_wrap, self.k_fracture,
                          self.k_matrix)
            self.q_wrap = kw * (H[:, -1] - H[:, 0])  # last column -> first
        else:
            self.q_wrap = np.zeros(nz)

        # Per-cell inflows, split by where they come from. Vertical flow is
        # strictly downward (checked in prototypes/probe_d_darcy.py), so the
        # only vertical inflow to a cell is from the row above it.
        self._in_above = np.zeros((nz, nx))
        self._in_above[0, :] = self.infiltration * dx
        self._in_above[1:, :] = np.maximum(self.q_v, 0.0)
        self._in_left = np.zeros((nz, nx))          # from the cell to the left
        self._in_left[:, 1:] = np.maximum(self.q_h, 0.0)
        self._in_right = np.zeros((nz, nx))         # from the cell to the right
        self._in_right[:, :-1] = np.maximum(-self.q_h, 0.0)
        # Across the seam, column nx-1 is "left of" column 0 and vice versa.
        self._in_left[:, 0] += np.maximum(self.q_wrap, 0.0)
        self._in_right[:, -1] += np.maximum(-self.q_wrap, 0.0)
        self.q = self._in_above + self._in_left + self._in_right
        # What leaves the domain through the base, per bottom-row cell.
        self.q_out_base = self._k_base * (H[-1, :] - self._h_base)
        return self.q

    # ---- lifecycle

    def initialize(self):
        """Set the rock fresh and route the flow. Run once, before update()."""
        self.M = np.ones((self.nz, self.nx))
        self.c = np.zeros((self.nz, self.nx))
        self.t = 0.0
        self.solve_flow()
        return self

    def update(self, dt=None):
        """
        Advance the rock state by one step, and return the step actually taken.

        The solute balance in a cell is steady -- transport is fast compared
        with the rock -- so what leaves equals what entered plus what dissolved:

            Q_i (1 + beta_i) c_i  -  sum(lateral inflow) c_j
                =  Q_i beta_i  +  (solute arriving from above)

        where ``Q_i`` is the through-flux and ``beta = expm1(dx / L_eq)``. That
        choice of beta is not an approximation: substituted into the balance it
        reproduces the exact exponential ``c_out = 1 + (c_in - 1) exp(-dx/L_eq)``
        in the one-dimensional case, while keeping the equation linear in the
        two-dimensional one.

        Vertical flow is everywhere downward, so rows can be swept in order.
        Within a row, water moving sideways along a joint couples neighbouring
        cells, which makes each row a tridiagonal system rather than a plain
        pass -- and that lateral coupling is the whole point of solving for the
        head instead of routing water downhill by rule.
        """
        # L_eq scales with the local flux and with the mineral left; see
        # local_equilibration_length().
        L_eq = self.local_equilibration_length()
        beta = np.where(self.q > 0.0, np.expm1(self.dx / np.maximum(L_eq, 1e-300)),
                        0.0)

        c = np.zeros((self.nz, self.nx))
        Q = self.q
        solute_above = np.zeros(self.nx)

        for iz in range(self.nz):
            q_i = Q[iz, :]
            live = q_i > 0.0

            diag = np.where(live, q_i * (1.0 + beta[iz, :]), 1.0)
            # Coefficients on the left and right neighbours' concentrations.
            lower = np.where(live, -self._in_left[iz, :], 0.0)
            upper = np.where(live, -self._in_right[iz, :], 0.0)
            rhs = np.where(live, q_i * beta[iz, :] + solute_above, 0.0)

            ab = np.zeros((3, self.nx))
            ab[0, 1:] = upper[:-1]        # superdiagonal
            ab[1, :] = diag
            ab[2, :-1] = lower[1:]        # subdiagonal

            if self.network.periodic_x:
                # The seam makes the row CYCLIC: cell 0 draws on cell nx-1 and
                # vice versa, which puts entries in the two far corners.
                # Sherman-Morrison folds them back into a banded solve.
                # NB: not named beta -- that is the reaction array above.
                corner_tr = -self._in_left[iz, 0] if live[0] else 0.0
                corner_bl = -self._in_right[iz, -1] if live[-1] else 0.0
                gamma = -ab[1, 0]
                ab[1, 0] -= gamma
                ab[1, -1] -= corner_tr * corner_bl / gamma
                u = np.zeros(self.nx); u[0], u[-1] = gamma, corner_tr
                y = solve_banded((1, 1), ab, rhs)
                z = solve_banded((1, 1), ab, u)
                vy = y[0] + corner_bl / gamma * y[-1]
                vz = z[0] + corner_bl / gamma * z[-1]
                c[iz, :] = np.clip(y - z * vy / (1.0 + vz), 0.0, 1.0)
            else:
                c[iz, :] = np.clip(solve_banded((1, 1), ab, rhs), 0.0, 1.0)

            if iz < self.nz - 1:
                solute_above = self._in_above[iz + 1, :] * c[iz, :]

        # What the water picked up is what the rock lost, per unit volume.
        rate = Q * beta * (1.0 - c) / (self.tau * self.dx * self.dx)
        step = min(dt if dt is not None else self.dt_max,
                   self.dx_max / max(rate.max(), 1e-30))
        self.M = np.clip(self.M - rate * step, 0.0, 1.0)
        self.c = c
        self.t += step
        return step

    def run(self, years):
        """Advance to ``years`` of model time, initializing if needed."""
        if self.M is None:
            self.initialize()
        target = years * YEAR
        while self.t < target:
            self.update()
        return self

    def finalize(self):
        """Nothing to release; present so the lifecycle is complete."""
        return self

    # ---- output

    @property
    def dissolved_fraction(self):
        """Fraction of the soluble phase that has dissolved, in [0, 1]."""
        return 1.0 - self.M

    @property
    def is_grus(self):
        """Cells that have lost enough of the soluble phase to fall apart."""
        return self.dissolved_fraction > self.x_grus

    @property
    def is_corestone(self):
        """Cells still effectively unaltered."""
        return self.dissolved_fraction < self.x_core

    @property
    def affinity(self):
        """The bracket, ``1 - C/C_eq``: how much capacity the water has left."""
        return 1.0 - self.c

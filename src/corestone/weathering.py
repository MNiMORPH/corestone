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
a solubility. The scale it sets is the **saturation length**

    saturation_length = q * C_eq / (k(T) * A)

the e-folding length of the approach to saturation -- *not* a distance at which
equilibrium is reached, because there is no equilibrium here. ``c`` approaches
1 asymptotically and never arrives.

Solute moves by advection **and by diffusion**:

    div(q c) - div(D grad c) = r (1 - c),      r = k A / C_eq

Without the diffusive term, rock off a flow path never weathers at all: a block
interior saturates and then sits at ``c = 1`` for ever, and the model gives
joints entirely dissolved beside blocks entirely untouched with nothing in
between. Diffusive export of solute toward a flushed joint keeps the interior
undersaturated, which is what forms a weathering rind -- and, because a corner
sheds solute to two faces and a face to one, what rounds corestones.

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

#: Seconds in a Julian year. Time is seconds internally; years appear only at
#: the input and output edges, and this is where the conversion is named.
YEAR = 365.25 * 24 * 3600.0

#: Column ordering for every sparse factorisation here. Both matrices -- the
#: conductance Laplacian for the head and the transport operator for the solute
#: -- are STRUCTURALLY SYMMETRIC, because every link contributes both (i, j) and
#: (j, i). That is the precondition for minimum degree on ``A + A.T``, and it is
#: asserted in the tests rather than assumed, since an upwind stencil that
#: reached only one way would silently break it. SuperLU's default is COLAMD,
#: which orders for ``A.T A`` and is the right choice only when the structure is
#: unsymmetric; here it roughly doubles the fill. Measured on the 22,650-cell
#: section: nnz(L + U) 2,070,424 -> 1,060,262, factorisation 302 -> 170 ms,
#: back-substitution 3.16 -> 2.74 ms.
ORDERING = "MMD_AT_PLUS_A"


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
        self.L_ref = 0.50                 # saturation length at T_ref, mean
                                          # infiltration, fresh rock [m]
        self.T_ref = 285.0                # reference temperature [K]
        self.E_a = 60.0e3                 # activation energy [J/mol]
                                          # feldspar-ish, UNVERIFIED: needs
                                          # Palandri & Kharaka (2004)
        self.delta_H_r = 25.0e3           # enthalpy of the dissolution
                                          # reaction [J/mol], van 't Hoff.
                                          # UNVERIFIED. Only the DIFFERENCE
                                          # (E_a - delta_H_r) sets the length
                                          # scale, so the two must not be
                                          # picked independently.
        self.R_gas = 8.314                # gas constant [J/mol/K]
        self.tau_ref = 6700.0             # M0/C_eq at T_ref: volumes of
                                          # saturated water per volume of rock
        self.D_molecular = 1.0e-9         # aqueous diffusivity [m2/s]
        self.tortuosity = 10.0            # matrix tortuosity [-]
        self.dispersivity = 0.05          # longitudinal dispersivity [m]
        self.x_grus = 0.50                # soluble fraction lost -> grus
        self.x_core = 0.05                # below this, effectively unaltered
        self.f_inert = 0.30               # quartz: never dissolves, stays sand
        self.dt_max = 2000.0 * YEAR       # step ceiling. 500 yr cost 4x the
        self.dx_max = 0.05                # steps for 0.8 % in the answer, back
                                          # when a step was cheap; see design/07
        self.krylov_tol = 1.0e-10         # convergence of the reused solve
        self.max_krylov_iterations = 15   # past this, refactorise instead

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
        self._T = None                    # constant part of the solute operator
        self._T_key = None                # the values it was built from
        self._A = None                    # the step matrix, reused in place
        self._diag = None                 # where its diagonal sits in .data
        self._lu = None                   # cached factorisation, reused
        self._x = None                    # last solute solution, unclipped
        self.factorisations = 0           # how many times it was rebuilt

    # ---- parameter setters (one per parameter; units in the docstring)

    def set_temperature(self, value):
        """Temperature [K]. Higher shrinks the equilibration length."""
        self.T = value

    def set_infiltration(self, value):
        """Recharge at the ground surface [m/s]."""
        self.infiltration = value

    def set_saturation_length(self, value):
        """Saturation length at the reference temperature [m]."""
        self.L_ref = value

    # ---- the physics

    @property
    def rate_factor(self):
        """``k(T) / k_ref``, Arrhenius."""
        return np.exp(-(self.E_a / self.R_gas)
                      * (1.0 / self.T - 1.0 / self.T_ref))

    @property
    def solubility_factor(self):
        """
        ``C_eq(T) / C_eq_ref``, van 't Hoff.

        Solubility is temperature dependent too, and treating it as constant
        was not a small error: in the transport-limited regime -- which is most
        of this model -- the amount dissolved scales with ``C_eq`` and not with
        the rate constant at all, so holding it fixed removed the half of the
        temperature dependence that dominates.
        """
        return np.exp(-(self.delta_H_r / self.R_gas)
                      * (1.0 / self.T - 1.0 / self.T_ref))

    @property
    def saturation_length(self):
        """
        The e-folding length of the approach to saturation [m], for water at
        the mean infiltration rate through fresh rock.

            saturation_length = q * C_eq / (k(T) * A)

        **Not** a distance at which equilibrium is reached. There is no
        equilibrium: ``c`` approaches 1 asymptotically, and after n of these
        lengths the undersaturation is ``exp(-n)`` of what it was. The system
        has an asymptote, and the model has a normalisation rather than a
        thermodynamics -- ``C_eq`` never appears alone.

        Because it goes as ``C_eq / k``, its temperature dependence is set by
        ``(E_a - delta_H_r)``, not by ``E_a`` alone.
        """
        return self.L_ref * self.solubility_factor / self.rate_factor

    @property
    def specific_reaction_coefficient(self):
        """
        ``r / M`` [1/s]: the reaction coefficient per unit mineral remaining.

        A SCALAR, and that is the point of naming it. ``r`` falls with the
        soluble mineral because the reactive surface area does, and it falls
        exactly in proportion -- so the quotient carries no ``M`` at all, and
        the ratio can be formed without ever dividing by a mineral content that
        is on its way to zero. It is what lets the rock be integrated exactly
        rather than tangentially; see :meth:`update`.
        """
        return ((self.infiltration / self.L_ref)
                * self.rate_factor / self.solubility_factor)

    @property
    def reaction_coefficient(self):
        """
        ``r = k A / C_eq`` [1/s]: the rate at which undersaturation is consumed.

        This is the flux-independent form. Since ``r = q / saturation_length``
        and the saturation length is itself proportional to ``q``, the flux
        cancels **explicitly** here rather than implicitly inside an exponent.
        It falls with the soluble mineral remaining, because the reactive
        surface area does.
        """
        return self.specific_reaction_coefficient * np.maximum(self.M, 0.0)

    @property
    def tau(self):
        """
        ``M0 / C_eq``: volumes of saturated water needed per volume of rock.

        Falls as solubility rises, so a warmer and more soluble fluid carries
        more away per unit volume. This is the second place ``C_eq`` enters.
        """
        return self.tau_ref / self.solubility_factor

    def local_saturation_length(self):
        """The saturation length cell by cell [m]: ``q / r``."""
        r = self.reaction_coefficient
        return (self.q / self.dx) / np.maximum(r, 1e-300)

    def transport_coefficients(self):
        """
        Solute transport coefficient on every link [m2/s].

            D = D_molecular / tortuosity  +  dispersivity * |v|

        Molecular diffusion plus hydrodynamic (mechanical) dispersion. The
        second is sometimes called turbulent diffusion and the operator is the
        same, but at these fluxes the Reynolds number is around 3e-5 -- laminar
        by five orders of magnitude. The distinction matters only because it is
        pore and aperture geometry that sets the dispersivity, not eddies.

        Without this term the model is pure advection, and rock that is not on
        a flow path never weathers at all: a block interior saturates and then
        sits at ``c = 1`` for ever. Diffusive export of solute toward a flushed
        joint is what keeps the interior undersaturated, and is therefore what
        lets a weathering rind form. It is also the geometric route to
        spheroidal rounding, since a corner sheds solute to two faces and an
        edge to one.
        """
        dm_v = np.where(self.network.link_v, self.D_molecular,
                        self.D_molecular / self.tortuosity)
        dm_h = np.where(self.network.link_h, self.D_molecular,
                        self.D_molecular / self.tortuosity)
        return (dm_v + self.dispersivity * np.abs(self.q_v) / self.dx,
                dm_h + self.dispersivity * np.abs(self.q_h) / self.dx)

    def _transport_operator(self):
        """
        The part of the solute operator that never changes: advection,
        diffusion and the outflow through the base.

        Built once. The flow field is static, so only the reaction term on the
        diagonal moves from step to step -- which is what makes a cached
        factorisation worth keeping.
        """
        # Cache on the values that built it. A bare `is not None` check made
        # changing D_molecular or the dispersivity silently do nothing, which a
        # test caught -- the same shape of defect as a docstring drifting from
        # its code, and one a cache invites.
        key = (self.D_molecular, self.tortuosity, self.dispersivity,
               id(self.q_v), float(np.sum(self.q_v)), float(np.sum(self.q_h)))
        if self._T is not None and self._T_key == key:
            return self._T
        self._lu = None                    # the factorisation went with it
        nz, nx, dx = self.nz, self.nx, self.dx
        idx = np.arange(nz * nx).reshape(nz, nx)
        D_v, D_h = self.transport_coefficients()
        rows, cols, vals = [], [], []

        def add(i, j, v):
            rows.append(np.asarray(i).ravel())
            cols.append(np.asarray(j).ravel())
            vals.append(np.asarray(v, dtype=float).ravel()
                        * np.ones(np.asarray(i).size))

        for a_, b_, f in ((idx[:-1, :], idx[1:, :], self.q_v),
                          (idx[:, :-1], idx[:, 1:], self.q_h)):
            fwd, rev = np.maximum(f, 0.0), np.maximum(-f, 0.0)
            add(a_, a_, fwd);  add(b_, a_, -fwd)
            add(b_, b_, rev);  add(a_, b_, -rev)
        for a_, b_, D in ((idx[:-1, :], idx[1:, :], D_v),
                          (idx[:, :-1], idx[:, 1:], D_h)):
            add(a_, a_, D);  add(a_, b_, -D)
            add(b_, b_, D);  add(b_, a_, -D)
        if self.network.periodic_x:
            fwd = np.maximum(self.q_wrap, 0.0)
            rev = np.maximum(-self.q_wrap, 0.0)
            Dw = np.where(self.network.link_wrap, self.D_molecular,
                          self.D_molecular / self.tortuosity) \
                 + self.dispersivity * np.abs(self.q_wrap) / dx
            l, rgt = idx[:, -1], idx[:, 0]
            add(l, l, fwd + Dw);     add(rgt, l, -(fwd + Dw))
            add(rgt, rgt, rev + Dw); add(l, rgt, -(rev + Dw))
        add(idx[-1, :], idx[-1, :], self.q_out_base)

        self._T = sp.coo_matrix(
            (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
            shape=(nz * nx, nz * nx)).tocsc()
        self._T_key = key

        # Where the diagonal sits inside the CSC value array, so the reaction
        # term can be written straight into it every step instead of building a
        # second matrix and adding. Every cell has a diagonal entry -- every
        # cell has at least one link, and a link puts a coefficient on both of
        # its cells' diagonals -- and that is asserted here rather than trusted,
        # because a missing one would silently drop that cell's reaction.
        n = nz * nx
        col_of = np.repeat(np.arange(n), np.diff(self._T.indptr))
        self._diag = np.nonzero(self._T.indices == col_of)[0]
        if self._diag.size != n:
            raise AssertionError("%d of %d cells have no diagonal entry"
                                 % (n - self._diag.size, n))
        self._A = self._T.copy()
        return self._T

    def _step_matrix(self, r):
        """
        The operator for one step: transport, plus the reaction on the diagonal.

        Only the diagonal moves from step to step, so this writes into a matrix
        that is allocated once rather than adding two sparse matrices and
        converting the result. The values are identical -- the same two floats
        added in the same order -- and it costs 0.105 ms in place against
        1.025 ms via ``sp.diags`` at 22,650 cells.
        """
        T = self._transport_operator()
        A = self._A
        A.data[:] = T.data
        A.data[self._diag] += np.broadcast_to(
            r * self.dx * self.dx, (self.nz, self.nx)).ravel()
        return A

    def solve_solute(self, r):
        """
        Steady advection-diffusion-reaction for the normalised concentration.

            sum_out f c_i - sum_in f c_j + sum_links D (c_i - c_j)
                + r dx^2 c_i  =  r dx^2

        One sparse solve for the whole field. Diffusion is not directional, so
        the row-by-row sweep that pure advection allowed is gone -- and with it
        the cyclic Sherman-Morrison seam solve and the ``Q * beta * (1 - c)``
        pickup, which cancelled catastrophically once ``beta`` grew past about
        1e5. The reaction term here is linear and bounded, so neither hazard
        remains.

        The solve STARTS FROM THE PREVIOUS STEP'S FIELD. The rock changes by at
        most ``dx_max`` of its mineral content in a step, so the concentration
        moves very little, and the old field is a far better guess than the
        zero the solver would otherwise assume. This costs no accuracy: the
        convergence test is on the residual ``||b - A x||``, which does not
        know where the iteration started.

        Water entering at the surface carries ``c = 0``, so it contributes
        nothing to the inflow sum: undersaturation enters only through the
        reaction term on the right.

        Cost. Only the diagonal changes between steps, so the LU factorisation
        is cached and reused as a preconditioner: at 22,650 cells a
        back-substitution is 2.7 ms against 170 ms to refactorise.

        WHY IT IS REFRESHED, corrected. This docstring used to say the
        preconditioner is rebuilt when the iteration count says so, and called
        that self-tuning. It is not what happens. Raising
        ``max_krylov_iterations`` from 15 to a million changes nothing at all:
        the same 15 factorisations, the same answer to the last bit. Every
        refresh but the first is triggered by ``info = -10`` -- a BiCGSTAB
        BREAKDOWN, after two iterations, on 14 of 108 steps at dx = 0.02. The
        algorithm divides by a quantity that has gone to zero, which happens
        here because the preconditioner is very good rather than because it is
        stale. The cap is a backstop that has never fired at its current value.

        Lowering the cap so that it DOES fire is worth about 15 %, and that is
        an untaken decision rather than an oversight: it is a tuned number,
        the optimum moves between 4 and 5 from case to case, and it buys
        little. Measured, min of three alternating runs, against cap = 15:

            cap     app 3 m   fine 3 m   45 deg   wide
              4       x1.28      x0.86    x1.09   x1.38
              5       x1.28      x1.14    x1.04   x1.15
              8       x1.17      x0.87    x0.94   x1.04

        Defect correction -- ``x <- x + LU^-1 (b - A x)``, which cannot break
        down and needs one back-substitution per iteration instead of two --
        was implemented and measured, and it is a WASH: x1.06, x0.97, x1.17 on
        the three cases, alternated in one process. It is not here because it
        was not better, and the 3.2x that first appeared to favour it was
        machine load drifting between sequential runs, not the solver.
        """
        dx = self.dx
        A = self._step_matrix(r)
        b = np.broadcast_to(r * dx * dx, (self.nz, self.nx)).ravel().copy()

        def direct():
            self._lu = spl.splu(A, permc_spec=ORDERING)
            self.factorisations += 1
            return self._lu.solve(b)

        if self._lu is None:
            x = direct()
        else:
            n_it = [0]
            P = spl.LinearOperator(A.shape, matvec=self._lu.solve)
            x, info = spl.bicgstab(
                A, b, x0=self._x, M=P, atol=0.0, tol=self.krylov_tol,
                callback=lambda xk: n_it.__setitem__(0, n_it[0] + 1))
            if info != 0 or n_it[0] > self.max_krylov_iterations:
                x = direct()                       # preconditioner has gone stale
        # Kept unclipped, and kept whether or not it was clipped on the way
        # out, because it is the next step's starting guess and not an answer.
        self._x = x
        return np.clip(x, 0.0, 1.0).reshape(self.nz, self.nx)

    def flow_operator(self):
        """
        The conductance matrix for the head, and the right-hand side.

        Assembled entirely as triplets, base boundary included. COO sums
        duplicate entries, so the boundary conductance is simply one more
        contribution to the diagonal -- there is no need to reach into an
        assembled matrix and modify a row, which needed a LIL round trip and
        cost 44 ms of a 180 ms initialize() at 22,650 cells.

        Returned rather than kept private so that a test can look at it: it is
        the second matrix ``ORDERING`` claims is structurally symmetric.
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
        for a_, b_, k in pairs:
            rows += [a_, a_, b_, b_]
            cols += [a_, b_, b_, a_]
            vals += [k, -k, k, -k]

        # Infiltration into the top row [m2/s per unit thickness].
        rhs = np.zeros(n)
        rhs[idx[0, :]] = self.infiltration * dx

        # Base: the drainage boundary, psi = 0, so H = -depth. Applied as a
        # conductance to an external fixed head rather than by overwriting the
        # row. Overwriting pins the head but destroys continuity IN that row,
        # and the transport step then treats those cells as ordinary ones --
        # which cost about half a percent in the solute balance while the water
        # balance stayed exact, because water is solved and solute is swept.
        self._k_base = np.where(self.network.cell[-1, :],
                                self.k_fracture, self.k_matrix)
        self._h_base = -(nz - 0.5) * dx - 0.5 * dx
        rows.append(idx[-1, :]); cols.append(idx[-1, :])
        vals.append(self._k_base)
        rhs[idx[-1, :]] += self._k_base * self._h_base

        A = sp.coo_matrix((np.concatenate(vals),
                           (np.concatenate(rows), np.concatenate(cols))),
                          shape=(n, n)).tocsc()
        return A, rhs

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
        kv = np.where(self.network.link_v, self.k_fracture, self.k_matrix)
        kh = np.where(self.network.link_h, self.k_fracture, self.k_matrix)

        A, b = self.flow_operator()
        H = spl.splu(A, permc_spec=ORDERING).solve(b).reshape(nz, nx)
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
        self._T = None
        self._T_key = None
        self._A = None
        self._diag = None
        self._lu = None
        self._x = None
        self.factorisations = 0
        self.M = np.ones((self.nz, self.nx))
        self.c = np.zeros((self.nz, self.nx))
        self.t = 0.0
        self.solve_flow()
        return self

    def update(self, dt=None):
        """
        Advance the rock state by one step, and return the step actually taken.

        Solute transport is solved at steady state -- it is fast next to the
        rock changing -- so there is no storage term and the model carries no
        porosity, and therefore no residence time.

        What the rock loses is what the water gains:

            d(M/M0)/dt = - r (1 - c) / tau

        and ``r`` is proportional to ``M``, because the reactive surface area
        is. Over a step in which ``c`` is held the equation is therefore linear
        in ``M``, and it integrates EXACTLY:

            M(t + dt) = M(t) exp(-lambda dt), lambda = (r / M) (1 - c) / tau

        Forward Euler stood here before, taking the tangent to that
        exponential. The tangent always undershoots, which is the only reason
        the step ever needed clipping at zero, and it was six times less
        accurate at identical cost -- the same solves, one np.exp added.
        ``lambda`` is formed from :attr:`specific_reaction_coefficient` and so
        never divides by a mineral content approaching zero.

        The step is still limited, because ``c`` is held across it while the
        rock underneath it moves. ``dx_max`` bounds that drift the same way it
        did before -- on the tangent, ``lambda M dt``, which is the old Euler
        rate exactly -- so the parameter still means what it used to mean.
        """
        c = self.solve_solute(self.reaction_coefficient)
        lam = self.specific_reaction_coefficient * (1.0 - c) / self.tau
        rate = lam * self.M                        # d(M/M0)/dt [1/s]
        step = min(dt if dt is not None else self.dt_max,
                   self.dx_max / max(rate.max(), 1e-30))
        # The clip is a guard, not a mechanism: lambda >= 0 because c <= 1, so
        # the exponential cannot leave (0, 1] on its own.
        self.M = np.clip(self.M * np.exp(-lam * step), 0.0, 1.0)
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

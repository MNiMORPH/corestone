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

Flow is steady gravity-driven descent: one sweep from the surface down, each
cell handing its water to the three cells below split by their conductance.
Gravity makes the grid a one-way cascade ordered by depth, so there is no
pressure solve. See ``design/02-teaching-scope.md`` for what that buys and what
it costs.

**Every parameter here is a placeholder.** None is measured. They are tabulated
in the design document, and no number from this module should be used as a
result.
"""

import numpy as np

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
        self.k_fracture = 1000.0          # routing conductance, joint cell
        self.k_matrix = 1.0               # routing conductance, intact rock
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
        self.c = None                     # normalised concentration, C/C_eq
        self.q = None                     # water flux per cell [m2/s]
        self._K = None                    # routing conductance per cell

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
        How far water travels before it is saturated, at the current
        temperature [m].

        ``L_eq`` goes as ``1 / k(T)``, so it *shrinks* as the rock gets more
        reactive. That is the counterintuitive part: a hotter system does less
        weathering per metre of flow path, not more, because the water runs out
        of capacity sooner.
        """
        return self.L_eq_ref * np.exp(
            (self.E_a / self.R_gas) * (1.0 / self.T - 1.0 / self.T_ref))

    def route_flow(self):
        """
        Steady gravity-driven descent, in one sweep from the surface down.

        Each cell hands its water to the three cells below, split by their
        conductance. Gravity orders the grid by depth, so this is a single pass
        and no pressure solve is needed. Conductance does not change as the
        rock weathers, so this is computed once.
        """
        self._K = np.where(self.network.cell, self.k_fracture, self.k_matrix)
        q = np.zeros((self.nz, self.nx))
        q[0, :] = self.infiltration * self.dx        # m2/s per unit thickness
        for iz in range(self.nz - 1):
            f_l, f_c, f_r = self._split(iz + 1)
            send = q[iz, :]
            q[iz + 1, :] += send * f_c
            q[iz + 1, :-1] += (send * f_l)[1:]
            q[iz + 1, 1:] += (send * f_r)[:-1]
        self.q = q
        return q

    def _split(self, iz):
        """Fractions of a cell's water going down-left, down and down-right."""
        below = self._K[iz, :]
        wl = np.concatenate([[0.0], below[:-1]])
        wr = np.concatenate([below[1:], [0.0]])
        tot = wl + below + wr
        return wl / tot, below / tot, wr / tot

    # ---- lifecycle

    def initialize(self):
        """Set the rock fresh and route the flow. Run once, before update()."""
        self.M = np.ones((self.nz, self.nx))
        self.c = np.zeros((self.nz, self.nx))
        self.t = 0.0
        self.route_flow()
        return self

    def update(self, dt=None):
        """
        Advance the rock state by one step, and return the step actually taken.

        The solute sweep follows the flow: water enters the top fresh, and each
        cell's outlet concentration is the exact integral of the rate law over
        the cell. What the water picked up is what the rock lost.
        """
        # Surface area falls with the mineral that is left, so L_eq grows.
        L_eq = self.equilibration_length / np.maximum(self.M, 1e-6)

        c_in = np.zeros((self.nz, self.nx))
        dissolved = np.zeros((self.nz, self.nx))
        carry_q = np.zeros(self.nx)
        carry_qc = np.zeros(self.nx)

        for iz in range(self.nz):
            qi = self.q[iz, :]
            if iz == 0:
                ci = np.zeros(self.nx)
            else:
                ci = np.where(qi > 0.0,
                              carry_qc / np.maximum(carry_q, 1e-300), 0.0)
            c_out = 1.0 + (ci - 1.0) * np.exp(-self.dx / L_eq[iz, :])
            dissolved[iz, :] = qi * (c_out - ci) / self.dx
            c_in[iz, :] = ci
            if iz == self.nz - 1:
                break
            carry_q = np.zeros(self.nx)
            carry_qc = np.zeros(self.nx)
            for w, shift in zip(self._split(iz + 1), (-1, 0, +1)):
                add_q, add_qc = qi * w, qi * w * c_out
                if shift == 0:
                    carry_q += add_q
                    carry_qc += add_qc
                elif shift == -1:
                    carry_q[:-1] += add_q[1:]
                    carry_qc[:-1] += add_qc[1:]
                else:
                    carry_q[1:] += add_q[:-1]
                    carry_qc[1:] += add_qc[:-1]

        rate = dissolved / (self.tau * self.dx)        # d(M/M0)/dt [1/s]
        step = min(dt if dt is not None else self.dt_max,
                   self.dx_max / max(rate.max(), 1e-30))
        self.M = np.clip(self.M - rate * step, 0.0, 1.0)
        self.c = c_in
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

"""
Probe L: should infiltration be prescribed, or limited by what the rock can take?

OPEN. Written 2026-09-22 after Andy asked whether the surface forcing should be
a constant head rather than a constant flux, so that water is not forced
through a rock that cannot pass it.

THE DEFECT, MEASURED FIRST. The model prescribes a flux on every surface cell,
``rhs[idx[0, :]] = infiltration * dx``. Total inflow is therefore the same
whatever the rock is like:

    joint spacing   total inflow
      0.5 m         2.852e-08 m2/s
      1.0 m         2.852e-08
      2.0 m         2.852e-08
      3.0 m         2.852e-08

Fracture the granite or seal it, the same rain goes in. So the joint-spacing
slider changes only WHERE water goes and never HOW MUCH, which removes one of
the two things fractures do.

With no fractures at all it stops being a limitation and becomes nonsense.
Intact matrix K_sat is 0.0128 m/yr against a prescribed 0.30, so Darcy needs a
gradient of 23 where gravity alone gives 1, and the solver duly reports +67 m
of head at the land surface -- water standing 67 m deep on the outcrop. The
model does not fail. It returns that, quietly.

THE FIX, WHICH IS STANDARD. Rain falls at a rate; the rock accepts it up to its
infiltration capacity; the rest is runoff. As a boundary condition that is a
SWITCH, cell by cell:

    prescribed flux   while the head it demands stays below the land surface
    prescribed head   h = 0, once it does not

A cell in the second state has ponded, and the water it cannot take runs off.
The set of ponded cells is not known before the solve, so the flow problem
becomes nonlinear and is iterated: solve, look for surface cells above h = 0,
switch those to fixed head, solve again, and also release any fixed-head cell
that is trying to take MORE than the rainfall. Two or three passes settle it.

This is the ordinary atmospheric boundary condition of unsaturated-flow
modelling, not an invention, and it is the same conductance-to-an-external-head
device this model already uses at its base.

ANSWERED, 2026-09-22. Three things, and the third is the one that decides
whether this is worth building.

1.  IT NEVER BINDS AT THE SETTINGS THE DEMO OFFERS. Zero ponded cells at every
    joint spacing from 0.3 to 3 m, 100 % of the rain taken, and the highest
    surface head anywhere is -2.85 m -- nearly three metres below the ground.
    So the change costs nothing at present settings and every measured number
    on the page stands.

2.  IT FIXES THE NO-FRACTURE CASE EXACTLY. Inflow falls from the prescribed
    0.30 m/yr to 0.0126, against an intact matrix K_sat of 0.0128 -- which is
    the check that the physics is right, because ponded water over a
    gravity-driven column gives q -> K_sat at unit gradient. Surface head goes
    from +67.33 m to -0.05. Runoff, which the model could not previously
    represent at all, is 0.2874 m/yr: 96 % of the rain.

    Water is conserved to 1.8e-13 with the constraint active.

3.  BUT THE TRANSITION IS BINARY, NOT GRADED. Ponded cells against rainfall
    and joint spacing:

        rain \ spacing   0.3 m   1.0 m   3.0 m   none
          0.30 m/yr         0       0       0      60
          1.00 m/yr         0       0       0      60
          3.00 m/yr         0       0       0      60
         10.00 m/yr         0       0      23      60

    A single joint has enormous capacity -- K_sat_fracture is 30,000 times the
    matrix -- so ANY joint takes all the rain the demo can offer, and the
    slider's effect on inflow is off or on rather than a curve. The demo tops
    out at 1.00 m/yr, so even the 3 m-spacing, 10 m/yr corner is out of reach.

WHAT THAT MEANS. The constraint is worth having, and for one reason rather
than three: it makes the no-fracture endpoint POSSIBLE. Without it that
endpoint returns 67 m of standing water and the model has to refuse the
setting. With it, the endpoint is the physically interesting case Andy asked
for -- little water, all of it near the surface, a shallow weathering horizon
-- and it arrives as a consequence of Darcy rather than as a special rule.

It is NOT worth having as a correction to the fractured cases, because there
it does nothing at all.

Run:  PYTHONPATH=src python3 prototypes/probe_l_ponding_boundary.py
"""

import numpy as np
import scipy.sparse as sp

from corestone import FractureNetwork, Weathering, orthogonal_grid, YEAR


class Ponding(Weathering):
    """Infiltration capped by what the rock will take; the rest runs off."""

    #: Head at which a surface cell has ponded. Zero is the land surface.
    h_pond = 0.0

    #: Most passes of the switching iteration. It converges in two or three;
    #: this is a backstop, not a tuning parameter.
    max_ponding_passes = 12

    def __init__(self, network):
        Weathering.__init__(self, network)
        self.ponded = None          # boolean mask over surface cells
        self.runoff = 0.0           # [m2/s] rain the rock would not take
        self.ponding_passes = 0

    def flow_operator(self, ponded=None):
        """The base operator, with ponded surface cells switched to fixed head.

        A flux cell keeps ``rhs = q dx``. A ponded cell instead gets a
        conductance to an external head at the land surface, exactly as the
        base boundary does -- which pins the head without overwriting the row
        and so keeps continuity intact.
        """
        # Default to the CONVERGED mask, because the base class's solve_flow
        # calls this with no argument and would otherwise rebuild the
        # unconstrained operator and discard the whole iteration.
        if ponded is None:
            ponded = self.ponded
        A, rhs = Weathering.flow_operator(self)
        if ponded is None or not ponded.any():
            return A, rhs
        idx = np.arange(self.nz * self.nx).reshape(self.nz, self.nx)
        top = idx[0, :][ponded]
        K_sat_top = np.where(self.network.cell[0, :], self.K_sat_fracture,
                         self.k_matrix_at_T)[ponded]
        rhs[top] = K_sat_top * self.h_pond          # replaces the prescribed flux
        # Add the conductance on the diagonal as a separate COO term rather
        # than indexing into the assembled matrix: scipy will not add a scalar
        # into a sparse slice, and summing duplicates is how the rest of this
        # operator is built anyway.
        n = self.nz * self.nx
        extra = sp.coo_matrix((K_sat_top, (top, top)), shape=(n, n))
        return (A + extra).tocsc(), rhs

    def solve_flow(self):
        """Solve, switching ponded cells, until the set of them stops moving."""
        ponded = np.zeros(self.nx, dtype=bool)
        for n in range(self.max_ponding_passes):
            self._flow_lu = None                # the operator changed
            self._H_prev = None
            A, b = self.flow_operator(ponded)
            H = self._solve_head(A, b).reshape(self.nz, self.nx)

            # a flux cell that has risen above the surface must pond...
            add = (~ponded) & (H[0, :] > self.h_pond)
            # ...and a ponded cell drawing more than the rain must be released
            K_sat_top = np.where(self.network.cell[0, :], self.K_sat_fracture,
                             self.k_matrix_at_T)
            taking = K_sat_top * (self.h_pond - H[0, :])
            drop = ponded & (taking > self.rainfall * self.dx)
            if not add.any() and not drop.any():
                break
            ponded = (ponded | add) & ~drop
        self.ponded = ponded
        self.ponding_passes = n + 1

        # The base class now does the rest -- link fluxes, per-cell inflows,
        # the base outflow -- and its call to flow_operator() picks the
        # converged mask up from self.ponded.
        self._flow_lu = None
        self._H_prev = None
        out = Weathering.solve_flow(self)

        # The base class sets the surface inflow to the PRESCRIBED rate, which
        # is what a ponded cell no longer takes. Recompute it from the head it
        # actually reached, and rebuild the per-cell inflow that depends on it.
        if self.ponded is not None and self.ponded.any():
            K_sat_top = np.where(self.network.cell[0, :], self.K_sat_fracture,
                             self.k_matrix_at_T)
            taken = np.where(self.ponded,
                             K_sat_top * (self.h_pond - self.H[0, :]),
                             self.rainfall * self.dx)
            self._in_above[0, :] = np.maximum(taken, 0.0)
            self.q = self._in_above + self._in_left + self._in_right

        rain = self.rainfall * self.dx * self.nx
        took = float(self._in_above[0, :].sum())
        self.runoff = rain - took
        return out


def build(cls, spacing=1.0, fractured=True, q=0.30, tC=11.85):
    net = FractureNetwork(60, 60, 0.05, periodic_x=True)
    if fractured:
        net = net.seed(sets=orthogonal_grid(spacing),
                       rng=np.random.default_rng(12345))
    else:
        net.link_v = np.zeros((59, 60), dtype=bool)
        net.link_h = np.zeros((60, 59), dtype=bool)
        net.link_wrap = np.zeros(60, dtype=bool)
        net.cell = np.zeros((60, 60), dtype=bool)
    m = cls(net)
    m.set_rainfall(q / YEAR)
    m.set_temperature(tC + 273.15)
    m.initialize()
    return m


if __name__ == "__main__":
    rain = 0.30

    print("1. Does the constraint bind at the settings the demo offers?\n")
    print("   spacing   ponded cells   inflow [m/yr]   of %.2f   h_top max"
          % rain)
    for spac in (0.3, 0.5, 1.0, 2.0, 3.0):
        m = build(Ponding, spacing=spac)
        took = float(m._in_above[0, :].sum()) / m.dx / m.nx * YEAR
        print("    %.1f m      %3d of %d      %.4f        %5.1f %%   %+7.2f m"
              % (spac, m.ponded.sum(), m.nx, took, 100 * took / rain,
                 m.H[0].max()))

    print("\n2. And with no fractures at all:\n")
    for cls, lab in ((Weathering, "as shipped "), (Ponding, "with ponding")):
        m = build(cls, fractured=False)
        took = float(m._in_above[0, :].sum()) / m.dx / m.nx * YEAR
        extra = ""
        if isinstance(m, Ponding):
            extra = "   ponded %d of %d, %d passes" % (
                m.ponded.sum(), m.nx, m.ponding_passes)
        print("   %s  inflow %.4f m/yr   h_top %+8.2f m%s"
              % (lab, took, m.H[0].max(), extra))

    print("\n3. Does inflow now fall as the rock seals up?\n")
    print("   spacing   inflow [m/yr]   runoff [m/yr]")
    for spac in (0.3, 1.0, 3.0):
        m = build(Ponding, spacing=spac)
        took = float(m._in_above[0, :].sum()) / m.dx / m.nx * YEAR
        print("    %.1f m      %.4f          %.4f" % (spac, took, rain - took))
    m = build(Ponding, fractured=False)
    took = float(m._in_above[0, :].sum()) / m.dx / m.nx * YEAR
    print("    none       %.4f          %.4f" % (took, rain - took))

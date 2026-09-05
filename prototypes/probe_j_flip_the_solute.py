"""
Probe J: flip the solute from product to reactant, and look at it.

OPEN. Step 2 of design/08-BUILD.md, written as a subclass so that the picture
can be judged before 200 lines of docstring are rewritten around it.

WHAT CHANGES, and it is three lines of arithmetic:

    now   div(q c) - div(D grad c) + r c = r          inlet c = 0, c -> 1
    then  div(q c) - div(D grad c) + r c = q_in dx     inlet c = 1, c -> 0

The operator is IDENTICAL -- the reaction still sits on the diagonal as r dx^2
and the advection and diffusion are untouched. Only the right-hand side moves:
the volumetric source r dx^2 goes away, and the surface inflow, which
contributed nothing while it carried c = 0, becomes the only source. The rock
law loses its complement:

    now   dM/dt = -r (1 - c) / tau          rock stops where water is saturated
    then  dM/dt = -r c / tau_O2             rock stops where oxygen is gone

and the exponential integrator survives untouched, because lambda is still
independent of M.

WHAT IT FOUND: THE FLIP IS DECISIVELY BETTER, AND THE FIRST STATISTIC LIED

The first measurement taken here was the ratio of mean extent in joint cells
to mean extent in block interiors, at matched mean extent. It said the flip
was a disaster:

    silica  joints 0.8795  interiors 0.0159  ratio 55.3
    oxygen  joints 0.3839  interiors 0.1049  ratio  3.7

Two fixes were then tried against that number and neither worked. Gating k(M)
and tortuosity(M) on design 08's cracking criterion, so that nothing opens
below x_c = 0.10, moved it from 3.7 to 4.5. Raising tortuosity_fresh across
its whole measured range -- 1e4 to 3e4 to 1e5, all inside the 1e3-1e5 the
through-diffusion experiments give -- moved it 3.7 to 4.9 to 4.7,
non-monotone. Three levers, no effect: which was the signal that the number
was measuring the wrong thing.

THE PICTURE SHOWED IT AT ONCE. Under silica the joints are sharp dark lines
and the blocks are flat and untouched; under oxygen the blocks carry a smooth
graded rind with ROUNDED CORNERS. The second is spheroidal weathering. The
first is a grid.

The statistic that says so is extent against distance from the nearest joint:

    silica  0.00m 0.879   0.05m 0.087   0.10m 0.002   0.15m 0.001  ...
    oxygen  0.00m 0.384   0.05m 0.245   0.10m 0.153   0.15m 0.102
            0.20m 0.074   0.25m 0.056   0.30m 0.045   0.35m 0.037

Silica falls by a factor of four hundred in two cells: its rind is ONE CELL
WIDE, which is to say it is not resolved and its width is the grid's, not the
rock's. Oxygen decays with an e-folding near 10 cm and is still at 0.029 in
the block core -- a rind of 15 cm and more, resolved over six cells and
upward.

Measured rindlet zones in granite run 20 to 60 cm. Oxygen lands in that range;
silica is four to twelve times too thin. And the joint-to-interior ratio was
high for silica precisely BECAUSE the rind was unresolved -- every bit of
reaction sat in the joint cell itself, so the contrast being measured was
between a joint and everything else, not between a rind and a core.

This is the standing lesson in FRAME, met in the wild: every real defect in
this model has been caught by looking at a picture or reading an equation
against its code, never by an aggregate statistic.

Run:  PYTHONPATH=src python3 prototypes/probe_j_flip_the_solute.py
      PYTHONPATH=src python3 prototypes/probe_j_flip_the_solute.py --figure OUT.png
"""

import numpy as np
import scipy.sparse.linalg as spl

from corestone import FractureNetwork, Weathering, orthogonal_grid, YEAR
from corestone.weathering import ORDERING, _RTOL


class Oxidation(Weathering):
    """Weathering paced by oxidation of biotite Fe(II) by dissolved O2."""

    @property
    def reaction_coefficient(self):
        """``r = k_ox A`` [1/s], falling with the biotite remaining."""
        return self.specific_oxidation_coefficient * np.maximum(self.M, 0.0)

    def transport_coefficients(self):
        """
        The base method with D_aqueous -> D_O2_aqueous.

        Caught by looking at the numbers: the first run of this probe moved
        OXYGEN at the diffusivity of dissolved SILICA, because that is what
        the base class carries. O2 is about twice as mobile, and the
        penetration depth goes as its square root, so the error was 1.45x on
        the one length that decides whether a block interior is sheltered.
        """
        D = self.D_O2_aqueous
        tv, th = self._tort if self._tort is not None else self.link_tortuosity()
        dm_v = np.where(self.network.link_v, D, D / tv)
        dm_h = np.where(self.network.link_h, D, D / th)
        return (dm_v + self.dispersivity * np.abs(self.q_v) / self.dx,
                dm_h + self.dispersivity * np.abs(self.q_h) / self.dx)

    def _transport_operator(self):
        """Base method, but the seam term must use D_O2 too; the base builds
        it inline from ``self.D_aqueous``. Swapped for the call's duration."""
        saved = Weathering.D_aqueous
        try:
            Weathering.D_aqueous = property(lambda s: s.D_O2_aqueous)
            return Weathering._transport_operator(self)
        finally:
            Weathering.D_aqueous = saved

    def solve_solute(self, r):
        """
        The base method with ONE line changed: the right-hand side.

        Copied rather than refactored because this is a probe. If it survives,
        the base class grows a ``_solute_source`` and this override goes away.
        """
        dx = self.dx
        A = self._step_matrix(r)
        b = np.zeros(self.nz * self.nx)
        # Rain arrives air-saturated, so c = 1 at the surface by construction
        # -- there is no parameter here. An inflow with a known concentration
        # is a source, and it is the only one now.
        b.reshape(self.nz, self.nx)[0, :] = self.infiltration * dx

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
                A, b, x0=self._x, M=P, atol=0.0,
                callback=lambda xk: n_it.__setitem__(0, n_it[0] + 1),
                **{_RTOL: self.krylov_tol})
            if info != 0 or n_it[0] > self.max_krylov_iterations:
                x = direct()
        self._x = x
        return np.clip(x, 0.0, 1.0).reshape(self.nz, self.nx)

    def update(self, dt=None, dt_limit=None):
        """The base method, with ``(1 - c) / tau`` -> ``c / tau_O2``."""
        saved = (Weathering.specific_reaction_coefficient,
                 Weathering.tau)
        try:
            Weathering.specific_reaction_coefficient = property(
                lambda s: s.specific_oxidation_coefficient)
            Weathering.tau = property(lambda s: s.tau_oxidation)
            return self._update_with_flipped_rock_law(dt, dt_limit)
        finally:
            Weathering.specific_reaction_coefficient, Weathering.tau = saved

    def _update_with_flipped_rock_law(self, dt, dt_limit):
        """
        ``lambda = (r / M) c / tau_O2``.

        The base ``update`` forms ``(1 - c_held)``. Rather than copy sixty
        lines of step control to change one sign, the probe hands it a c that
        has already been complemented, and complements the stored field back
        afterwards so that ``self.c`` still means dissolved oxygen.
        """
        if self._c_held is None:
            self._c_held = self.solve_solute(self.reaction_coefficient)
        real_solve = self.solve_solute
        self.solve_solute = lambda r: 1.0 - real_solve(r)
        self._c_held = 1.0 - self._c_held
        try:
            step = Weathering.update(self, dt=dt, dt_limit=dt_limit)
        finally:
            self.solve_solute = real_solve
            if self._c_held is not None:
                self._c_held = 1.0 - self._c_held
            self.c = 1.0 - self.c
        return step


def build(cls, nx=60, nz=60, dx=0.05, spacing=1.0, tC=11.85):
    net = FractureNetwork(nz, nx, dx, periodic_x=True).seed(
        sets=orthogonal_grid(spacing), rng=np.random.default_rng(12345))
    m = cls(net)
    m.set_infiltration(0.30 / YEAR)
    m.set_temperature(tC + 273.15)
    m.initialize()
    return m


def profile(m, name):
    x = 1.0 - m.M
    joint = m.network.cell
    print("  %-10s mean %.4f  max %.4f  | joints %.4f  interiors %.4f"
          % (name, x.mean(), x.max(), x[joint].mean(), x[~joint].mean()))
    print("             c: surface %.4f  base %.4f  min %.4f  max %.4f"
          % (m.c[0, :].mean(), m.c[-1, :].mean(), m.c.min(), m.c.max()))
    rows = x.mean(axis=1)
    print("             by depth: " +
          "  ".join("%.0fm %.3f" % (i * m.dx, rows[i])
                    for i in range(0, m.nz, max(m.nz // 6, 1))))


def rind_profile(m, steps=8):
    """Extent of reaction against distance from the nearest joint [m].

    The measure that matches what the picture shows, and the one the
    joint-versus-interior ratio cannot express: a ratio is high when the rind
    is UNRESOLVED, because then all the reaction sits in the joint cell.
    """
    from scipy import ndimage
    x = 1.0 - m.M
    dist = ndimage.distance_transform_edt(~m.network.cell) * m.dx
    out = []
    for k in range(steps):
        sel = np.isclose(dist, k * m.dx, atol=1e-9)
        if sel.any():
            out.append((k * m.dx, float(x[sel].mean())))
    return out


def to_extent(cls, target, cap_years=3.0e7):
    """Advance until the mean extent of reaction reaches ``target``."""
    m = build(cls)
    t = 0.0
    while (1.0 - m.M).mean() < target and t < cap_years:
        t += m.update() / YEAR
    return m, t


def figure(path, target=0.10):
    """The side-by-side that settled it. Needs matplotlib; the probe does not."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    runs = [(name,) + to_extent(cls, target) for cls, name in
            ((Weathering, "silica (today)"), (Oxidation, "oxygen (design 08)"))]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.6), constrained_layout=True)
    m0 = runs[0][1]
    ext = [0, m0.nx * m0.dx, m0.nz * m0.dx, 0]
    for col, (name, m, t) in enumerate(runs):
        x = 1.0 - m.M
        j, i = x[m.network.cell].mean(), x[~m.network.cell].mean()
        a = axes[0, col]
        im = a.imshow(x, extent=ext, cmap="Oranges", vmin=0, vmax=1,
                      interpolation="nearest")
        a.set_title("%s\n%.0f kyr  -  mean %.2f  -  joints/interiors = %.1f"
                    % (name, t / 1e3, x.mean(), j / max(i, 1e-12)), fontsize=11)
        a.set_ylabel("Depth [m]")
        fig.colorbar(im, ax=a, label="extent of reaction")
        b = axes[1, col]
        im2 = b.imshow(m.c, extent=ext, cmap="Blues", vmin=0, vmax=1,
                       interpolation="nearest")
        b.set_xlabel("Distance [m]")
        b.set_ylabel("Depth [m]")
        b.set_title("solute: %s" % ("dissolved silica, c = C/C_eq" if col == 0
                                    else "dissolved O2, c = C/C_sat"),
                    fontsize=10)
        fig.colorbar(im2, ax=b, label="saturated ->" if col == 0
                     else "<- oxygen used up")
    fig.suptitle("The same rock at the same mean extent of reaction (%.0f %%), "
                 "driven two ways" % (100 * target), fontsize=13)
    fig.savefig(path, dpi=130)
    return path


if __name__ == "__main__":
    import sys
    if "--figure" in sys.argv:
        print("wrote", figure(sys.argv[sys.argv.index("--figure") + 1]))
        for cls, name in ((Weathering, "silica"), (Oxidation, "oxygen")):
            m, t = to_extent(cls, 0.10)
            print("  %-7s %.0f kyr, mean %.3f -- extent vs distance from joint:"
                  % (name, t / 1e3, (1.0 - m.M).mean()))
            print("     " + "  ".join("%.2fm:%.3f" % p for p in rind_profile(m)))
        raise SystemExit

    print("t = 0, fresh rock: where does the solute sit?\n")
    for cls, name in ((Weathering, "silica"), (Oxidation, "oxygen")):
        m = build(cls)
        m.c = m.solve_solute(m.reaction_coefficient)
        print("  %-8s c surface %.4f  base %.4f   (silica c rises with depth,"
              % (name, m.c[0, :].mean(), m.c[-1, :].mean()))
        print("           %sDamkohler %.4f)"
              % (" " * 8, m.damkohler if cls is Weathering
                 else m.oxidation_damkohler))
    print()

    for kyr in (100, 1000):
        print("after %d kyr:" % kyr)
        for cls, name in ((Weathering, "silica"), (Oxidation, "oxygen")):
            m = build(cls)
            m.run(years=kyr * 1e3)
            profile(m, name)
        print()

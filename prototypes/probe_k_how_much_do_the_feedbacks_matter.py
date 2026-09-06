"""
Probe K: how much do k(M) and tortuosity(M) actually do?

ANSWERED. They are the most important thing measured in this model, and they
matter far more to DISSOLUTION than to oxidation -- which decides where a
cracking gate may be applied and where it must not be.

SUPERSEDES probe_h_evolving_permeability.py, whose numbers predate three
corrections (a505892, c0d7749, 7cbd0a7) and are not to be trusted; it also
only tested the conductivity, and the tortuosity turns out to matter more.

METHOD. Hold k(M), tortuosity(M), or both, at their FRESH values by handing
the interpolation an M of all ones, and compare against the model as shipped
at matched model time. Three things are reported: the mean extent of
reaction, the rind -- extent against distance from the nearest joint, in
cells -- and the depth profile as the ratio of the top third of the section
to the bottom third.

    dissolution, 1000 kyr
      k(M)  tort(M)   mean     rind: 0, 1, 2, 3 cells out       top/bottom
      yes   yes      0.3072   0.9998  0.8537  0.3863  0.0900      1.485
      FIXED yes      0.2478   0.9997  0.7997  0.1089  0.0129      1.270
      yes   FIXED    0.2058   1.0000  0.4659  0.1959  0.0396      1.670
      FIXED FIXED    0.1039   1.0000  0.0403  0.0004  0.0002      1.029

    oxidation, 300 kyr
      yes   yes      0.3808   0.6303  0.5331  0.4365  0.3643      1.209
      FIXED yes      0.3057   0.6307  0.5295  0.4146  0.3042      1.137
      yes   FIXED    0.3291   0.6308  0.4567  0.3590  0.3061      1.264
      FIXED FIXED    0.1979   0.6318  0.3745  0.2230  0.1449      1.126

WHAT IT SAYS

1.  **Under dissolution the feedbacks ARE the rind.** Frozen, extent one cell
    from a joint falls from 0.85 to 0.04, and two cells out from 0.39 to
    0.0004 -- a factor of a thousand. The model becomes binary: joints
    entirely dissolved, blocks entirely untouched, nothing in between. The
    thin, grid-scale rind that made the dissolution driver look weak next to
    oxidation is what the model collapses TO when these are removed, not what
    it does.

2.  **They are also the depth profile.** Frozen, top-to-bottom is 1.03 --
    flat. Live, 1.49. Water opens the rock it passes through and arrives
    saturated below, and without that the section weathers uniformly with
    depth.

3.  **Tortuosity matters more than conductivity, and they are super-additive.**
    Freezing k alone costs 19 % of the mean; freezing tortuosity alone costs
    33 %; freezing both costs 66 %. Opening the rock to flow and opening it to
    diffusion reinforce each other, which is the reactive-infiltration
    feedback doing what it does everywhere else.

    (Freezing tortuosity RAISES the depth ratio, to 1.67. Shallow rock
    weathers first, and without the matrix opening behind it that head start
    is never given back.)

4.  **Oxidation depends on them far less.** Both frozen costs 48 % of the mean
    against dissolution's 66 %, and the joint cell does not move at all
    (0.6303 to 0.6318). Oxygen already penetrates 4.5 cm into INTACT rock by
    diffusion, so it does not need the rock opened for it; silica does.

THE CONSEQUENCE FOR DESIGN 10

Gating k(M) and tortuosity(M) on the cracking criterion is the obvious next
step and it must NOT be applied to the dissolution driver. The physics is the
reason, and it is clean:

  * Dissolving, MASS LEAVES THE ROCK. Porosity opens directly, with no
    fracture required, and gating it behind a cracking threshold would
    suppress the mechanism that makes the rind -- turning the in-class model
    back into the binary picture in row four above.
  * Oxidising, NOTHING LEAVES. The iron is oxidised in place, and Goodfellow
    et al.'s own line is "major changes in rock properties can occur with only
    minor element leaching". The rock opens only because it CRACKS, so the
    gate belongs here and nowhere else.

That is also why Goodfellow's conductivity series is measured against Fe(III)
rather than against leaching: they were measuring the oxidation pathway.

Run:  PYTHONPATH=src python3 prototypes/probe_k_how_much_do_the_feedbacks_matter.py
"""

import numpy as np
from scipy import ndimage

from corestone import FractureNetwork, Weathering, orthogonal_grid, YEAR


class Frozen(Weathering):
    """Weathering with k(M) and/or tortuosity(M) held at their fresh values."""

    freeze_k = False
    freeze_tort = False

    def _fresh(self, method):
        saved, self.M = self.M, np.ones((self.nz, self.nx))
        try:
            return method(self)
        finally:
            self.M = saved

    def link_conductivity(self):
        if not self.freeze_k:
            return Weathering.link_conductivity(self)
        return self._fresh(Weathering.link_conductivity)

    def link_tortuosity(self):
        if not self.freeze_tort:
            return Weathering.link_tortuosity(self)
        return self._fresh(Weathering.link_tortuosity)

    def link_tortuosity_wrap(self):
        if not self.freeze_tort:
            return Weathering.link_tortuosity_wrap(self)
        return self._fresh(Weathering.link_tortuosity_wrap)


def run(driver, years, freeze_k=False, freeze_tort=False):
    net = FractureNetwork(60, 60, 0.05, periodic_x=True).seed(
        sets=orthogonal_grid(1.0), rng=np.random.default_rng(12345))
    m = Frozen(net)
    m.freeze_k, m.freeze_tort = freeze_k, freeze_tort
    m.set_driver(driver)
    m.set_infiltration(0.30 / YEAR)
    m.set_temperature(285.0)
    m.initialize()
    m.run(years=years)
    return m


def describe(m, cells=4):
    """Mean extent, the rind cell by cell, and the depth ratio."""
    x = 1.0 - m.M
    dist = ndimage.distance_transform_edt(~m.network.cell) * m.dx
    rind = [float(x[np.isclose(dist, k * m.dx)].mean()) for k in range(cells)]
    rows = x.mean(axis=1)
    third = max(m.nz // 3, 1)
    return (float(x.mean()), rind,
            float(rows[:third].mean() / max(rows[-third:].mean(), 1e-12)))


if __name__ == "__main__":
    for driver, years in (("dissolution", 1000e3), ("oxidation", 300e3)):
        print("=== %s, %.0f kyr ===" % (driver, years / 1e3))
        print("  k(M)  tort(M)   mean     rind: 0, 1, 2, 3 cells out"
              "       top/bottom")
        for fk, ft, lab in ((False, False, "yes   yes  "),
                            (True, False, "FIXED yes  "),
                            (False, True, "yes   FIXED"),
                            (True, True, "FIXED FIXED")):
            mean, rind, ratio = describe(run(driver, years, fk, ft))
            print("  %s  %.4f   %s      %.3f"
                  % (lab, mean, "  ".join("%.4f" % r for r in rind), ratio))
        print()

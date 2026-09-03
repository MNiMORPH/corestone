"""
Probe H: should the flow field evolve as the rock weathers?

The model solves the head ONCE, from a conductance that depends only on the
joint network. Dissolving rock opens porosity, so weathered rock should
conduct better than fresh rock, which draws more water into the weathered
zone, which dissolves it faster still. That is a positive feedback -- the one
that makes wormholes in limestone and reactive-infiltration fingers in every
dissolving porous medium -- and leaving it out is a simplification nobody had
weighed.

The question is not whether the feedback is real. It is whether it changes
what this demo shows, and by how much, because the answer decides between
"document the assumption" and "rewrite the flow step".

The closure. Conductivity is interpolated geometrically -- linearly in its
logarithm, which is how conductivity varies -- between intact granite at
M = 1 and grus at M = 0:

    log k = M log k_matrix + (1 - M) log k_grus

on the mean of the two cells a link joins. A jointed link keeps k_fracture:
an open joint is an open joint whatever the rock beside it has done.

k_grus is a NEW PARAMETER and the only one here. Set equal to k_fracture, so
that fully dissolved rock conducts as well as a joint and no new order of
magnitude is invented; the joint network in effect grows into the weathered
zone. Every other parameter in this model is a placeholder and so is this.

Run:  PYTHONPATH=src python3 prototypes/probe_h_evolving_permeability.py
"""

import time

import numpy as np

from corestone import FractureNetwork, Weathering, orthogonal_grid, YEAR


class EvolvingK(Weathering):
    """Weathering with a conductivity that follows the rock."""

    k_grus = 1.0e-5                 # = k_fracture. PLACEHOLDER, see above.

    def link_conductivity(self):
        net = self.network
        lo, hi = np.log(self.k_matrix), np.log(self.k_grus)

        def k_of(m):
            """Geometric interpolation between intact granite and grus."""
            return np.exp(m * lo + (1.0 - m) * hi)

        M = np.clip(self.M, 0.0, 1.0)
        kv = np.where(net.link_v, self.k_fracture,
                      k_of(0.5 * (M[:-1, :] + M[1:, :])))
        kh = np.where(net.link_h, self.k_fracture,
                      k_of(0.5 * (M[:, :-1] + M[:, 1:])))
        if net.periodic_x:
            kw = np.where(net.link_wrap, self.k_fracture,
                          k_of(0.5 * (M[:, -1] + M[:, 0])))
        else:
            kw = np.zeros(self.nz)
        return kv, kh, kw


def build(cls, nx=60, dx=0.05, spacing=1.0):
    net = FractureNetwork(nx, nx, dx, periodic_x=True).seed(
        sets=orthogonal_grid(spacing), rng=np.random.default_rng(12345))
    m = cls(net)
    m.set_infiltration(0.30 / YEAR)
    m.c_drift_max = 0.01
    m.initialize()
    return m


def run(cls, kyr, resolve_every=None):
    """Advance to ``kyr``, re-solving the flow every ``resolve_every`` steps."""
    m = build(cls)
    target = kyr * 1e3 * YEAR
    n = flows = 0
    t0 = time.perf_counter()
    while m.t < target - 1e-9 * YEAR:
        m.update(dt_limit=target - m.t)
        n += 1
        if resolve_every and n % resolve_every == 0:
            m.solve_flow()
            flows += 1
    return time.perf_counter() - t0, n, flows, m


def report(name, elapsed, n, flows, m, ref=None):
    d = m.dissolved_fraction
    line = ("%-26s %6.2fs %4d steps %4d flows | grus %5.1f %% | "
            "mean %.4f  max %.4f" % (name, elapsed, n, flows,
                                     100 * m.is_grus.mean(), d.mean(), d.max()))
    if ref is not None:
        line += "  | max|diff| %.4f" % np.abs(d - ref).max()
    print(line)
    return d


if __name__ == "__main__":
    KYR = 100

    print("Does re-solving the flow change the answer?  %d kyr, 3 m at 5 cm\n"
          % KYR)

    t, n, f, m = run(Weathering, KYR)
    static = report("static k (as shipped)", t, n, f, m)

    for every in (10, 5, 1):
        t, n, f, m = run(EvolvingK, KYR, resolve_every=every)
        report("k(M), re-solved every %2d" % every, t, n, f, m, static)

    # Where does it differ? The feedback should focus flow into the weathered
    # zone, so the joints should run further ahead of the block interiors.
    t, n, f, m = run(EvolvingK, KYR, resolve_every=1)
    net = m.network
    joint = net.cell
    print("\nblock interiors vs joint cells, dissolved fraction:")
    for label, d in (("static", static), ("evolving", m.dissolved_fraction)):
        print("  %-9s joints %.4f   interior %.4f   ratio %.3f"
              % (label, d[joint].mean(), d[~joint].mean(),
                 d[joint].mean() / max(d[~joint].mean(), 1e-12)))

    # And the cost of the flow solve itself, which is what decides how often
    # it can be afforded.
    m = build(Weathering)
    t0 = time.perf_counter()
    for _ in range(10):
        m.solve_flow()
    print("\none flow solve: %.1f ms   (a weathering step is ~2 ms here)"
          % (1e2 * (time.perf_counter() - t0)))

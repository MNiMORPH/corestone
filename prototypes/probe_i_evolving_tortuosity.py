"""
Probe I: should the matrix tortuosity evolve with M, as the conductivity does?

THE ASYMMETRY. link_conductivity interpolates k geometrically between intact
granite and fully weathered rock as the soluble phase goes. The tortuosity
does not: it is fixed at 10, which is a WEATHERED value -- saprolite at ~30 %
porosity gives D_eff/D_0 ~ 0.1. Intact granite is 2e-14 to 1.3e-12 m2/s
against a free-water 1e-9, so a tortuosity of 1e3 to 1e5, centre 1e4.

So fresh rock in this model diffuses about a thousand times too freely, and it
does so exactly when every cell is fresh -- at t = 0, when the rind is forming.
That flatters the feature the model exists to show.

Dissolving rock opens porosity to diffusion as surely as it opens it to flow,
so the fix is the same interpolation:

    tortuosity(M) = tortuosity_fresh^M * tortuosity_weathered^(1 - M)

WHAT THIS PROBE IS FOR. Two things could make the change not worth making, and
both are measurable before touching src/:

  1. Dispersion may simply take over. The dispersive term in matrix cells is
     ~4.75e-13 m2/s. If molecular drops to 1e-13 then dispersion wins ~5:1
     where today molecular wins 200:1 -- and the dispersivity is a bulk
     aquifer number (Gelhar et al. 1992) with no business describing intact
     granite at 1e-11 m/s. That would move the error rather than remove it.

  2. The answer may not be constrained. The literature spans 1e3 to 1e5. If
     t90 swings wildly across that range we would be trading a known cheat for
     an unconstrained parameter, which is a worse trade.

And a prediction to be judged against: the model already runs at 2.5 m/Myr
against 4-7 measured for temperate granite regoliths, so slowing early
weathering should make the validation WORSE, not better.

Run:  PYTHONPATH=src python3 prototypes/probe_i_evolving_tortuosity.py
"""

import time

import numpy as np

from corestone import FractureNetwork, Weathering, orthogonal_grid, YEAR

TORT_WEATHERED = 10.0


class EvolvingTortuosity(Weathering):
    """Tortuosity interpolated between fresh and weathered, like k(M)."""

    tortuosity_fresh = 1.0e4

    def transport_coefficients(self):
        D = self.D_aqueous
        M = np.clip(self.M, 0.0, 1.0)
        lo, hi = np.log(self.tortuosity_fresh), np.log(TORT_WEATHERED)

        def tort(m):
            return np.exp(m * lo + (1.0 - m) * hi)

        tv = tort(0.5 * (M[:-1, :] + M[1:, :]))
        th = tort(0.5 * (M[:, :-1] + M[:, 1:]))
        dm_v = np.where(self.network.link_v, D, D / tv)
        dm_h = np.where(self.network.link_h, D, D / th)
        return (dm_v + self.dispersivity * np.abs(self.q_v) / self.dx,
                dm_h + self.dispersivity * np.abs(self.q_h) / self.dx)


def build(cls, tC=12.0, **kw):
    net = FractureNetwork(60, 60, 0.05, periodic_x=True).seed(
        sets=orthogonal_grid(1.0), rng=np.random.default_rng(12345))
    m = cls(net)
    for k, v in kw.items():
        setattr(m, k, v)
    m.set_infiltration(0.30 / YEAR)
    m.set_temperature(tC + 273.15)
    m.c_drift_max = 0.01
    m.flow_tolerance = 0.05
    m.dt_max = 1000.0 * YEAR
    return m.initialize()


def to_90(m, cap=40000.0):
    t0 = time.perf_counter()
    n = 0
    while float(m.dissolved_fraction.mean()) < 0.90 and m.t < cap * 1e3 * YEAR:
        m.update()
        n += 1
    return m.t / (1e3 * YEAR), n, time.perf_counter() - t0


def rind(m):
    """Sharpness: how much of the section sits between 20 % and 80 % gone."""
    d = m.dissolved_fraction
    return float(((d > 0.2) & (d < 0.8)).mean())


if __name__ == "__main__":
    print("1. WHICH TRANSPORT TERM RULES THE MATRIX, fresh vs weathered\n")
    m = build(Weathering)
    v = m.darcy_speed
    matrix = ~m.network.cell
    disp = m.dispersivity * float(np.median(v[matrix]))
    print("   dispersive term in matrix cells      %.3e m2/s" % disp)
    for name, tort in (("as shipped (10, weathered)", 10.0),
                       ("fresh granite (1e4)", 1.0e4)):
        mol = m.D_aqueous / tort
        print("   molecular at tortuosity %-9s %.3e   molecular/dispersive %8.3f"
              % (name.split("(")[1].rstrip(")"), mol, mol / disp))

    print("\n2. HOW MUCH IS THE ANSWER CONSTRAINED? sweep the fresh value\n")
    base_kyr, base_n, base_s = to_90(build(Weathering))
    print("   %-34s t90 %7.0f kyr   %4d steps  %5.1f s   front %.2f m/Myr"
          % ("as shipped (fixed 10)", base_kyr, base_n, base_s, 3.0 / (base_kyr / 1e3)))
    for tf in (1.0e3, 1.0e4, 1.0e5):
        mm = build(EvolvingTortuosity, tortuosity_fresh=tf)
        kyr, n, s = to_90(mm)
        print("   %-34s t90 %7.0f kyr   %4d steps  %5.1f s   front %.2f m/Myr"
              % ("evolving, fresh = %.0e" % tf, kyr, n, s, 3.0 / (kyr / 1e3)))

    print("\n3. DOES THE RIND SHARPEN? fraction of the section part-dissolved")
    print("   (measured at 50 % overall dissolution, so like against like)\n")
    for label, cls, kw in (("as shipped", Weathering, {}),
                           ("evolving, fresh = 1e4", EvolvingTortuosity,
                            dict(tortuosity_fresh=1.0e4))):
        mm = build(cls, **kw)
        while float(mm.dissolved_fraction.mean()) < 0.50 and mm.t < 40000e3 * YEAR:
            mm.update()
        print("   %-24s %5.1f %% of cells between 20 %% and 80 %% gone"
              % (label, 100 * rind(mm)))

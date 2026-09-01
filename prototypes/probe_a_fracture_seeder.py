"""
Probe A: how far is granite from the nearest fracture?

FIRST VERSION ASKED THE WRONG QUESTION. It measured connected components of
cells with fractured links removed -- i.e. it treated fractures as barriers
partitioning the rock into blocks. They are not. A joint is a conduit; the rock
is mechanically and hydraulically continuous across it, so removing fractured
links can never disconnect anything, and every test "leaked" including a single
fracture spanning the whole domain.

A corestone is not rock fenced off by fractures. It is rock the water never
reaches: rock whose DISTANCE from the fracture network exceeds how far
weathering can penetrate from a fracture wall in the time available. So the
quantity that matters is the distance transform from the network, and the
question this probe answers is whether the seeded network leaves distances of
the order of the joint half-spacing -- or leaves most of the domain metres from
any water.

Ugly on purpose. Run:  python3 prototypes/probe_a_fracture_seeder.py
"""
import numpy as np
from scipy.ndimage import distance_transform_edt

# ---- domain (see design/01-fracture-seeding.md) -------------------------------
LX, LZ, DX = 20.0, 15.0, 0.05
NX, NZ = int(round(LX / DX)), int(round(LZ / DX))

# ---- joint sets: ALL OF THESE ARE PROPOSED VALUES, not measured ---------------
SETS = [
    dict(name="J1", dip_deg=+75.0, kappa=20.0, spacing=1.5),
    dict(name="J2", dip_deg=-75.0, kappa=20.0, spacing=1.5),
    # Sub-horizontal (sheeting) set. Two near-vertical conjugate sets cannot
    # bound a block in the vertical, no matter how persistent -- see the
    # numbers in design/01-fracture-seeding.md. This set is what closes it.
    dict(name="SH", dip_deg=  5.0, kappa=40.0, spacing=1.5),
]
LEN_MIN, LEN_EXP = 4.0, 2.0          # trace length ~ power law, p(L) ~ L^-LEN_EXP
PAD = 3.0                            # seed this far outside the domain, then clip
RNG = np.random.default_rng(12345)


def sample_lengths(n):
    """Bounded power law, exponent LEN_EXP, min LEN_MIN, capped at the diagonal."""
    u = RNG.random(n)
    L = LEN_MIN * (1.0 - u) ** (-1.0 / (LEN_EXP - 1.0))
    return np.minimum(L, np.hypot(LX, LZ))


def segments_for_set(spec):
    """Center-lines offset along the set normal at drawn spacings."""
    theta = np.deg2rad(spec["dip_deg"])
    d = np.array([np.cos(theta), np.sin(theta)])      # along-fracture (x, depth)
    n = np.array([-d[1], d[0]])                       # normal

    # Seed over a padded box and let rasterize() clip. Without the pad no
    # fracture can be centred just outside the domain, so edge cells come out
    # artificially far from the network -- an edge artifact, not geology.
    corners = np.array([[-PAD, -PAD], [LX + PAD, -PAD],
                        [-PAD, LZ + PAD], [LX + PAD, LZ + PAD]], dtype=float)
    s_c, t_c = corners @ d, corners @ n
    segs, t = [], t_c.min()
    while t < t_c.max():
        # perturb this fracture's orientation about the set mean (von Mises)
        th = np.deg2rad(spec["dip_deg"]) + RNG.vonmises(0.0, spec["kappa"])
        dd = np.array([np.cos(th), np.sin(th)])
        L = sample_lengths(1)[0]
        s_mid = RNG.uniform(s_c.min(), s_c.max())
        mid = s_mid * d + t * n
        segs.append((mid - 0.5 * L * dd, mid + 0.5 * L * dd))
        t += RNG.lognormal(np.log(spec["spacing"]), 0.35)
    return segs


def rasterize(segs):
    """Mark the cells a fracture passes through, and the links along its trace.

    Links along the trace are the ones that carry fracture flow: a vertically
    running joint enhances the VERTICAL links between the cells it threads.
    Cells are what the distance transform needs.
    """
    v = np.zeros((NZ - 1, NX), dtype=bool)     # link between [iz,ix] and [iz+1,ix]
    h = np.zeros((NZ, NX - 1), dtype=bool)     # link between [iz,ix] and [iz,ix+1]
    cell = np.zeros((NZ, NX), dtype=bool)      # cells the trace passes through
    trace_length = 0.0
    for p0, p1 in segs:
        L = np.hypot(*(p1 - p0))
        ns = max(int(np.ceil(L / (0.25 * DX))), 2)
        pts = p0 + np.linspace(0, 1, ns)[:, None] * (p1 - p0)
        ix = np.floor(pts[:, 0] / DX).astype(int)
        iz = np.floor(pts[:, 1] / DX).astype(int)
        ok = (ix >= 0) & (ix < NX) & (iz >= 0) & (iz < NZ)
        ix, iz = ix[ok], iz[ok]
        if ix.size < 2:
            continue
        trace_length += L * ok.mean()
        cell[iz, ix] = True
        dix, diz = np.diff(ix), np.diff(iz)
        for k in np.nonzero((dix != 0) | (diz != 0))[0]:
            if diz[k] != 0:                                   # crossed vertically
                v[min(iz[k], iz[k + 1]), ix[k]] = True
            if dix[k] != 0:                                   # crossed horizontally
                h[iz[k], min(ix[k], ix[k + 1])] = True
    return v, h, cell, trace_length


segs = []
for spec in SETS:
    s = segments_for_set(spec)
    print(f"{spec['name']}: {len(s):4d} fractures, mean dip {spec['dip_deg']:+.0f} deg, "
          f"mean spacing {spec['spacing']} m")
    segs += s

v, h, cell, tl = rasterize(segs)
area = LX * LZ

# Distance from every cell to the nearest fracture [m].
dist = distance_transform_edt(~cell, sampling=DX)

print()
print(f"grid                {NZ} x {NX} = {NZ*NX} cells at dx = {DX} m "
      f"({LZ} x {LX} m)")
print(f"fractures           {len(segs)}")
print(f"P21 intensity       {tl/area:.3f} m/m2   (trace length per area)")
print(f"links fractured     vertical {v.mean()*100:.2f} %, "
      f"horizontal {h.mean()*100:.2f} %")
print(f"distance to nearest fracture [m]:")
print(f"  median {np.median(dist):.2f}   p90 {np.percentile(dist,90):.2f}   "
      f"max {dist.max():.2f}")
for r in (0.05, 0.10, 0.25, 0.50, 1.00):
    print(f"  within {r:4.2f} m of a fracture: {(dist <= r).mean()*100:5.1f} % of rock")

# ---- reference: the SAME measurement on a perfectly regular network ----------
# Comparing against an analytic guess invites comparing two different things.
# Build the regular case the same way and measure it the same way.
def regular_reference(spacing):
    segs, diag = [], 2 * np.hypot(LX, LZ)
    for spec in SETS:
        th = np.deg2rad(spec["dip_deg"])
        d = np.array([np.cos(th), np.sin(th)])
        n = np.array([-d[1], d[0]])
        corners = np.array([[-PAD, -PAD], [LX + PAD, -PAD],
                            [-PAD, LZ + PAD], [LX + PAD, LZ + PAD]], dtype=float)
        t_c = corners @ n
        s_c = corners @ d
        for t in np.arange(t_c.min(), t_c.max(), spacing):
            mid = 0.5 * (s_c.min() + s_c.max()) * d + t * n
            segs.append((mid - 0.5 * diag * d, mid + 0.5 * diag * d))
    _, _, c, _ = rasterize(segs)
    return distance_transform_edt(~c, sampling=DX)

ref = regular_reference(SETS[0]["spacing"])
print()
print(f"REFERENCE -- the same three sets, perfectly regular, fully persistent, "
      f"spacing {SETS[0]['spacing']} m:")
print(f"  median {np.median(ref):.2f}   p90 {np.percentile(ref,90):.2f}   "
      f"max {ref.max():.2f}")
print()
print("The seeded network should match the reference in the median. A heavier "
      "tail is real:")
print("it is the lognormal spacing and finite trace length, and it is what "
      "makes large corestones.")

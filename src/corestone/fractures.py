#! /usr/bin/python3

"""
Seed a joint network into a 2D vertical cross section.

Granite weathers where water reaches it, and water reaches it along joints, so
the joint network decides -- before any chemistry exists -- which rock can
weather and which cannot. This module builds that network.

Joints are placed by *set*: a family of subparallel fractures with a
characteristic orientation and spacing. Real granite carries at least three:
two steep conjugate sets and a subhorizontal sheeting set. Two steep sets alone
cannot bound a block in the vertical at any persistence -- see
``design/01-fracture-seeding.md`` for the measurement.

What a fracture *is*, in this discretization: a conduit, not a barrier. Rock is
continuous across a joint. A fracture marks the cells it threads and raises the
conductance of the links along its trace; it never disconnects anything.
"""

import numpy as np
from scipy.ndimage import distance_transform_edt


class JointSet(object):
    """
    One family of subparallel joints.

    dip_deg  : mean orientation, degrees from horizontal, in the (x, depth)
               plane. Positive dips to the right, so +90 is vertical.
    kappa    : von Mises concentration of the orientation scatter. Larger is
               tighter; 20 is roughly +/- 13 degrees.
    spacing  : mean normal spacing between joints of this set [m].
    """

    def __init__(self, name, dip_deg, kappa, spacing,
                 spacing_sigma=0.35, length_min=4.0, length_exponent=2.0):
        self.name = name
        self.dip_deg = dip_deg            # mean dip [deg from horizontal]
        self.kappa = kappa                # von Mises concentration [-]
        self.spacing = spacing            # mean normal spacing [m]
        self.spacing_sigma = spacing_sigma  # lognormal sigma of the spacing
        self.length_min = length_min      # minimum trace length [m]
        self.length_exponent = length_exponent  # p(L) ~ L**-exponent


#: The default granite case: two steep conjugate sets plus a sheeting set.
#: PLACEHOLDER VALUES -- replace with field joint spacings and orientations.
GRANITE_SETS = [
    JointSet("J1", dip_deg=+75.0, kappa=20.0, spacing=1.5),
    JointSet("J2", dip_deg=-75.0, kappa=20.0, spacing=1.5),
    JointSet("SH", dip_deg=+5.0, kappa=40.0, spacing=1.5),
]


class FractureNetwork(object):
    """
    A seeded joint network on a cell/link grid.

    Cells are `(nz, nx)`; `iz` increases downward from the ground surface.
    Fracture state lives on the links between cell centres: `(nz-1, nx)`
    vertical and `(nz, nx-1)` horizontal.
    """

    def __init__(self, nz, nx, dx):
        self.nz = nz                      # rows, increasing downward
        self.nx = nx                      # columns, increasing rightward
        self.dx = dx                      # cell size [m], square cells

        self.cell = None                  # (nz, nx) bool: fracture threads it
        self.link_v = None                # (nz-1, nx) bool: vertical links
        self.link_h = None                # (nz, nx-1) bool: horizontal links
        self.segments = None              # list of (p0, p1) endpoints [m]
        self.trace_length = None          # total trace length in domain [m]

    # ---- geometry

    @property
    def lx(self):
        """Domain width [m]."""
        return self.nx * self.dx

    @property
    def lz(self):
        """Domain depth [m]."""
        return self.nz * self.dx

    @property
    def p21(self):
        """Fracture intensity: trace length per unit area [m/m2]."""
        return self.trace_length / (self.lx * self.lz)

    def distance_to_fracture(self):
        """
        Distance from every cell to the nearest fracture [m].

        The quantity that decides corestones: rock further from a fracture than
        the weathering can reach never alters, whatever its mineralogy.
        """
        return distance_transform_edt(~self.cell, sampling=self.dx)

    # ---- seeding

    def seed(self, sets=None, pad=3.0, rng=None):
        """
        Place the joints and mark the cells and links they occupy.

        `pad` seeds fractures this far outside the domain before clipping.
        Without it no fracture can be centred just beyond the edge, and edge
        cells come out artificially far from the network -- an artifact of the
        boundary, not geology.
        """
        sets = GRANITE_SETS if sets is None else sets
        rng = np.random.default_rng() if rng is None else rng

        self.segments = []
        for js in sets:
            self.segments += self._segments_for_set(js, pad, rng)
        self._rasterize()
        return self

    def _segments_for_set(self, js, pad, rng):
        """Centre-lines offset along the set normal at drawn spacings."""
        theta = np.deg2rad(js.dip_deg)
        d = np.array([np.cos(theta), np.sin(theta)])     # along-fracture
        n = np.array([-d[1], d[0]])                      # set normal

        box = np.array([[-pad, -pad], [self.lx + pad, -pad],
                        [-pad, self.lz + pad], [self.lx + pad, self.lz + pad]])
        s_c, t_c = box @ d, box @ n

        diag = np.hypot(self.lx, self.lz)
        segs, t = [], t_c.min()
        while t < t_c.max():
            th = theta + rng.vonmises(0.0, js.kappa)
            dd = np.array([np.cos(th), np.sin(th)])
            u = rng.random()
            L = min(js.length_min * (1.0 - u) ** (-1.0 / (js.length_exponent - 1.0)),
                    diag)
            mid = rng.uniform(s_c.min(), s_c.max()) * d + t * n
            segs.append((mid - 0.5 * L * dd, mid + 0.5 * L * dd))
            t += rng.lognormal(np.log(js.spacing), js.spacing_sigma)
        return segs

    def _rasterize(self):
        """Mark the cells each trace threads, and the links along it."""
        self.cell = np.zeros((self.nz, self.nx), dtype=bool)
        self.link_v = np.zeros((self.nz - 1, self.nx), dtype=bool)
        self.link_h = np.zeros((self.nz, self.nx - 1), dtype=bool)
        self.trace_length = 0.0

        for p0, p1 in self.segments:
            L = np.hypot(*(p1 - p0))
            ns = max(int(np.ceil(L / (0.25 * self.dx))), 2)
            pts = p0 + np.linspace(0.0, 1.0, ns)[:, None] * (p1 - p0)
            ix = np.floor(pts[:, 0] / self.dx).astype(int)
            iz = np.floor(pts[:, 1] / self.dx).astype(int)
            ok = (ix >= 0) & (ix < self.nx) & (iz >= 0) & (iz < self.nz)
            if ok.sum() < 2:
                continue
            self.trace_length += L * ok.mean()
            ix, iz = ix[ok], iz[ok]
            self.cell[iz, ix] = True

            dix, diz = np.diff(ix), np.diff(iz)
            for k in np.nonzero((dix != 0) | (diz != 0))[0]:
                if diz[k] != 0:
                    self.link_v[min(iz[k], iz[k + 1]), ix[k]] = True
                if dix[k] != 0:
                    self.link_h[iz[k], min(ix[k], ix[k + 1])] = True

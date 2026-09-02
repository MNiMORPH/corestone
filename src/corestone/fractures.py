#! /usr/bin/python3

"""
Seed a joint network into a 2D vertical cross section.

Granite weathers where water reaches it, and water reaches it along joints, so
the joint network decides -- before any chemistry exists -- which rock can
weather and which cannot. This module builds that network.

The generator follows the standard outcrop hierarchy rather than scattering
free-floating segments:

  - a **primary** set of *throughgoing* joints, each spanning the domain;
  - a **secondary** set that **abuts** the primary one, each trace running from
    one primary joint to the next and terminating there.

That abutting rule is what real joint sets do, and it is what makes a network
connected: it produces Y (T-shaped) nodes rather than I nodes (isolated tips).
An earlier generator here placed both sets as free segments drawn from a
power-law length distribution, and the network never linked up -- see
``design/01-fracture-seeding.md`` and ``design/03-throughgoing-joints.md``.

What a fracture *is*, in this discretization: a conduit, not a barrier. Rock is
continuous across a joint. A fracture marks the cells it threads and raises the
conductance of the links along its trace; it never disconnects anything.

Provenance of the method: the throughgoing-plus-abutting construction and the
X/Y/I topological description are standard in discrete fracture network work
(Sanderson & Nixon 2015 for the topology; ADFNE, Fadakar-Alghalandis 2017, and
FracSim2D for the generators). No existing Python generator was usable here --
they are MATLAB, C++, Python 2, or 3D-only -- so the algorithm is reimplemented
and the network is validated against ``fractopo``.
"""

import numpy as np
from scipy.ndimage import distance_transform_edt


class JointSet(object):
    """
    One family of subparallel joints.

    dip_deg  : mean orientation, degrees from horizontal, in the (x, depth)
               plane. 0 is horizontal, 90 is vertical.
    spacing  : mean normal spacing between joints of this set [m].
    kappa    : von Mises concentration of the orientation scatter. Larger is
               tighter; 25 is roughly +/- 11 degrees. ``None`` means no
               scatter at all -- every joint of the set exactly parallel.
    abuts    : name of the set these joints terminate against, or ``None`` for
               a throughgoing set that spans the domain.
    spans    : how many primary joints an abutting trace crosses before it
               stops. 1 means it terminates at the very next one.
    spacing_sigma : lognormal sigma of the spacing. ``0`` means exact,
               evenly spaced joints.
    density  : fraction of the available gaps between primary joints that
               carry a trace of this set. 1.0 fills every one, giving a clean
               orthogonal block system; lower values leave gaps and make the
               blocks more irregular.
    """

    def __init__(self, name, dip_deg, spacing, kappa=25.0,
                 spacing_sigma=0.30, abuts=None, spans=1, density=1.0):
        self.name = name
        self.dip_deg = dip_deg
        self.spacing = spacing
        self.kappa = kappa
        self.spacing_sigma = spacing_sigma
        self.abuts = abuts
        self.spans = spans
        self.density = density

    @property
    def is_throughgoing(self):
        return self.abuts is None


def conjugate_sets(dip_primary=90.0, dip_secondary=0.0, spacing=1.5,
                   kappa=25.0, spacing_sigma=0.30, spans=1, density=1.0):
    """
    A conjugate pair: a throughgoing set and a set that abuts it.

    The default is the orthogonal case -- vertical joints cut by horizontal
    ones, 90 degrees apart -- which is the simplest geometry that bounds a
    block on all four sides. For a symmetric shear pair instead, pass
    ``dip_primary=+45, dip_secondary=-45``: still 90 degrees apart, differently
    oriented relative to the surface.
    """
    return [
        JointSet("J1", dip_deg=dip_primary, spacing=spacing, kappa=kappa,
                 spacing_sigma=spacing_sigma),
        JointSet("J2", dip_deg=dip_secondary, spacing=spacing, kappa=kappa,
                 spacing_sigma=spacing_sigma, abuts="J1", spans=spans,
                 density=density),
    ]


def uniform_grid_shape(width, depth, dx, spacing):
    """
    Cell counts that let a regular joint set tile the domain exactly.

    A joint sits on each wall and every block between them is a full spacing
    across, which needs ``n * spacing / dx + 1`` cells along each axis. The odd
    ``+ 1`` also makes the count odd whenever the spacing is an even number of
    cells, which is what lets the pattern rasterise mirror-symmetrically -- see
    :meth:`FractureNetwork._rasterize`.

    Returns ``(nz, nx)`` for the largest such domain no larger than the size
    asked for.
    """
    per = spacing / dx
    nx = int(np.floor(width / dx / per)) * int(round(per)) + 1
    nz = int(np.floor(depth / dx / per)) * int(round(per)) + 1
    return nz, nx


def orthogonal_grid(spacing=1.5, density=1.0):
    """
    A perfectly regular vertical/horizontal joint grid.

    No orientation scatter and no spacing variability: every joint exactly
    vertical or exactly horizontal, evenly spaced, so the blocks are identical
    squares of side ``spacing``. The clearest object to teach from, and the
    case whose fracture intensity has an exact answer, P21 = 2 / spacing.
    """
    return conjugate_sets(dip_primary=90.0, dip_secondary=0.0, spacing=spacing,
                          kappa=None, spacing_sigma=0.0, density=density)


#: The default granite case. PLACEHOLDER SPACING -- replace with field values.
GRANITE_SETS = orthogonal_grid()


def _unit(dip_deg):
    """Unit vector along a joint of this dip, in (x, depth)."""
    th = np.deg2rad(dip_deg)
    return np.array([np.cos(th), np.sin(th)])


def _clip_to_box(origin, direction, lx, lz):
    """
    Clip an infinite line to the domain rectangle. Returns (p0, p1) or None.

    Liang-Barsky on the parametric line ``origin + t * direction``.
    """
    t0, t1 = -np.inf, np.inf
    for p, q in ((-direction[0], origin[0]), (direction[0], lx - origin[0]),
                 (-direction[1], origin[1]), (direction[1], lz - origin[1])):
        if p == 0.0:
            if q < 0.0:
                return None                      # parallel and outside
            continue
        t = q / p
        if p < 0.0:
            t0 = max(t0, t)
        else:
            t1 = min(t1, t)
    if t0 >= t1:
        return None
    return origin + t0 * direction, origin + t1 * direction


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
        self.segment_set = None           # name of the set each trace came from
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

    def seed(self, sets=None, rng=None):
        """
        Place the joints and mark the cells and links they occupy.

        Throughgoing sets are laid down first and span the domain. Sets that
        name one of them in ``abuts`` are then cut back so that each trace runs
        from one of its host joints to the next.
        """
        sets = GRANITE_SETS if sets is None else sets
        rng = np.random.default_rng() if rng is None else rng

        lines = {}                       # set name -> list of (origin, dir)
        self.segments = []
        self.segment_set = []

        for js in [s for s in sets if s.is_throughgoing]:
            lines[js.name] = self._lines_for_set(js, rng)
            for origin, d in lines[js.name]:
                seg = _clip_to_box(origin, d, self.lx, self.lz)
                if seg is not None:
                    self.segments.append(seg)
                    self.segment_set.append(js.name)

        for js in [s for s in sets if not s.is_throughgoing]:
            hosts = lines.get(js.abuts)
            if hosts is None:
                raise ValueError(
                    "set %r abuts %r, which is not a throughgoing set here"
                    % (js.name, js.abuts))
            for origin, d in self._lines_for_set(js, rng):
                for seg in self._abut(origin, d, hosts, js, rng):
                    self.segments.append(seg)
                    self.segment_set.append(js.name)

        self._rasterize()
        return self

    def _lines_for_set(self, js, rng):
        """
        Infinite lines of this set, offset along its normal at drawn spacings.

        Seeded across the whole domain diagonal so that no corner is starved --
        an edge with no joint beyond it reads as unfractured rock when it is
        only the boundary of the model.
        """
        d0 = _unit(js.dip_deg)
        n0 = np.array([-d0[1], d0[0]])
        # For a regular set the joints must land on the first and last CELLS,
        # not on the domain edges: a joint placed at x = lx rasterises one
        # column past the end and is clipped, which loses it and breaks the
        # mirror symmetry. Span the cell-centre box instead.
        h = 0.5 * self.dx if (js.spacing_sigma == 0.0 and js.kappa is None) \
            else 0.0
        box = np.array([[h, h], [self.lx - h, h],
                        [h, self.lz - h], [self.lx - h, self.lz - h]])
        t_c, s_c = box @ n0, box @ d0

        regular = js.spacing_sigma == 0.0 and js.kappa is None
        mid = 0.5 * (s_c.min() + s_c.max()) * d0

        if regular:
            # Joints run right out to the walls: the first and last sit ON the
            # boundary, so every block in between is a full spacing across and
            # the edge blocks are no different from the interior ones. Nothing
            # is left over and there is no margin to make symmetric.
            #
            # This is exactly uniform only when the domain is a whole number of
            # spacings. When it is not, the final block is short -- so choose
            # ``n * spacing / dx + 1`` cells across (see ``uniform_grid_shape``).
            span = t_c.max() - t_c.min()
            n = max(int(np.floor(span / js.spacing + 1e-9)), 0)
            offsets = t_c.min() + np.arange(n + 1) * js.spacing
            return [(mid + t * n0, d0) for t in offsets]

        out, t = [], t_c.min()
        while t < t_c.max():
            th = np.deg2rad(js.dip_deg) + rng.vonmises(0.0, js.kappa)
            d = np.array([np.cos(th), np.sin(th)])
            out.append((mid + t * n0, d))
            t += rng.lognormal(np.log(js.spacing), js.spacing_sigma)
        return out

    def _abut(self, origin, d, hosts, js, rng):
        """
        Cut a line into traces that run between host joints and stop there.

        This is the rule that makes a network out of a scatter of lines: a
        younger joint stops at an older one, giving a Y node instead of a free
        tip. Every gap between consecutive hosts is a candidate, occupied with
        probability ``js.density`` -- one trace per line would leave the block
        edges mostly unjointed and the blocks unbounded in this direction.

        Where the line crosses too few hosts to define a gap, it is clipped to
        the domain instead, so the edges are not left bare.
        """
        ts = []
        for h_origin, h_d in hosts:
            den = d[0] * h_d[1] - d[1] * h_d[0]
            if abs(den) < 1e-12:
                continue                                    # parallel
            w = h_origin - origin
            ts.append((w[0] * h_d[1] - w[1] * h_d[0]) / den)
        ts = np.sort(np.array(ts))

        if ts.size < js.spans + 1:
            seg = _clip_to_box(origin, d, self.lx, self.lz)
            return [] if seg is None else [seg]

        out = []
        for i in range(0, ts.size - js.spans, js.spans):
            if js.density < 1.0 and rng.random() > js.density:
                continue
            p0, p1 = origin + ts[i] * d, origin + ts[i + js.spans] * d
            if _clip_to_box(p0, p1 - p0, self.lx, self.lz) is not None:
                out.append((p0, p1))
        return out

    def segments_of(self, set_name):
        """The traces belonging to one named set."""
        return [seg for seg, nm in zip(self.segments, self.segment_set)
                if nm == set_name]

    def _rasterize(self):
        """
        Mark the cells each trace threads, and the links along it.

        A note on symmetry. A centred regular pattern is symmetric as geometry,
        but the raster is only mirror-symmetric if the cell count across the
        domain is **odd**: with an even count and a spacing that is an even
        number of cells, the mirror of the first joint column falls one cell
        short of the last, and every joint carries a one-cell offset. Measured
        on a 20 m section at 5 cm: mean mirror mismatch in the dissolved
        fraction 0.065 at nx = 400, and exactly 0 at nx = 401.
        """
        self.cell = np.zeros((self.nz, self.nx), dtype=bool)
        self.link_v = np.zeros((self.nz - 1, self.nx), dtype=bool)
        self.link_h = np.zeros((self.nz, self.nx - 1), dtype=bool)
        self.trace_length = 0.0

        for p0, p1 in self.segments:
            L = np.hypot(*(p1 - p0))
            if L <= 0.0:
                continue
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

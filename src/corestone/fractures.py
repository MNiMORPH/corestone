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


def periodic_grid_shape(width, depth, dx, spacing):
    """
    Cell counts for a section that wraps left-to-right onto itself.

    Periodic in x, so the width must be a whole number of spacings with no
    joint repeated at the seam: ``n * spacing / dx`` columns. Depth is not
    periodic and keeps the ``+ 1`` that puts a joint on the top and bottom.

    A no-flow wall is not neutral. It forces the lateral flow to vanish there,
    which in a section with subhorizontal joints manufactures a domain-scale
    circulation and a drainage divide down the middle: measured, the centre
    block weathered a third as much as the blocks two in from the walls, and
    the effect grew with the width of the section rather than staying near the
    edges. Wrapping removes the walls entirely.
    """
    per = spacing / dx
    if abs(per - round(per)) > 1e-9:
        raise ValueError(
            "a periodic tiling needs the spacing to be a whole number of "
            "cells; spacing / dx = %.4f. Choose dx = %.4f or %.4f."
            % (per, spacing / np.ceil(per), spacing / np.floor(per)))
    per = int(round(per))
    nx = max(int(np.floor(width / dx / per)), 1) * per
    nz = max(int(np.floor(depth / dx / per)), 1) * per + 1
    return nz, nx


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


def orthogonal_grid(spacing=1.5, density=1.0, rotation=0.0):
    """
    A perfectly regular orthogonal joint grid, optionally rotated.

    No orientation scatter and no spacing variability: every joint exactly
    vertical or exactly horizontal, evenly spaced, so the blocks are identical
    squares of side ``spacing``. The clearest object to teach from, and the
    case whose fracture intensity has an exact answer, P21 = 2 / spacing.

    ``rotation`` turns the whole pair, keeping them perpendicular: 0 gives
    vertical joints cut by horizontal ones, 45 gives a diamond lattice at
    +45/-45.

    A rotated set only tiles a periodic domain exactly when its period ALONG X,
    ``spacing / |n_x|``, divides the width in whole cells. At 45 degrees that
    period is ``spacing * sqrt(2)``, so a spacing of 1.5 m needs a width that
    is a multiple of 2.121 m -- which is not a whole number of 5 cm cells.
    :func:`rotated_grid_shape` picks a spacing and width that do work.
    """
    return conjugate_sets(dip_primary=90.0 - rotation,
                          dip_secondary=-rotation, spacing=spacing,
                          kappa=None, spacing_sigma=0.0, density=density)


def tiling_angles(nx, max_index=None):
    """
    Rotations at which an orthogonal joint pair tiles a periodic domain exactly.

    A pair rotated by theta has along-x periods ``S/cos(theta)`` and
    ``S/sin(theta)``. Both divide the width only when ``tan(theta) = b/a`` for
    integers a and b -- the joint directions have to be rational-slope lattice
    directions -- and, in cells, only when a and b each divide ``nx``. The
    joint spacing that tiles is then ``lx / (k * hypot(a, b))`` for integer k.

    Off these angles the pattern does not close on itself and leaves a seam
    where the joints fail to line up. It is a local defect rather than the
    domain-scale circulation that no-flow walls produce, but it is visible, and
    snapping the angle removes it.

    Returns ``[(degrees, a, b), ...]`` sorted by angle, for 0 to 45 degrees.
    """
    from math import atan2, degrees, gcd
    limit = nx if max_index is None else min(nx, max_index)
    divisors = [d for d in range(1, limit + 1) if nx % d == 0]
    out = {}
    for a in divisors:
        for b in divisors + [0]:
            if b > a or gcd(a, b) != 1:
                continue
            ang = degrees(atan2(b, a))
            if 0.0 <= ang <= 45.0 + 1e-9:
                out.setdefault(round(ang, 6), (ang, a, b))
    return [out[k] for k in sorted(out)]


def tiling_spacings(lx, a, b, low, high, count=8):
    """
    Every tiling spacing at one lattice angle, coarsest first.

    The interactive demos need the whole snapped set to fill a slider, not the
    single nearest value :func:`tiling_spacing` returns. ``low`` and ``high``
    bound which of them are offered and have no defaults on purpose: what
    counts as a sensible joint spacing depends on the section being shown, and
    it is the caller's decision rather than this function's.
    """
    import math
    h = math.hypot(a, b)
    out = []
    for k in range(1, count + 1):
        s = tiling_spacing(lx, a, b, lx / (k * h))
        if low <= s <= high and not any(abs(s - v) < 1e-9 for v in out):
            out.append(s)
    return sorted(out, reverse=True) or [lx / h]


def tiling_spacing(lx, a, b, target):
    """
    The tiling joint spacing nearest ``target`` for the lattice angle (a, b).

    Only ``lx / (k * hypot(a, b))`` tiles, so the spacing is quantised too --
    which is why the spacing slider snaps as well as the angle one.
    """
    import math
    h = math.hypot(a, b)
    k = max(int(round(lx / (target * h))), 1)
    return lx / (k * h)


def rotated_grid_shape(width, depth, dx, x_period, rotation=45.0):
    """
    Cell counts and joint spacing for a rotated grid that still tiles.

    ``x_period`` is the along-x repeat distance, which must be a whole number
    of cells; the joint spacing that produces it is ``x_period * |n_x|``. At
    45 degrees, an x-period of 2.0 m at dx = 0.05 m means 40 cells and a joint
    spacing of 1.414 m.

    Returns ``(nz, nx, spacing)``.
    """
    per = x_period / dx
    if abs(per - round(per)) > 1e-9:
        raise ValueError("x_period must be a whole number of cells; "
                         "x_period/dx = %.4f" % per)
    per = int(round(per))
    nx = max(int(np.floor(width / dx / per)), 1) * per
    nz = max(int(np.floor(depth / dx / per)), 1) * per + 1
    n_x = abs(np.sin(np.deg2rad(90.0 - rotation)))
    return nz, nx, x_period * n_x


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

    def __init__(self, nz, nx, dx, periodic_x=False):
        self.nz = nz                      # rows, increasing downward
        self.nx = nx                      # columns, increasing rightward
        self.dx = dx                      # cell size [m], square cells
        self.periodic_x = periodic_x      # wrap the left and right walls

        self.cell = None                  # (nz, nx) bool: fracture threads it
        self.link_v = None                # (nz-1, nx) bool: vertical links
        self.link_h = None                # (nz, nx-1) bool: horizontal links
        self.link_wrap = None             # (nz,) bool: the wrap-around link
                                          # joining column nx-1 to column 0,
                                          # used only when periodic_x
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

    # ---- construction from an arbitrary network

    @classmethod
    def from_masks(cls, link_v, link_h, dx, periodic_x=False, link_wrap=None,
                   segments=None, segment_set=None):
        """
        Build a network from link masks rather than by seeding one.

        The model reads only ``cell``, ``link_v``, ``link_h``, ``link_wrap``
        and ``periodic_x`` from a network, so a fracture network here is a data
        structure and not a particular algorithm. Anything can supply it: a
        different generator, a traced outcrop photograph, a discrete fracture
        network from another package, a hand-drawn array.

        ``link_v`` is ``(nz-1, nx)`` and ``link_h`` is ``(nz, nx-1)``: True
        where the link between those two cells conducts. A cell counts as
        fractured if any link touching it does.

        The grid must still be a uniform raster of square cells -- ``dx`` is a
        scalar. Non-uniform or unstructured grids would be a different model.

        ``segments`` is only used for drawing. Leave it out and plots simply
        will not overlay the joint traces; nothing in the physics needs it.
        """
        link_v = np.asarray(link_v, dtype=bool)
        link_h = np.asarray(link_h, dtype=bool)
        nz, nx = link_v.shape[0] + 1, link_v.shape[1]
        if link_h.shape != (nz, nx - 1):
            raise ValueError(
                "link_v is (nz-1, nx) = %s, so link_h must be (nz, nx-1) = %s, "
                "not %s" % (link_v.shape, (nz, nx - 1), link_h.shape))

        net = cls(nz, nx, dx, periodic_x=periodic_x)
        net.link_v, net.link_h = link_v, link_h

        cell = np.zeros((nz, nx), dtype=bool)
        cell[:-1, :] |= link_v
        cell[1:, :] |= link_v
        cell[:, :-1] |= link_h
        cell[:, 1:] |= link_h
        net.cell = cell

        if link_wrap is None:
            net.link_wrap = (cell[:, 0] & cell[:, -1]) if periodic_x \
                else np.zeros(nz, dtype=bool)
        else:
            net.link_wrap = np.asarray(link_wrap, dtype=bool)
        if periodic_x:
            net.cell[:, 0] |= net.link_wrap
            net.cell[:, -1] |= net.link_wrap

        # Derived, not measured: without traces there is no trace length, so
        # P21 is estimated from the marked links and should be read as such.
        net.trace_length = float(link_v.sum() + link_h.sum()) * dx
        net.segments = [] if segments is None else list(segments)
        net.segment_set = ([] if segment_set is None
                           else list(segment_set))
        if segments is not None and segment_set is None:
            net.segment_set = ["J1"] * len(net.segments)
        return net

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
            # A set wraps if its normal has any component along x. The period
            # along that normal is lx * |n_x|, NOT lx: only for joints whose
            # normal lies along x (that is, vertical joints) are the two the
            # same. Assuming lx here put the wrong number of joints on any
            # rotated set.
            # The closing-edge special case applies ONLY when the normal lies
            # along x -- vertical joints -- where the first and last joint are
            # the same joint seen from both walls. For any other orientation
            # the extremes of the normal range are different physical joints
            # (they are corners of the domain), and the domain extent along
            # the normal is LARGER than the x-period: at 45 degrees on a 4 m
            # section the extent is 5.66 m against a period of 2.83 m. Counting
            # from the period there left out half the joints.
            n_x = abs(float(n0[0]))
            wraps = self.periodic_x and n_x > 1.0 - 1e-9
            if wraps:
                n = max(int(round(self.lx * n_x / js.spacing)), 1)
                offsets = t_c.min() + np.arange(n) * js.spacing
            else:
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

        if self.periodic_x and ts.size:
            # The span between the outermost hosts runs across the seam, where
            # there is no wall to terminate against. It is one span of the
            # tiling like any other, seen as two pieces because the section is
            # drawn cut open: from the last host to the right edge, and from
            # the left edge to the first host.
            if js.density >= 1.0 or rng.random() <= js.density:
                for a, b in ((ts[-1], None), (None, ts[0])):
                    seg = _clip_to_box(
                        origin + (a if a is not None else ts[0]) * d, d,
                        self.lx, self.lz)
                    if seg is None:
                        continue
                    q0, q1 = seg
                    out.append((q0, origin + b * d) if a is None
                               else (origin + a * d, q1))
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
                # A DIAGONAL step needs a connected route, not two stubs. The
                # grid has no diagonal link, so the trace has to turn a corner:
                # go down first, then across, through the intermediate cell.
                # Marking the outgoing down-link and the outgoing across-link
                # from the SAME cell instead leaves the next trace cell
                # unreachable from either, so a joint at any angle off the axes
                # was a chain of disconnected stubs. Measured at 45 degrees: the
                # joints carried 2x the matrix flux instead of ~1000x.
                if diz[k] != 0:
                    self.link_v[min(iz[k], iz[k + 1]), ix[k]] = True
                if dix[k] != 0:
                    row = iz[k + 1] if diz[k] != 0 else iz[k]
                    self.link_h[row, min(ix[k], ix[k + 1])] = True
                    if diz[k] != 0:
                        self.cell[iz[k + 1], ix[k]] = True   # the corner cell

        # The wrap link is fractured where a joint reaches both walls, which
        # for a periodic tiling means the sub-horizontal set.
        self.link_wrap = self.cell[:, 0] & self.cell[:, -1]

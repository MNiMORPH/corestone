"""
Marching squares, because the demo needs one contour and cannot import a
library that draws them.

``contourpy`` is what bokeh's own ``figure.contour`` uses, and it is not in
the browser build -- the app ships numpy, scipy, bokeh, panel and this
package, and adding a compiled dependency to draw one line is a poor trade.
This is the algorithm, in about forty lines.

WHY INTERPOLATE RATHER THAN OUTLINE THE CELLS. A cell-resolution boundary
could be had for nothing, by eroding the mask and drawing the cells that
survive. It would be wrong in a way worth naming: it draws the contour as a
staircase jumping at cell edges, which asserts that the field is piecewise
constant and that the boundary really does step. The field is a continuum and
the grid is only where it was sampled, so the honest reading between two
samples is a straight line between them. That is what this does.
"""

import numpy as np

#: Which edges each of the sixteen corner-sign patterns joins. Edges are
#: numbered 0 top, 1 right, 2 bottom, 3 left; the case index has bit 0 for the
#: top-left corner and runs clockwise. Cases 5 and 10 are the saddles, where
#: two segments cross the square and the pairing is genuinely ambiguous; they
#: are resolved on the mean of the four corners, below.
_CASES = {
    0: (), 1: ((3, 0),), 2: ((0, 1),), 3: ((3, 1),),
    4: ((1, 2),), 5: ((3, 0), (1, 2)), 6: ((0, 2),), 7: ((3, 2),),
    8: ((2, 3),), 9: ((2, 0),), 10: ((0, 1), (2, 3)), 11: ((2, 1),),
    12: ((1, 3),), 13: ((1, 0),), 14: ((0, 3),), 15: (),
}
#: The saddles again, paired the other way, for when the centre says so.
_SADDLE_FLIP = {5: ((0, 1), (2, 3)), 10: ((3, 0), (1, 2))}


def contour_segments(field, level, dx, dz, x0=0.0, z0=0.0):
    """
    Line segments of ``field == level``, interpolated between cell centres.

    ``field`` is ``(nz, nx)`` on cell centres, the first row at the top. The
    centre of cell ``[i, j]`` sits at ``x0 + (j + 0.5) dx`` and
    ``z0 + (i + 0.5) dz``, matching how the demo places its image.

    Returns ``(xs0, zs0, xs1, zs1)``, four 1-D arrays, one entry per segment.
    Empty arrays if the level is not crossed. Cells that are not finite are
    treated as not crossed, so a masked or divergent field cannot fabricate a
    line.
    """
    g = np.asarray(field, dtype=float) - float(level)
    if g.ndim != 2 or g.shape[0] < 2 or g.shape[1] < 2:
        z = np.empty(0)
        return z, z.copy(), z.copy(), z.copy()
    g = np.where(np.isfinite(g), g, np.nan)

    # The four corners of every square of adjacent cell centres.
    tl, tr = g[:-1, :-1], g[:-1, 1:]
    bl, br = g[1:, :-1], g[1:, 1:]
    ok = np.isfinite(tl) & np.isfinite(tr) & np.isfinite(bl) & np.isfinite(br)

    case = ((tl > 0).astype(np.uint8)
            | ((tr > 0).astype(np.uint8) << 1)
            | ((br > 0).astype(np.uint8) << 2)
            | ((bl > 0).astype(np.uint8) << 3))

    def _cross(a, b):
        """Where along an edge the sign changes, in [0, 1]."""
        d = b - a
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(d != 0.0, -a / d, 0.5)
        return np.clip(t, 0.0, 1.0)

    ii, jj = np.nonzero(ok & (case != 0) & (case != 15))
    if ii.size == 0:
        z = np.empty(0)
        return z, z.copy(), z.copy(), z.copy()

    # Edge crossings, as (x, z) in model units, for the squares that matter.
    xc = x0 + (jj + 0.5) * dx
    zc = z0 + (ii + 0.5) * dz
    t_top = _cross(tl[ii, jj], tr[ii, jj])
    t_bot = _cross(bl[ii, jj], br[ii, jj])
    t_lft = _cross(tl[ii, jj], bl[ii, jj])
    t_rgt = _cross(tr[ii, jj], br[ii, jj])
    pts = {
        0: (xc + t_top * dx, zc),
        1: (xc + dx, zc + t_rgt * dz),
        2: (xc + t_bot * dx, zc + dz),
        3: (xc, zc + t_lft * dz),
    }

    centre = 0.25 * (tl[ii, jj] + tr[ii, jj] + bl[ii, jj] + br[ii, jj])
    cs = case[ii, jj]
    xs0, zs0, xs1, zs1 = [], [], [], []
    for c in np.unique(cs):
        sel = cs == c
        if c in _SADDLE_FLIP:
            # The centre's sign says which way the two branches connect.
            flip = centre[sel] > 0
            for arr, use in ((~flip, _CASES[int(c)]), (flip, _SADDLE_FLIP[int(c)])):
                idx = np.nonzero(sel)[0][arr]
                for a, b in use:
                    xs0.append(pts[a][0][idx]); zs0.append(pts[a][1][idx])
                    xs1.append(pts[b][0][idx]); zs1.append(pts[b][1][idx])
            continue
        idx = np.nonzero(sel)[0]
        for a, b in _CASES[int(c)]:
            xs0.append(pts[a][0][idx]); zs0.append(pts[a][1][idx])
            xs1.append(pts[b][0][idx]); zs1.append(pts[b][1][idx])
    cat = lambda v: np.concatenate(v) if v else np.empty(0)
    return cat(xs0), cat(zs0), cat(xs1), cat(zs1)

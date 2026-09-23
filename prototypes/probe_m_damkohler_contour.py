"""
Probe M: is a Da = 1 contour worth drawing on the demo?

Andy's idea: a student toggle that draws where the Damkohler number crosses 1,
separating rock the water crosses still fresh from rock where it saturates on
arrival. The question this probe settles is whether that line SAYS anything
the figure does not already show.

The worry is specific. Da = spacing / L and L = q / (dM/dt), so L inherits the
flow field, which the joints make wildly uneven -- measured at 1 m spacing it
runs about 7 m in cells a joint touches and 0.24 mm in the matrix between
them. If the contour sits within a cell or two of the joints it is just the
joint trace redrawn, and the demo already draws that.

What would make it worth having is MOVEMENT: weathering raises K_sat, which
spreads the flow, which lengthens L away from the joints, so the Da = 1 line
should migrate outward as the rock opens up. That is a thing to watch, and it
is the answer to "where is the water still doing work".

Reports, at several times: what fraction of the section has Da < 1, and how
far the boundary sits from the nearest joint.
"""

import numpy as np
from scipy.ndimage import distance_transform_edt

import corestone as cs

SPACING, DX, N = 1.0, 0.05, 60


def main():
    net = cs.FractureNetwork(N, N, DX).seed(
        sets=cs.orthogonal_grid(spacing=SPACING))
    w = cs.Weathering(net)
    w.initialize()
    # distance from every cell to the nearest joint-touching cell [m]
    touch = np.zeros((N, N), bool)
    touch[:-1, :] |= net.link_v          # cells either side of a vertical link
    touch[1:, :] |= net.link_v
    touch[:, :-1] |= net.link_h          # ...and of a horizontal one. Measuring
    touch[:, 1:] |= net.link_h           # to vertical joints alone reports a
                                         # cell beside a horizontal joint as far
                                         # from any joint, which it is not.
    dist = distance_transform_edt(~touch, sampling=DX)
    print('%-9s %8s %10s %12s %9s' % (
        'time', 'mean M', 'Da<1 frac', 'boundary [m]', 'cells'))
    last = 0.0
    for years in (0.0, 30e3, 100e3, 300e3):
        if years > last:
            w.run(years - last)
            last = years
        L = w.local_saturation_length()
        fresh = (SPACING / np.maximum(L, 1e-300)) < 1.0
        edge = dist[fresh].max() if fresh.any() else float('nan')
        print('%-9s %8.3f %9.1f%% %12.3f %9.1f' % (
            '%g kyr' % (years / 1e3), float(w.M.mean()),
            100.0 * fresh.mean(), edge, edge / DX))


if __name__ == "__main__":
    main()

"""
The seeded joint network puts water where the design says it should.
"""

import numpy as np
import pytest

from corestone import FractureNetwork, JointSet, GRANITE_SETS


def _net(seed=12345, sets=None, nz=75, nx=100, dx=0.2):
    return FractureNetwork(nz, nx, dx).seed(
        sets=sets, rng=np.random.default_rng(seed))


def test_seeding_is_reproducible_from_the_generator():
    a, b = _net(7), _net(7)
    assert np.array_equal(a.cell, b.cell)
    assert np.array_equal(a.link_v, b.link_v)
    assert len(a.segments) == len(b.segments)


def test_link_arrays_have_the_cell_link_shapes():
    n = _net()
    assert n.cell.shape == (n.nz, n.nx)
    assert n.link_v.shape == (n.nz - 1, n.nx)      # between vertical neighbours
    assert n.link_h.shape == (n.nz, n.nx - 1)      # between horizontal ones


def test_distance_is_zero_on_fractures_and_positive_off_them():
    n = _net()
    d = n.distance_to_fracture()
    assert d[n.cell].max() == 0.0
    assert d[~n.cell].min() > 0.0


def test_a_sheeting_set_brings_the_rock_closer_to_water():
    """
    Two steep conjugate sets cannot bound a block in the vertical, at any
    persistence. Adding a subhorizontal set is what closes it -- the finding in
    design/01-fracture-seeding.md, pinned here so it cannot silently regress.
    """
    steep = [js for js in GRANITE_SETS if abs(js.dip_deg) > 45.0]
    assert len(steep) == 2

    d_steep = _net(sets=steep).distance_to_fracture()
    d_all = _net(sets=GRANITE_SETS).distance_to_fracture()

    assert np.median(d_all) < np.median(d_steep)
    assert np.percentile(d_all, 90) < np.percentile(d_steep, 90)


def test_closer_spacing_puts_rock_nearer_a_fracture():
    wide = [JointSet("A", 75., 20., 3.0), JointSet("B", -75., 20., 3.0)]
    tight = [JointSet("A", 75., 20., 0.75), JointSet("B", -75., 20., 0.75)]
    assert np.median(_net(sets=tight).distance_to_fracture()) \
         < np.median(_net(sets=wide).distance_to_fracture())


def test_padding_does_not_leave_the_edges_starved():
    """
    Seeding only inside the domain leaves edge cells artificially far from any
    fracture, because no fracture can be centred just outside it.
    """
    n = _net()
    d = n.distance_to_fracture()
    edge = np.concatenate([d[0, :], d[-1, :], d[:, 0], d[:, -1]])
    interior = d[1:-1, 1:-1]
    # The edge should not be dramatically worse than the interior.
    assert np.median(edge) < 2.0 * np.median(interior)


def test_p21_is_a_sane_intensity():
    n = _net()
    assert 0.0 < n.p21 < 10.0
    assert n.trace_length > 0.0

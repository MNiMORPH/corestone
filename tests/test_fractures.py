"""
The seeded joint network puts water where the design says it should.
"""

import numpy as np
import pytest

from corestone import (FractureNetwork, JointSet, conjugate_sets,
                       orthogonal_grid, uniform_grid_shape, GRANITE_SETS)

NZ, NX, DX = 75, 100, 0.20


def _net(seed=12345, sets=None):
    return FractureNetwork(NZ, NX, DX).seed(
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


def test_the_default_pair_is_conjugate_at_ninety_degrees():
    assert len(GRANITE_SETS) == 2
    a, b = GRANITE_SETS
    assert abs((a.dip_deg - b.dip_deg) % 180.0) == pytest.approx(90.0)
    assert a.is_throughgoing
    assert b.abuts == a.name


def _on_boundary(pt, n, tol=1e-6):
    return (pt[0] <= tol or pt[0] >= n.lx - tol
            or pt[1] <= tol or pt[1] >= n.lz - tol)


def test_throughgoing_traces_run_boundary_to_boundary():
    """
    A primary joint crosses the whole section -- both of its tips are on the
    domain edge, never in intact rock. That is what makes the network connect:
    an earlier generator drew both sets as free segments from a power-law
    length distribution, and the joints came out shorter than the gaps between
    them, so nothing linked up.

    Note it is boundary-to-boundary, not top-to-bottom: with orientation
    scatter a near-vertical joint seeded near a side exits through that side.
    """
    n = _net(sets=[JointSet("J1", dip_deg=90.0, spacing=1.5)])
    assert len(n.segments) >= 5
    for p0, p1 in n.segments:
        assert _on_boundary(p0, n) and _on_boundary(p1, n)


def test_abutting_traces_terminate_on_a_host_joint():
    """
    The abutting rule is what turns a scatter of lines into a network: a
    younger joint stops at an older one, giving a Y node rather than a free
    tip. Checked here by geometry -- each secondary endpoint should sit on a
    primary trace, or on the domain edge where it was clipped.
    """
    n = _net(sets=conjugate_sets(90.0, 0.0))
    primary = n.segments_of("J1")
    assert len(primary) >= 3

    def on_a_primary_or_edge(pt, tol=1e-6):
        if _on_boundary(pt, n, tol):
            return True
        for p0, p1 in primary:
            d = p1 - p0
            L2 = d @ d
            t = np.clip(((pt - p0) @ d) / L2, 0.0, 1.0)
            if np.hypot(*(p0 + t * d - pt)) < 1e-6:
                return True
        return False

    secondary = n.segments_of("J2")
    assert len(secondary) >= 3
    for p0, p1 in secondary:
        assert on_a_primary_or_edge(p0) and on_a_primary_or_edge(p1)


def test_a_conjugate_pair_beats_a_single_set():
    """
    One set of parallel joints cannot bound a block: it leaves corridors
    unbroken in its own direction. The conjugate partner is what closes them.
    """
    one = [JointSet("J1", dip_deg=90.0, spacing=1.5)]
    d_one = _net(sets=one).distance_to_fracture()
    d_two = _net(sets=conjugate_sets(90.0, 0.0, spacing=1.5)).distance_to_fracture()
    assert np.percentile(d_two, 90) < np.percentile(d_one, 90)


def test_closer_spacing_puts_rock_nearer_a_fracture():
    tight = _net(sets=conjugate_sets(90.0, 0.0, spacing=0.8))
    wide = _net(sets=conjugate_sets(90.0, 0.0, spacing=3.0))
    assert np.percentile(tight.distance_to_fracture(), 90) \
         < np.percentile(wide.distance_to_fracture(), 90)
    assert tight.p21 > wide.p21


def test_abutting_a_set_that_is_not_there_is_an_error():
    bad = [JointSet("J2", dip_deg=0.0, spacing=1.5, abuts="nope")]
    with pytest.raises(ValueError, match="not a throughgoing set"):
        _net(sets=bad)


def test_p21_is_a_sane_intensity():
    n = _net()
    assert 0.0 < n.p21 < 10.0
    assert n.trace_length > 0.0


def test_density_controls_how_many_block_edges_are_jointed():
    """
    One abutting trace per line leaves the block edges mostly unjointed, so the
    blocks are bounded in one direction only. Every gap is a candidate, taken
    with probability `density`.
    """
    full = _net(sets=conjugate_sets(90.0, 0.0, density=1.0))
    sparse = _net(sets=conjugate_sets(90.0, 0.0, density=0.3))
    n_full = sum(1 for s in full.segment_set if s == "J2")
    n_sparse = sum(1 for s in sparse.segment_set if s == "J2")
    assert n_full > 3 * n_sparse
    assert full.p21 > sparse.p21
    # Bounding the blocks in both directions brings the far corners closer.
    assert full.distance_to_fracture().max() \
         < sparse.distance_to_fracture().max()


def test_a_full_orthogonal_grid_has_the_analytic_intensity():
    """
    P21 for a grid of lines at spacing S in two orthogonal directions is 2/S.
    A check on the generator that does not depend on any of our own machinery.
    """
    S = 1.5
    n = _net(sets=conjugate_sets(90.0, 0.0, spacing=S, density=1.0))
    assert n.p21 == pytest.approx(2.0 / S, rel=0.15)


def test_the_orthogonal_grid_is_exactly_axis_aligned():
    """No orientation scatter: every joint exactly vertical or exactly
    horizontal, so one component of every trace vector is zero."""
    n = _net(sets=orthogonal_grid(spacing=1.5))
    for p0, p1 in n.segments:
        d = p1 - p0
        assert min(abs(d[0]), abs(d[1])) < 1e-9


def test_the_orthogonal_grid_is_evenly_spaced():
    """Exact spacing, not lognormal: consecutive parallel joints are one
    spacing apart to machine precision."""
    S = 1.5
    n = _net(sets=orthogonal_grid(spacing=S))
    xs = np.unique(np.round([p0[0] for p0, p1 in n.segments_of("J1")], 9))
    assert len(xs) >= 5
    assert np.allclose(np.diff(xs), S)


def test_the_orthogonal_grid_does_not_depend_on_the_generator():
    """With no scatter, no spacing variability and every gap filled, the
    network is fully determined -- the random generator is not consulted."""
    a = _net(seed=1, sets=orthogonal_grid(spacing=1.5))
    b = _net(seed=999, sets=orthogonal_grid(spacing=1.5))
    assert np.array_equal(a.cell, b.cell)
    assert len(a.segments) == len(b.segments)


def test_a_regular_grid_is_uniform_right_out_to_the_walls():
    """
    A joint sits on the first and last cell, and every block between them is a
    full spacing across -- no odd strip at the edges. An earlier version
    stepped from one edge until the domain ran out and piled the whole
    remainder on the far side: at 1.5 m spacing across 20 m the margins came
    out 1.25 m and 0.75 m, and the flow inherited that asymmetry.
    """
    nz, nx = uniform_grid_shape(20.0, 15.0, 0.05, 1.5)
    n = FractureNetwork(nz, nx, 0.05).seed(sets=orthogonal_grid(1.5),
                                           rng=np.random.default_rng(1))
    xs = np.array(sorted({round(p0[0], 6) for p0, _ in n.segments_of("J1")}))
    assert np.allclose(np.diff(xs), 1.5)              # uniform spacing
    assert xs.min() == pytest.approx(0.5 * n.dx)      # a joint on each wall
    assert xs.max() == pytest.approx(n.lx - 0.5 * n.dx)
    assert np.array_equal(n.cell, n.cell[:, ::-1])
    assert np.array_equal(n.cell, n.cell[::-1, :])


def test_a_periodic_tiling_refuses_a_spacing_that_is_not_whole_cells():
    """
    1.5 m at 20 cm is 7.5 cells. Rounding it silently would tile at 1.6 m while
    the caller believed 1.5 m -- exactly the kind of quiet substitution that
    turns into a wrong result nobody can trace.
    """
    from corestone import periodic_grid_shape
    with pytest.raises(ValueError, match="whole number of cells"):
        periodic_grid_shape(20.0, 15.0, 0.20, 1.5)
    assert periodic_grid_shape(20.0, 15.0, 0.05, 1.5) == (301, 390)

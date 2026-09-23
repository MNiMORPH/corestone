"""
The contour is drawn by us, so it is tested by us.

Each test names the claim it makes about :func:`corestone.contour.contour_segments`.
The analytic cases matter more than they look: a marching-squares table with a
transposed row still produces a plausible-looking tangle of line segments, and
the only cheap way to catch that is to contour a field whose answer is known
exactly and check the coordinates rather than the picture.
"""

import numpy as np
import pytest

from corestone.contour import contour_segments

DX = DZ = 0.1
N = 10


def _linear_x(nz=N, nx=N, dx=DX):
    """``f = x`` sampled at cell centres."""
    return np.tile((np.arange(nx) + 0.5) * dx, (nz, 1))


def test_a_linear_field_contours_to_a_straight_line_at_the_exact_level():
    x0, z0, x1, z1 = contour_segments(_linear_x(), 0.5, DX, DZ)
    assert x0.size > 0
    assert np.allclose(x0, 0.5, atol=0.0, rtol=0.0)
    assert np.allclose(x1, 0.5, atol=0.0, rtol=0.0)


def test_the_line_lands_between_samples_and_not_on_one():
    """
    0.47 falls between the cell centres at 0.45 and 0.55, and the contour has
    to land there. A staircase implementation would put it on one or the
    other, which is the whole reason this is interpolated.
    """
    x0, _, x1, _ = contour_segments(_linear_x(), 0.47, DX, DZ)
    assert np.allclose(x0, 0.47) and np.allclose(x1, 0.47)


def test_a_level_outside_the_field_gives_nothing():
    for level in (-1.0, 99.0):
        x0, z0, x1, z1 = contour_segments(_linear_x(), level, DX, DZ)
        assert x0.size == z0.size == x1.size == z1.size == 0


def test_a_diagonal_field_contours_along_the_diagonal():
    """Exercises the corner cases the axis-aligned test never reaches."""
    i = (np.arange(N) + 0.5) * DZ
    f = i[:, None] + i[None, :]
    x0, z0, x1, z1 = contour_segments(f, 1.0, DX, DZ)
    assert x0.size > 0
    assert np.allclose(x0 + z0, 1.0) and np.allclose(x1 + z1, 1.0)


def test_non_finite_cells_do_not_fabricate_a_line():
    f = _linear_x()
    f[4:6, :] = np.nan
    x0, z0, _, _ = contour_segments(f, 0.5, DX, DZ)
    # nothing is drawn in the rows whose corners are not all finite
    assert not np.any((z0 > 0.35) & (z0 < 0.65))


def test_the_saddle_is_resolved_and_produces_two_branches():
    """
    A 2x2 with diagonal corners high is the ambiguous case. It must give two
    segments, not one and not four.
    """
    f = np.array([[1.0, -1.0], [-1.0, 1.0]])
    x0, _, _, _ = contour_segments(f, 0.0, DX, DZ)
    assert x0.size == 2


def test_the_offset_places_the_first_sample_at_half_a_cell():
    """The image starts at the domain corner, so a cell centre is half in."""
    x0, _, _, _ = contour_segments(_linear_x(), 0.5, DX, DZ, x0=10.0)
    assert np.allclose(x0, 10.5)


@pytest.mark.parametrize("shape", [(1, 5), (5, 1), (1, 1)])
def test_a_field_too_small_to_have_a_square_gives_nothing(shape):
    x0, _, _, _ = contour_segments(np.ones(shape), 0.5, DX, DZ)
    assert x0.size == 0

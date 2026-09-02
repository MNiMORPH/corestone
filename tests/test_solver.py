"""
The solver's shortcuts, and the properties they depend on.

Every routine here is about *cost*, not physics -- but each one is fast only
because the matrix has some property, and a property assumed is a property that
will eventually stop holding. So each shortcut is paired with a test of the
thing it relies on, and with a test that it still gives the same answer as the
plain method it replaced.
"""

import numpy as np
import pytest
import scipy.sparse as sp

from corestone import FractureNetwork, Weathering, orthogonal_grid, YEAR


def _model(nz=41, nx=40, dx=0.10, spacing=0.8, periodic=True):
    net = FractureNetwork(nz, nx, dx, periodic_x=periodic).seed(
        sets=orthogonal_grid(spacing), rng=np.random.default_rng(12345))
    m = Weathering(net)
    m.set_infiltration(0.30 / YEAR)
    return m.initialize()


def _structurally_symmetric(A):
    """True if the sparsity PATTERN is symmetric; values need not be."""
    S = (A.tocsr() != 0).astype(np.int8)
    return (S - S.T).nnz == 0


@pytest.mark.parametrize("periodic", [True, False])
def test_the_transport_operator_is_structurally_symmetric(periodic):
    """
    The precondition for ordering by minimum degree on ``A + A.T``.

    Every link contributes in both directions -- upwind advection sends the
    flux to one neighbour and takes the diagonal from the other, and diffusion
    is symmetric outright -- so the pattern is symmetric even though the values
    are not. An upwind stencil that reached only downstream would break this
    silently, costing fill rather than correctness, which is exactly the kind
    of regression no physics test would catch.
    """
    m = _model(periodic=periodic)
    assert _structurally_symmetric(m._transport_operator())


def test_the_transport_operator_is_not_symmetric_in_its_values():
    """
    The companion to the test above, so that "symmetric" is never read as more
    than it is: advection is directional, so ``MMD_ATA`` and a symmetric solver
    would both be wrong here.
    """
    A = _model()._transport_operator().tocsr()
    assert (A - A.T).nnz > 0


def test_every_cell_carries_a_diagonal_entry():
    """
    The reaction term is added to the diagonal of the transport operator every
    step. That is only possible in place if the diagonal is already present in
    the pattern -- true because every cell has at least one link, and every
    link puts a positive coefficient on both of its cells' diagonals.
    """
    m = _model()
    A = m._transport_operator().tocsc()
    n = m.nz * m.nx
    assert (A.diagonal() != 0.0).all()
    assert A.getnnz() >= n


def test_the_in_place_step_matrix_equals_the_sparse_addition_exactly():
    """
    ``_step_matrix`` writes the reaction term into a matrix allocated once,
    rather than adding ``sp.diags(r dx^2)`` to the transport operator. Same two
    floats, added in the same order, so this is an equality and not a
    tolerance -- if it ever needs a tolerance, something else has changed.
    """
    m = _model()
    r = m.reaction_coefficient
    fast = m._step_matrix(r).tocsc()
    slow = (m._transport_operator()
            + sp.diags((r * m.dx * m.dx).ravel())).tocsc()
    assert (fast - slow).nnz == 0
    assert np.array_equal(fast.toarray(), slow.toarray())


def test_the_step_matrix_is_rebuilt_when_the_transport_coefficients_change():
    """
    The reused container is invalidated with the operator it mirrors. A cache
    that outlived its inputs is a defect this model has already had once --
    setting ``D_molecular = 0`` silently did nothing -- so the second cache
    gets the same test as the first.
    """
    m = _model()
    r = m.reaction_coefficient
    before = m._step_matrix(r).copy()
    m.D_molecular = 0.0
    m.dispersivity = 0.0
    after = m._step_matrix(r)
    assert (before - after).nnz > 0
    assert after.shape == before.shape

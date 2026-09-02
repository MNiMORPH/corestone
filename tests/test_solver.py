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


def test_the_reused_solution_does_not_change_the_answer():
    """
    The warm start is a guess, not an approximation. Whatever it starts from,
    the iteration is required to drive ``||b - A x||`` below the tolerance, so
    the converged field must match a direct factorise-and-solve of the very
    same system. Checked after the rock has weathered a while, when the guess
    and the answer have had time to separate.
    """
    import scipy.sparse.linalg as spl
    m = _model()
    m.run(years=40e3)
    r = m.reaction_coefficient
    warm = m.solve_solute(r)
    assert m._x is not None                       # the guess was actually kept

    A = m._step_matrix(r)
    b = np.broadcast_to(r * m.dx * m.dx, (m.nz, m.nx)).ravel()
    cold = np.clip(spl.splu(A.tocsc()).solve(b.copy()), 0.0, 1.0)
    assert np.abs(warm.ravel() - cold).max() < 1e-9


def test_the_residual_of_the_returned_field_meets_the_tolerance():
    """
    Stated as the solver's own contract, rather than inferred from agreement
    with another solver: the field it hands back solves the system it was
    given, to the tolerance it was given.
    """
    m = _model()
    m.run(years=40e3)
    r = m.reaction_coefficient
    x = m.solve_solute(r)
    A = m._step_matrix(r)
    b = np.broadcast_to(r * m.dx * m.dx, (m.nz, m.nx)).ravel()
    res = np.linalg.norm(b - A @ x.ravel()) / np.linalg.norm(b)
    assert res < 1e-8


@pytest.mark.parametrize("periodic", [True, False])
def test_the_flow_operator_is_structurally_symmetric(periodic):
    """
    The second matrix ``ORDERING`` speaks for. This one is symmetric in its
    values too -- it is a conductance Laplacian, and the base boundary adds to
    the diagonal only -- so the claim is stronger here than for transport.
    """
    A, _ = _model(periodic=periodic).flow_operator()
    assert _structurally_symmetric(A)
    assert np.abs((A - A.T).data).max() == 0.0 if (A - A.T).nnz else True


def test_the_base_conductance_reaches_the_matrix_through_the_triplets():
    """
    The base boundary is assembled as one more triplet and summed by COO,
    rather than written into an already-built matrix. The diagonal of the last
    row must therefore exceed what the links alone would put there, by exactly
    the base conductance.
    """
    m = _model()
    A, rhs = m.flow_operator()
    d = A.diagonal().reshape(m.nz, m.nx)
    # links only: rebuild the same row without the base term
    net = m.network
    kv = np.where(net.link_v, m.k_fracture, m.k_matrix)
    kh = np.where(net.link_h, m.k_fracture, m.k_matrix)
    links = kv[-1, :].copy()                       # from the row above
    links[:-1] += kh[-1, :]
    links[1:] += kh[-1, :]
    if net.periodic_x:
        kw = np.where(net.link_wrap, m.k_fracture, m.k_matrix)
        links[-1] += kw[-1]
        links[0] += kw[-1]
    assert np.allclose(d[-1, :] - links, m._k_base, rtol=1e-12)
    assert np.allclose(rhs.reshape(m.nz, m.nx)[-1, :],
                       m._k_base * m._h_base, rtol=1e-12)


def test_run_lands_on_the_time_it_was_asked_for():
    """
    ``run`` used to step past its target by up to one step. Harmless in a
    single run and poisonous in any comparison: two settings would be compared
    at two different model times, and the difference read as the error of the
    coarser one. It put a floor of ~1e-2 under a convergence study that should
    have gone to zero, and made that study look non-monotone.

    Checked with long steps, since that is when the overshoot was large: with
    a slack limiter a run to 30 kyr used to stop at 32.5 kyr.
    """
    for dx_max in (0.5, 0.05):
        m = _model()
        m.dt_max = 50000.0 * YEAR
        m.dx_max = dx_max
        m.run(years=30e3)
        assert m.t == pytest.approx(30e3 * YEAR, rel=1e-12)

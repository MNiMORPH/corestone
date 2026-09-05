"""
The solver's shortcuts, and the properties they depend on.

Every routine here is about *cost*, not physics -- but each one is fast only
because the matrix has some property, and a property assumed is a property that
will eventually stop holding. So each shortcut is paired with a test of the
thing it relies on, and with a test that it still gives the same answer as the
plain method it replaced.
"""

import inspect

import numpy as np
import scipy.sparse.linalg as spl

from corestone.weathering import ORDERING
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
    m.grain_size = 0.0
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
    # The model's own right-hand side: the claim here is about the SOLVER
    # -- that the field it returns solves the system it was handed -- and
    # that system's source depends on which reaction is driving.
    b = m._solute_source(r)
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
    # The model's own right-hand side: the claim here is about the SOLVER
    # -- that the field it returns solves the system it was handed -- and
    # that system's source depends on which reaction is driving.
    b = m._solute_source(r)
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
    # k_matrix_at_T, not k_matrix: the conductivities carry the viscosity at
    # the working temperature, and a fresh model has M = 1 everywhere so
    # k(M) is exactly that. Using the uncorrected value here left this
    # comparison wrong by k_matrix - k_matrix_at_T, which is 9.45e-11 -- and
    # np.allclose's default atol of 1e-8 swallowed it whole.
    kv = np.where(net.link_v, m.k_fracture, m.k_matrix_at_T)
    kh = np.where(net.link_h, m.k_fracture, m.k_matrix_at_T)
    links = kv[-1, :].copy()                       # from the row above
    links[:-1] += kh[-1, :]
    links[1:] += kh[-1, :]
    if net.periodic_x:
        kw = np.where(net.link_wrap, m.k_fracture, m.k_matrix_at_T)
        links[-1] += kw[-1]
        links[0] += kw[-1]
    assert np.allclose(d[-1, :] - links, m._k_base, rtol=1e-12, atol=0.0)
    assert np.allclose(rhs.reshape(m.nz, m.nx)[-1, :],  # atol below
                       m._k_base * m._h_base, rtol=1e-12, atol=0.0)


def test_the_step_control_makes_the_error_a_dial():
    """
    The property that justifies replacing dx_max: tightening ``c_drift_max``
    must reduce the error, every time, on the way to the converged answer.

    dx_max could not do this. It bounds the change in M in whichever cell is
    dissolving fastest, and which cell that is jumps about, so the error is not
    monotone in it -- measured on the 3 m section at dt_max = 50 kyr, 0.05 gave
    1.18e-4, 0.10 gave 2.66e-5 and 0.20 gave 2.27e-4. A knob you cannot turn
    predictably is not a control, and no error budget can be set against one.
    """
    # 300 kyr, not 30. At 30 kyr the two coarsest budgets give bit-identical
    # answers under the oxidation driver, because another cap -- dx_max or
    # dt_max -- binds before the drift ever reaches 0.03, so the comparison is
    # degenerate rather than non-monotone. Measured: at 30 kyr the errors are
    # 6.84e-03, 6.84e-03, 3.16e-03, 4.88e-04, and at 300 kyr 1.19e-01,
    # 1.04e-01, 2.87e-02, 6.95e-03. A control can only be shown to be a dial
    # over a range where it is the thing doing the controlling.
    YEARS = 300e3

    def at(drift):
        m = _model()
        m.c_drift_max = drift
        m.run(years=YEARS)
        return m.M

    m = _model()
    m.c_drift_max = 3e-4
    m.run(years=YEARS)
    ref = m.M

    errors = [np.abs(at(d) - ref).max() for d in (0.10, 0.03, 0.01, 0.003)]
    assert all(a > b for a, b in zip(errors, errors[1:])), errors
    # ...and it is a real range. The threshold was 20 while the model
    # dissolved; under oxidation the same four budgets span 16.9x
    # (1.086e-01, 9.971e-02, 2.854e-02, 6.420e-03 against the 3e-4
    # reference),
    # so the number is loosened to what the claim actually needs -- that
    # turning the dial moves the error by an order of magnitude -- rather
    # than kept at a value calibrated on the other reaction.
    assert errors[0] > 10.0 * errors[-1], errors


def test_an_explicit_step_overrides_the_drift_control():
    """``update(dt=...)`` means that step, whatever the controller wants."""
    m = _model()
    m.run(years=5e3)
    want = 137.0 * YEAR
    assert m.update(dt=want) == want


def test_run_lands_on_the_time_it_was_asked_for():
    """
    ``run`` used to step past its target by up to one step. Harmless in a
    single run and poisonous in any comparison: two settings would be compared
    at two different model times, and the difference read as the error of the
    coarser one. It put a floor of ~1e-2 under a convergence study that should
    have gone to zero, and made that study look non-monotone.

    Checked with long steps, since that is when the overshoot was large: at
    30 kyr with a slack drift budget the model used to stop at 32.5 kyr.
    """
    for drift in (0.03, 0.003):
        m = _model()
        m.c_drift_max = drift
        m.run(years=30e3)
        assert m.t == pytest.approx(30e3 * YEAR, rel=1e-12)


def test_a_step_shortened_to_land_on_the_target_does_not_shrink_the_next_one():
    """
    The control remembers the step it WANTED, not the one it was allowed. If a
    short final step were fed back, ``run`` called twice in a row would crawl
    where a single call would not, and the answer would depend on how the run
    was chopped into calls.
    """
    m = _model()
    m.run(years=20e3)
    tiny = 10.0 * YEAR
    step = m.update(dt_limit=tiny)                # capped hard by the caller
    assert step == tiny
    assert m._dt > 100.0 * tiny                   # the control is not fooled


def test_the_iterative_tolerance_is_named_the_way_this_scipy_names_it():
    """
    SciPy called it ``tol`` until 1.12, ``rtol`` from 1.12, and removed ``tol``
    in 1.14. Pyodide ships 1.14. A hardcoded ``tol=`` therefore worked on the
    workstation (SciPy 1.11) and raised TypeError in the browser -- inside a
    web worker, where the traceback never reached the page console, so the demo
    simply refused to advance and said nothing.

    Two assertions, because either alone would have missed it: the name the
    module resolved is one this SciPy actually accepts, and a real solve gets
    past the first step, which is where the iterative path is first taken.
    """
    import scipy.sparse.linalg as spl
    from corestone.weathering import _RTOL

    assert _RTOL in inspect.signature(spl.bicgstab).parameters

    A = sp.eye(8, format="csc") * 2.0
    b = np.ones(8)
    x, info = spl.bicgstab(A, b, atol=0.0, **{_RTOL: 1e-10})
    assert info == 0 and np.allclose(x, 0.5, rtol=1e-9, atol=0.0)

    m = _model()
    m.update()                                 # first step: direct
    m.update()                                 # second: the iterative path
    assert m.t > 0.0


def test_the_flow_matrix_is_symmetric_so_cg_applies():
    """Warm-started conjugate gradients is only legitimate because the
    conductance Laplacian is symmetric with a positive diagonal. If an edit
    ever breaks that -- an asymmetric boundary term, say -- CG would converge
    to the wrong thing quietly."""
    m = _model()
    A, _ = m.flow_operator()
    assert abs(A - A.T).max() == 0.0
    assert (A.diagonal() > 0).all()


def test_the_head_factorisation_is_kept_across_solves():
    """
    The point of the warm start. Re-solving the head is half the run at a
    converged flow_tolerance, almost all of it factorisation; keeping the
    previous factorisation as a preconditioner replaces most of those with
    back-substitutions, and that is what makes flow_tolerance = 0.01
    affordable.

    Bites: pinning max_head_iterations to 0 forces a refactorisation every
    time, which is the behaviour this replaced.
    """
    import scipy.sparse.linalg as spl_
    def count(warm):
        n = {"i": 0}
        orig = spl_.splu
        spl_.splu = lambda *a, **k: (n.__setitem__("i", n["i"] + 1), orig(*a, **k))[1]
        try:
            m = _model()
            if not warm:
                m.max_head_iterations = 0
            m.initialize()
            for _ in range(6):
                m.solve_flow()
        finally:
            spl_.splu = orig
        return n["i"]
    warm, always = count(True), count(False)
    assert warm < always, (warm, always)
    assert always >= 7, always            # one per solve, plus initialize


def test_the_warm_started_head_matches_a_direct_solve():
    """A preconditioner may be stale; the ANSWER may not be."""
    m = _model().initialize()
    A, b = m.flow_operator()
    direct = spl.splu(A, permc_spec=ORDERING).solve(b)
    m.solve_flow()
    for _ in range(4):
        m.M *= 0.97
        m.solve_flow()
    A, b = m.flow_operator()
    direct = spl.splu(A, permc_spec=ORDERING).solve(b)
    warm = m._solve_head(A, b)
    assert np.abs(warm - direct).max() / np.abs(direct).max() < 1e-9

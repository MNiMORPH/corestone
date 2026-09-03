"""
No equation may sit in a package docstring without a check that it is the one
the code implements.

``test_stated_equations.py`` verifies the equations we thought to verify. This
file is what stops that set from being the set we happened to remember. It
reads every docstring in ``src/corestone/``, pulls out every displayed
equation, and fails if one of them is not registered below against a test that
exists.

The ledger is deliberately verbatim: **rewording an equation breaks this
test**. That is friction on purpose. A changed equation is one whose check must
be read again, and the whole reason this file exists is that a docstring
drifted away from its code and nothing noticed for six revisions.

Design documents are out of scope. They record what was measured against the
code of the day, not a live claim about the code as it stands.
"""

import ast
import pathlib
import importlib.util
import re
import sys
import textwrap
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parent.parent / "src" / "corestone"

#: A displayed block counts as an equation if it contains a bare relation.
_RELATION = re.compile(r"(?<![<>=!])=(?!=)|->")

#: Pins the extractor itself. If this number moves without the ledger moving,
#: the extractor has broken and would otherwise report a quiet all-clear.
EXPECTED_BLOCKS = 11
EXPECTED_DISTINCT = 9

LEDGER = {
    "R = k(T) * A * (1 - C / C_eq)":
        "test_the_dissolution_rate_per_unit_volume_does_not_depend_on_the_flux",

    "saturation_length = q * C_eq / (k(T) * A)":
        "test_the_saturation_length_is_proportional_to_the_local_flux",

    "div(q c) - div(D grad c) = r (1 - c), r = k A / C_eq":
        "test_the_solved_concentration_satisfies_the_stated_cell_balance",

    "div( K grad H ) = 0, H = psi - d (d is depth, positive down)":
        "test_the_head_field_satisfies_the_darcy_equation_cell_by_cell",

    "D = D_molecular / tortuosity + dispersivity * |v|":
        "test_the_transport_coefficient_is_molecular_plus_dispersive",

    "sum_out f c_i - sum_in f c_j + sum_links D (c_i - c_j) + r dx^2 c_i "
    "= r dx^2":
        "test_the_solved_concentration_satisfies_the_stated_cell_balance",

    "d(M/M0)/dt = - r (1 - c) / tau":
        "test_what_the_rock_loses_is_what_the_water_carries_out_of_the_base",

    "M(t + dt) = M(t) exp(-lambda dt), lambda = (r / M) (1 - c) / tau":
        "test_the_rock_is_integrated_exactly_over_a_step_with_c_held",

    "k(M) = k_matrix^M * k_weathered^(1 - M)":
        "test_the_matrix_conducts_better_as_it_dissolves",
}


def _docstrings(path):
    """(qualified name, docstring) for every module, class and function."""
    tree = ast.parse(path.read_text())
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef,
                             ast.AsyncFunctionDef)):
            doc = ast.get_docstring(node, clean=False)
            if doc:
                yield getattr(node, "name", "<module>"), doc


def _display_blocks(doc):
    """Runs of consecutive lines indented past the docstring's own margin."""
    lines = textwrap.dedent(doc).splitlines()
    margin = min((len(l) - len(l.lstrip()) for l in lines if l.strip()),
                 default=0)
    block = []
    for line in lines:
        if line.strip() and (len(line) - len(line.lstrip())) > margin:
            block.append(line.strip())
        else:
            if block:
                yield block
            block = []
    if block:
        yield block


def stated_equations():
    """[(file, owner, equation text)] for the whole package."""
    found = []
    for path in sorted(SRC.glob("*.py")):
        for owner, doc in _docstrings(path):
            for block in _display_blocks(doc):
                if any(_RELATION.search(line) for line in block):
                    found.append((path.name, owner,
                                  " ".join(" ".join(block).split())))
    return found


def _stated_equations_module():
    """Load the companion test module by path, not by whatever is on sys.path."""
    import importlib.util
    path = pathlib.Path(__file__).with_name("test_stated_equations.py")
    spec = importlib.util.spec_from_file_location("_stated_equations", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_extractor_still_finds_the_equations_it_is_supposed_to():
    """
    A checker that silently finds nothing passes forever. Pin the count, so
    that breaking the extractor is a failure and not a quiet all-clear.
    """
    found = stated_equations()
    assert len(found) == EXPECTED_BLOCKS, [f[2] for f in found]
    assert len({eq for _, _, eq in found}) == EXPECTED_DISTINCT


@pytest.mark.parametrize("where,owner,equation",
                         stated_equations(),
                         ids=lambda v: str(v)[:60])
def test_every_stated_equation_is_registered_against_a_test(where, owner,
                                                            equation):
    """
    An equation in a docstring is a claim about the code beneath it. Either it
    is checked, or it should not be stated as though it were.
    """
    assert equation in LEDGER, (
        "%s (%s) states an equation with no check registered for it:\n"
        "    %s\n"
        "Write the check in tests/test_stated_equations.py and add the "
        "equation to LEDGER in this file. If the equation has merely been "
        "reworded, update the LEDGER key -- and re-read the check while you "
        "are there, which is the point of making you come here."
        % (where, owner, equation))


@pytest.mark.parametrize("equation,test_name", sorted(LEDGER.items()),
                         ids=lambda v: str(v)[:60])
def test_every_registered_test_exists_and_names_its_equation(equation,
                                                             test_name):
    """
    The other direction: a ledger entry pointing at a test that was renamed or
    deleted would leave the equation unguarded while still looking guarded.
    """
    stated = _stated_equations_module()
    assert hasattr(stated, test_name), (
        "LEDGER maps\n    %s\nto %s(), which does not exist in "
        "tests/test_stated_equations.py" % (equation, test_name))
    assert callable(getattr(stated, test_name))


def test_no_ledger_entry_is_stale():
    """A registered equation that no longer appears in the source is a sign the
    prose moved and the check did not follow it."""
    present = {eq for _, _, eq in stated_equations()}
    stale = sorted(set(LEDGER) - present)
    assert not stale, (
        "LEDGER registers equations that no longer appear in any package "
        "docstring: %s" % stale)

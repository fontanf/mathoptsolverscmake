import bz2
import gc
import gzip
import os

import numpy as np
import pytest

import mathoptmodel as mom

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "data", "miplib3")
P0033 = os.path.join(DATA_DIR, "p0033.mps")

RANGES_MPS = """NAME          RANGETEST
ROWS
 N  COST
 L  LIM1
 G  LIM2
 E  LIM3
COLUMNS
    MARKER                 'MARKER'                 'INTORG'
    X1        COST         1.0   LIM1         1.0
    X1        LIM2         1.0   LIM3         1.0
    MARKER                 'MARKER'                 'INTEND'
    X2        COST         2.0   LIM1         3.0
RHS
    RHS       LIM1        10.0   LIM2         3.0
    RHS       LIM3         5.0
RANGES
    RNG       LIM1         4.0   LIM2         6.0
    RNG       LIM3         2.0
BOUNDS
 UP BND       X1           1.0
 MI BND       X2
ENDATA
"""

LP = """Maximize
 obj: x + 2 y + 3 z
Subject To
 c1: x + y + z <= 10
 c2: x - y >= -3
Bounds
 0 <= z <= 5
General
 z
Binary
 x
End
"""


def test_read_mps_from_bytes():
    model = mom.read_mps_bytes(RANGES_MPS.encode())
    assert model.number_of_variables == 2
    assert model.number_of_constraints == 3
    assert model.objective_direction == mom.ObjectiveDirection.Minimize
    assert model.variables_names == ["X1", "X2"]
    assert model.constraints_names == ["LIM1", "LIM2", "LIM3"]
    # MPS RANGES semantics: L-row [rhs-|r|, rhs], G-row [rhs, rhs+|r|], E-row (r>=0) [rhs, rhs+r].
    np.testing.assert_array_equal(model.constraints_lower_bounds, [6, 3, 5])
    np.testing.assert_array_equal(model.constraints_upper_bounds, [10, 9, 7])
    np.testing.assert_array_equal(model.variables_lower_bounds, [0, -np.inf])
    np.testing.assert_array_equal(model.variables_upper_bounds, [1, np.inf])
    assert list(model.variables_types) == [mom.VariableType.Binary, mom.VariableType.Continuous]
    assert model.constraint_sense(0) == mom.ConstraintSense.Range


def test_constraint_rows():
    model = mom.read_mps_bytes(RANGES_MPS.encode())
    begin, end = model.constraints_starts[0], model.constraint_end(0)
    assert list(model.elements_variables[begin:end]) == [0, 1]
    assert list(model.elements_coefficients[begin:end]) == [1.0, 3.0]
    np.testing.assert_array_equal(model.objective_coefficients, [1, 2])


def test_read_lp_variable_types_and_classes():
    model = mom.read_lp_bytes(LP.encode())
    assert model.objective_direction == mom.ObjectiveDirection.Maximize
    types = dict(zip(model.variables_names, model.variables_types))
    assert types == {"x": mom.VariableType.Binary, "y": mom.VariableType.Continuous,
                     "z": mom.VariableType.Integer}
    assert len(model.constraints_classes) == 2
    assert model.constraints_classes[1] == model.constraint_class(1)


def test_arrays_are_read_only_views_that_outlive_the_model():
    model = mom.read_mps_bytes(RANGES_MPS.encode())
    lower = model.variables_lower_bounds
    assert not lower.flags.writeable
    with pytest.raises(ValueError):
        lower[0] = 42
    del model
    gc.collect()
    np.testing.assert_array_equal(lower, [0, -np.inf])


def test_parse_error_raises_value_error():
    with pytest.raises(ValueError):
        mom.read_nl_bytes(b"not an nl file\n")


@pytest.mark.parametrize("opener,suffix", [(open, ""), (gzip.open, ".gz"), (bz2.open, ".bz2")])
def test_read_handles_compression(tmp_path, opener, suffix):
    path = tmp_path / ("model.lp" + suffix)
    with opener(path, "wt") as f:
        f.write(LP)
    model = mom.read(path)
    assert model.number_of_variables == 3


def test_detect_format():
    assert mom.detect_format("a/b/model.MPS") == "mps"
    assert mom.detect_format("model.nl.bz2") == "nl"
    assert mom.detect_format("model.txt") is None
    with pytest.raises(ValueError):
        mom.read("model.txt")


@pytest.mark.parametrize("format", ["mps", "lp", "json", "nl"])
def test_write_then_read_round_trip(tmp_path, format):
    model = mom.read_lp_bytes(LP.encode())
    path = tmp_path / f"model.{format}"
    mom.write(model, path)
    again = mom.read(path)
    assert again.number_of_variables == model.number_of_variables
    assert again.number_of_constraints == model.number_of_constraints
    assert again.number_of_elements == model.number_of_elements
    assert again.objective_direction == model.objective_direction
    assert sorted(again.variables_types) == sorted(model.variables_types)


@pytest.mark.parametrize("format", ["mps", "lp", "json", "nl"])
@pytest.mark.parametrize("compression", ["", ".gz", ".bz2"])
def test_write_file_matches_string(tmp_path, format, compression):
    # Uncompressed files are written by the C++ writers directly, compressed
    # ones through the string writers: both must give the same contents.
    model = mom.read_lp_bytes(LP.encode())
    path = tmp_path / f"model.{format}{compression}"
    mom.write(model, path)
    opener = {"": open, ".gz": gzip.open, ".bz2": bz2.open}[compression]
    with opener(path, "rt") as f:
        assert f.read() == getattr(mom, f"write_{format}_string")(model)


def test_nonlinear_expression_tree():
    # minimize exp(x0) + x0 * x1, written as an AMPL .nl text file.
    nl = "\n".join([
        "g0", "2 0 1 0 0", "0 1", "0 0", "0 2 0", "0 0", "0 0 0 0 0", "0 0", "0 0",
        "0 0 0 0 0",
        "b", "3", "3",
        "O0 0", "o0", "o44", "v0", "o2", "v0", "v1",
    ]) + "\n"
    model = mom.read_nl_bytes(nl.encode())
    assert model.has_nonlinear()
    assert model.objective_nonlinear_elements_operators == b"+ev*vv"
    assert list(model.objective_nonlinear_elements_number_of_children) == [2, 1, 0, 2, 0, 0]
    assert model.evaluate_objective([1.0, 2.0]) == pytest.approx(np.e + 2.0)


@pytest.mark.skipif(not os.path.exists(P0033), reason="repository data not available")
def test_read_p0033():
    model = mom.read(P0033)
    assert (model.number_of_variables, model.number_of_constraints) == (33, 16)
    assert all(t == mom.VariableType.Binary for t in model.variables_types)
    assert model.is_milp() and not model.is_lp()


COMMON_SUBEXPRESSIONS_NL = "\n".join([
    # V2 = 2 x0 + x0^2, V3 = V2 * x1; minimize V2 s.t. c0: V3 <= 100.
    "g0", "2 1 1 0 0", "1 1", "0 0", "2 1 1", "0 0", "0 0 0 0 0", "0 0", "0 0", "2 0 0 0 0",
    "b", "3", "3", "r", "1 100",
    "V2 1 0", "0 2", "o5", "v0", "n2",
    "V3 0 0", "o2", "v2", "v1",
    "C0", "v3", "O0 0", "v2",
]) + "\n"


def test_common_subexpressions_are_defined_variables():
    model = mom.read_nl_bytes(COMMON_SUBEXPRESSIONS_NL.encode())
    assert model.has_defined_variables() and model.check()
    assert model.number_of_defined_variables == 2
    assert model.nonlinear_elements_operators == b"d"
    x = [3.0, 2.0]
    assert model.evaluate_objective(x) == pytest.approx(15.0)
    assert model.evaluate_constraint(x, 0) == pytest.approx(30.0)
    expanded = mom.expand_defined_variables(model)
    assert not expanded.has_defined_variables()
    assert expanded.evaluate_constraint(x, 0) == pytest.approx(30.0)


def test_abs_tanh_log10_min_max():
    # minimize abs(x0) + tanh(x1) + log10(x1) + min(x0, x1) + max(x0, x1)
    nl = "\n".join([
        "g0", "2 0 1 0 0", "0 1", "0 0", "0 2 0", "0 0", "0 0 0 0 0", "0 0", "0 0", "0 0 0 0 0",
        "b", "3", "3", "O0 0", "o54", "5", "o15", "v0", "o37", "v1", "o42", "v1",
        "o11", "2", "v0", "v1", "o12", "2", "v0", "v1",
    ]) + "\n"
    model = mom.read_nl_bytes(nl.encode())
    x0, x1 = -2.0, 5.0
    expected = abs(x0) + np.tanh(x1) + np.log10(x1) + min(x0, x1) + max(x0, x1)
    assert model.evaluate_objective([x0, x1]) == pytest.approx(expected)


def test_quadratic_mps_sections():
    mps = "\n".join([
        "NAME QP", "ROWS", " N obj", " L c1", "COLUMNS",
        "    x obj 1 c1 1", "    y obj 1 c1 1", "RHS", "    RHS c1 10",
        "BOUNDS", " FR BND x", " FR BND y",
        # objective += 1/2 (2 x^2 + 2 x y + 4 y^2); c1 += x^2 + x y + 3 y^2.
        "QUADOBJ", "    x x 2", "    x y 1", "    y y 4",
        "QCMATRIX c1", "    x x 1", "    x y 0.5", "    y x 0.5", "    y y 3",
        "ENDATA",
    ]) + "\n"
    model = mom.read_mps_bytes(mps.encode())
    assert model.has_quadratic()
    x, y = 1.5, -2.0
    assert model.evaluate_objective([x, y]) == pytest.approx(x + y + x * x + x * y + 2 * y * y)
    assert model.evaluate_constraint([x, y], 0) == pytest.approx(x + y + x * x + x * y + 3 * y * y)
    assert model.constraint_class(0) == mom.ConstraintClass.GeneralQuadratic


@pytest.fixture
def comma_locale():
    """Set the C locale to one with a decimal comma (as the JVM, or a program
    calling setlocale(LC_ALL, ""), does with such an environment)."""
    import locale
    saved = locale.setlocale(locale.LC_ALL)
    for name in ("fr_FR.UTF-8", "fr_FR.utf8", "de_DE.UTF-8", "de_DE.utf8", "French_France.1252"):
        try:
            locale.setlocale(locale.LC_ALL, name)
        except locale.Error:
            continue
        if locale.localeconv()["decimal_point"] == ",":
            break
    else:
        locale.setlocale(locale.LC_ALL, saved)
        pytest.skip("no locale with a decimal comma installed")
    yield
    locale.setlocale(locale.LC_ALL, saved)


@pytest.mark.parametrize("format", ["mps", "lp", "json", "nl"])
def test_readers_ignore_locale(comma_locale, format):
    model = mom.Model()
    model.add_variables(lower_bounds=0.0, upper_bounds=[2.5, 4.125], objective_coefficients=[1.5, 2.25])
    model.objective_direction = mom.ObjectiveDirection.Maximize
    model.add_constraints(upper_bounds=3.5, starts=np.array([0], dtype=np.int32),
                          variables=np.array([0, 1], dtype=np.int32), coefficients=[0.5, 1.75])
    text = getattr(mom, f"write_{format}_string")(model)
    read_back = getattr(mom, f"read_{format}_bytes")(text.encode())
    assert list(read_back.objective_coefficients) == [1.5, 2.25]
    assert list(read_back.variables_upper_bounds) == [2.5, 4.125]
    assert list(read_back.constraints_upper_bounds) == [3.5]
    assert list(read_back.elements_coefficients) == [0.5, 1.75]


def test_read_nl_names(tmp_path):
    model = mom.Model()
    model.add_variables(2, lower_bounds=0.0, upper_bounds=4.0, objective_coefficients=[1.0, 2.0])
    model.add_constraints(upper_bounds=3.0, starts=np.array([0], dtype=np.int32),
                          variables=np.array([0, 1], dtype=np.int32), coefficients=[1.0, 1.0])
    mom.write(model, tmp_path / "model.nl")
    assert mom.read(tmp_path / "model.nl").variables_names == []

    (tmp_path / "model.col").write_text("x[1]\ny[1]\n")
    (tmp_path / "model.row").write_text("capacity\nobj\n")
    read_back = mom.read(tmp_path / "model.nl")
    assert read_back.variables_names == ["x[1]", "y[1]"]
    assert read_back.constraints_names == ["capacity"]

    # Compressed .nl: the names files are next to it, uncompressed.
    with open(tmp_path / "model.nl", "rb") as f, gzip.open(tmp_path / "other.nl.gz", "wb") as g:
        g.write(f.read())
    (tmp_path / "other.col").write_text("a\nb\n")
    assert mom.read(tmp_path / "other.nl.gz").variables_names == ["a", "b"]

    # A names file that doesn't match the model is ignored.
    (tmp_path / "model.col").write_text("x\ny\nz\n")
    assert mom.read(tmp_path / "model.nl").variables_names == []
    assert mom.read(tmp_path / "model.nl").constraints_names == ["capacity"]

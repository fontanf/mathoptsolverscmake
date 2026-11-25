import math

import numpy as np
import pytest

import mathoptmodel as mom

INF = math.inf


def small_milp():
    """max x + 2 y + 3 z  s.t.  x + y + z <= 10,  x - y >= -3,  x in [0, 10],
    y binary, z integer in [0, 5]."""
    model = mom.Model()
    first = model.add_variables(
        lower_bounds=0.0,
        upper_bounds=np.array([10.0, 1.0, 5.0]),
        types=np.array([mom.VariableType.Continuous, mom.VariableType.Binary,
                        mom.VariableType.Integer], dtype=np.int8),
        names=["x", "y", "z"])
    assert first == 0
    model.set_objective(mom.ObjectiveDirection.Maximize, coefficients=np.array([1.0, 2.0, 3.0]))
    first = model.add_constraints(
        lower_bounds=np.array([-INF, -3.0]),
        upper_bounds=np.array([10.0, INF]),
        starts=np.array([0, 3], dtype=np.int32),
        variables=np.array([0, 1, 2, 0, 1], dtype=np.int32),
        coefficients=np.array([1.0, 1.0, 1.0, 1.0, -1.0]),
        names=["c1", "c2"])
    assert first == 0
    return model


def test_empty_model():
    model = mom.Model()
    assert model.number_of_variables == 0
    assert model.number_of_constraints == 0
    assert model.objective_direction == mom.ObjectiveDirection.Minimize
    assert model.check()


def test_build_milp_matches_lp_reader():
    built = small_milp()
    assert built.check()
    read = mom.read_lp_bytes(mom.write_lp_string(built).encode())
    assert read.objective_direction == mom.ObjectiveDirection.Maximize
    assert read.variables_names == ["x", "y", "z"]
    assert read.constraints_names == ["c1", "c2"]
    for field in ("variables_lower_bounds", "variables_upper_bounds", "variables_types",
                  "objective_coefficients", "constraints_lower_bounds",
                  "constraints_upper_bounds", "constraints_starts",
                  "elements_variables", "elements_coefficients"):
        np.testing.assert_array_equal(getattr(read, field), getattr(built, field), err_msg=field)
    assert built.evaluate_objective([1.0, 1.0, 2.0]) == 9.0


def test_defaults_and_scalars():
    model = mom.Model()
    assert model.add_variables(3) == 0
    np.testing.assert_array_equal(model.variables_lower_bounds, [-INF] * 3)
    np.testing.assert_array_equal(model.variables_upper_bounds, [INF] * 3)
    np.testing.assert_array_equal(model.variables_types, [0] * 3)
    np.testing.assert_array_equal(model.objective_coefficients, [0.0] * 3)
    assert model.add_variables(2, lower_bounds=0, upper_bounds=1,
                               types=mom.VariableType.Binary, objective_coefficients=4) == 3
    np.testing.assert_array_equal(model.variables_lower_bounds[3:], [0.0, 0.0])
    np.testing.assert_array_equal(model.variables_types[3:], [1, 1])
    np.testing.assert_array_equal(model.objective_coefficients[3:], [4.0, 4.0])
    # Constraints with no terms.
    assert model.add_constraints(2, upper_bounds=1.0) == 0
    np.testing.assert_array_equal(model.constraints_starts, [0, 0])
    np.testing.assert_array_equal(model.constraints_lower_bounds, [-INF, -INF])
    assert model.check()


def test_lists_and_other_dtypes_are_converted():
    model = mom.Model()
    model.add_variables(lower_bounds=[0, 0], upper_bounds=(1, 2), types=[0, 2])
    model.add_constraints(upper_bounds=np.array([3], dtype=np.int64),
                          starts=np.array([0], dtype=np.int64),
                          variables=np.array([0, 1], dtype=np.int64),
                          coefficients=np.array([1, 1], dtype=np.float32))
    # A non-contiguous array.
    model.set_objective(coefficients=np.arange(4.0)[::2])
    np.testing.assert_array_equal(model.objective_coefficients, [0.0, 2.0])
    np.testing.assert_array_equal(model.variables_types, [0, 2])
    assert model.check()


def test_batches_are_appended():
    model = mom.Model()
    model.add_variables(2, lower_bounds=0.0)
    model.add_constraints(upper_bounds=[1.0], starts=[0], variables=[0, 1], coefficients=[1.0, 2.0])
    assert model.add_variables(1, lower_bounds=0.0) == 2
    assert model.add_constraints(upper_bounds=[5.0, 6.0], starts=[0, 1],
                                 variables=[2, 0, 2], coefficients=[3.0, 4.0, 5.0]) == 1
    np.testing.assert_array_equal(model.constraints_starts, [0, 2, 3])
    np.testing.assert_array_equal(model.elements_variables, [0, 1, 2, 0, 2])
    assert model.constraint_end(1) == 3
    assert model.check()
    solution = [1.0, 1.0, 1.0]
    assert [model.evaluate_constraint(solution, i) for i in range(3)] == [3.0, 3.0, 9.0]


def test_names_are_completed():
    model = mom.Model()
    model.add_variables(2)
    model.add_variables(1, names=["z"])
    model.add_variables(1)
    assert model.variables_names == ["x0", "x1", "z", "x3"]
    model.add_constraints(1)
    model.add_constraints(1, names=["capacity"])
    assert model.constraints_names == ["c0", "capacity"]
    assert model.check()


def test_initial_values():
    model = mom.Model()
    model.add_variables(2)
    model.add_variables(initial_values=[1.5, 2.5])
    model.add_variables(1)
    np.testing.assert_array_equal(model.variables_initial_values, [0.0, 0.0, 1.5, 2.5, 0.0])
    assert model.check()


def test_quadratic():
    model = mom.Model()
    model.add_variables(2, lower_bounds=0.0)
    # min x^2 + x y  s.t.  [] , x^2 + y^2 <= 4 ; x + y >= 1 (no quadratic terms).
    model.set_objective(quadratic_variables_1=[0, 0], quadratic_variables_2=[0, 1],
                        quadratic_coefficients=[1.0, 1.0])
    model.add_constraints(upper_bounds=[4.0, INF], lower_bounds=[-INF, 1.0],
                          starts=[0, 0], variables=[0, 1], coefficients=[1.0, 1.0],
                          quadratic_starts=[0, 2], quadratic_variables_1=[0, 1],
                          quadratic_variables_2=[0, 1], quadratic_coefficients=[1.0, 1.0])
    assert model.has_quadratic()
    assert model.check()
    np.testing.assert_array_equal(model.quadratic_elements_constraints_starts, [0, 2])
    assert model.evaluate_objective([1.0, 2.0]) == 3.0
    assert model.evaluate_constraint([1.0, 2.0], 0) == 5.0
    assert model.evaluate_constraint([1.0, 2.0], 1) == 3.0
    # A later batch without quadratic terms keeps the starts consistent.
    model.add_constraints(1)
    np.testing.assert_array_equal(model.quadratic_elements_constraints_starts, [0, 2, 2])
    assert model.check()
    read = mom.read_json_bytes(mom.write_json_string(model).encode())
    assert read.evaluate_constraint([1.0, 2.0], 0) == 5.0


def test_nonlinear_and_defined_variables():
    model = mom.Model()
    model.add_variables(2, lower_bounds=0.0)
    model.add_constraints(1, upper_bounds=10.0, starts=[0], variables=[0], coefficients=[1.0])
    # c1: exp(x0) + x1 ^ 2 <= 4 ; c2: (no tree) ; c3: log(x0)
    model.add_constraints(
        upper_bounds=[4.0, 1.0, 0.0],
        nonlinear_starts=[0, 6, 6],
        nonlinear_operators=b"+evpvk" + b"lv",
        nonlinear_number_of_children=[2, 1, 0, 2, 0, 0, 1, 0],
        nonlinear_variables=[-1, -1, 0, -1, 1, -1, -1, 0],
        nonlinear_values=[0, 0, 0, 0, 0, 2, 0, 0])
    assert model.has_nonlinear()
    np.testing.assert_array_equal(model.nonlinear_elements_constraints_starts, [0, 0, 6, 6])
    assert model.check()
    solution = [1.0, 2.0]
    assert model.evaluate_constraint(solution, 1) == pytest.approx(math.e + 4.0)
    assert model.evaluate_constraint(solution, 2) == 0.0
    assert model.evaluate_constraint(solution, 3) == 0.0
    assert model.constraint_class(1) == mom.ConstraintClass.GeneralNonlinear

    # Defined variables: d0 = x0 * x1, d1 = d0 + x0 (operators as a uint8 array).
    assert model.add_defined_variables(
        starts=[0, 3],
        operators=np.frombuffer(b"*vv+dv", dtype=np.uint8),
        number_of_children=[2, 0, 0, 2, 0, 0],
        variables=[-1, 0, 1, -1, 0, 0]) == 0
    model.set_objective(nonlinear_operators=b"d", nonlinear_number_of_children=[0],
                        nonlinear_variables=[1])
    assert model.check()
    assert model.evaluate_objective([2.0, 3.0]) == 8.0
    expanded = mom.expand_defined_variables(model)
    assert not expanded.has_defined_variables()
    assert expanded.evaluate_objective([2.0, 3.0]) == 8.0
    # The .nl round trip keeps the values (.nl stores the nonlinear constraints
    # first, hence the sort).
    read = mom.read_nl_bytes(mom.write_nl_string(model).encode())
    assert read.evaluate_objective([2.0, 3.0]) == 8.0

    def values(m):
        return sorted(m.evaluate_constraint(solution, i) for i in range(m.number_of_constraints))

    assert values(read) == pytest.approx(values(model))

    # Removing the objective's tree.
    model.set_objective(nonlinear_operators=b"")
    assert model.evaluate_objective([2.0, 3.0]) == 0.0


def test_set_objective_keeps_unspecified_parts():
    model = small_milp()
    model.set_objective(direction=mom.ObjectiveDirection.Minimize)
    np.testing.assert_array_equal(model.objective_coefficients, [1.0, 2.0, 3.0])
    model.set_objective(coefficients=0.0)
    np.testing.assert_array_equal(model.objective_coefficients, [0.0, 0.0, 0.0])
    model.objective_direction = mom.ObjectiveDirection.Maximize
    assert model.objective_direction == mom.ObjectiveDirection.Maximize


@pytest.mark.parametrize("kwargs, message", [
    (dict(starts=[0], variables=[3], coefficients=[1.0]), "not a variable index"),
    (dict(starts=[0], variables=[0]), "given together"),
    (dict(variables=[0], coefficients=[1.0], number=1), "starts is required"),
    (dict(starts=[1], variables=[0], coefficients=[1.0]), "starts must begin at 0"),
    (dict(starts=[0, 2], variables=[0], coefficients=[1.0]), "starts must begin at 0"),
    (dict(number=2, upper_bounds=[1.0]), "were expected"),
    (dict(), "number of constraints is unknown"),
    (dict(quadratic_starts=[0], quadratic_variables_1=[0], quadratic_variables_2=[0]),
     "given together"),
    (dict(nonlinear_starts=[0], nonlinear_operators=b"+v", nonlinear_number_of_children=[2, 0],
          nonlinear_variables=[0, 0]), "missing 1 node"),
    (dict(nonlinear_starts=[0], nonlinear_operators=b"vv", nonlinear_number_of_children=[0, 0],
          nonlinear_variables=[0, 0]), "single root"),
    (dict(nonlinear_starts=[0], nonlinear_operators=b"ev", nonlinear_number_of_children=[2, 0],
          nonlinear_variables=[0, 0]), "has 2 children"),
    (dict(nonlinear_starts=[0], nonlinear_operators=b"Z", nonlinear_number_of_children=[0]),
     "unknown operator"),
    (dict(nonlinear_starts=[0], nonlinear_operators=b"v", nonlinear_number_of_children=[0]),
     "references variable -1"),
    (dict(nonlinear_starts=[0], nonlinear_operators=b"d", nonlinear_number_of_children=[0],
          nonlinear_variables=[0]), "defined variable 0"),
    (dict(nonlinear_starts=[0], nonlinear_operators=b"v", nonlinear_variables=[0]),
     "number_of_children is required"),
])
def test_invalid_constraints_leave_the_model_unchanged(kwargs, message):
    model = small_milp()
    before = mom.write_json_string(model)
    with pytest.raises(ValueError, match=message):
        model.add_constraints(**kwargs)
    assert mom.write_json_string(model) == before
    assert model.check()


def test_invalid_variables_and_objective():
    model = small_milp()
    with pytest.raises(ValueError, match="not a VariableType"):
        model.add_variables(types=[5])
    with pytest.raises(ValueError, match="number of variables is unknown"):
        model.add_variables()
    with pytest.raises(ValueError, match="3 variables"):
        model.set_objective(coefficients=[1.0])
    with pytest.raises(ValueError, match="not a variable index"):
        model.set_objective(quadratic_variables_1=[0], quadratic_variables_2=[7],
                            quadratic_coefficients=[1.0])
    with pytest.raises(ValueError, match="only reference lower indices"):
        model.add_defined_variables(starts=[0], operators=b"d", number_of_children=[0],
                                    variables=[0])
    with pytest.raises(ValueError, match="empty expression tree"):
        model.add_defined_variables(starts=[0, 1], operators=b"v", number_of_children=[0],
                                    variables=[0])
    assert model.number_of_variables == 3
    assert model.number_of_defined_variables == 0
    assert model.check()


def test_large_batch():
    rng = np.random.default_rng(0)
    n, m, k = 20_000, 10_000, 4
    model = mom.Model()
    model.add_variables(n, lower_bounds=0.0, upper_bounds=1.0, types=mom.VariableType.Binary)
    variables = np.concatenate([rng.choice(n, size=k, replace=False) for _ in range(m)])
    model.add_constraints(upper_bounds=1.0, starts=np.arange(0, m * k, k),
                          variables=variables, coefficients=np.ones(m * k))
    assert model.number_of_elements == m * k
    assert model.check()
    assert (model.constraints_classes == int(mom.ConstraintClass.SetPacking)).all()


def test_semi_variables():
    model = mom.Model()
    model.add_variables(lower_bounds=[1.0, 2.0], upper_bounds=[10.0, 8.0],
                        types=np.array([mom.VariableType.SemiContinuous, mom.VariableType.SemiInteger],
                                       dtype=np.int8))
    model.add_constraints(lower_bounds=0.5, starts=np.array([0], dtype=np.int32),
                          variables=np.array([0, 1], dtype=np.int32), coefficients=[1.0, 1.0])
    assert list(model.variables_types) == [int(mom.VariableType.SemiContinuous),
                                           int(mom.VariableType.SemiInteger)]
    assert model.check_solution([0.0, 0.0]) is False  # violates the constraint
    assert model.check_solution([0.0, 3.0]) is True
    assert model.check_solution([0.0, 3.5]) is False  # fractional semi-integer
    assert model.check_solution([0.6, 0.0]) is False  # neither 0 nor in [1, 10]
    for format in ("mps", "lp", "json"):
        text = getattr(mom, f"write_{format}_string")(model)
        read_back = getattr(mom, f"read_{format}_bytes")(text.encode())
        assert list(read_back.variables_types) == list(model.variables_types), format
    with pytest.raises(ValueError, match="semi"):
        mom.write_nl_string(model)


def test_format_expressions():
    model = mom.Model()
    model.add_variables(2, lower_bounds=0.0, upper_bounds=4.0, objective_coefficients=[1.0, -2.0])
    # 2 x0 + exp(x1) <= 3
    model.add_constraints(upper_bounds=3.0, starts=np.array([0], dtype=np.int32),
                          variables=np.array([0], dtype=np.int32), coefficients=[2.0],
                          nonlinear_starts=np.array([0], dtype=np.int32), nonlinear_operators=b"ev",
                          nonlinear_variables=np.array([-1, 1], dtype=np.int32),
                          nonlinear_number_of_children=np.array([1, 0], dtype=np.int32))
    assert model.format_constraint_expression(0) == "2 * x0 + exp(x1)"
    assert model.format_constraint(0) == "c0: 2 * x0 + exp(x1) <= 3\n"
    assert model.format_objective_expression() == "x0 - 2 * x1"
    assert model.constraint_degree(0) == math.inf
    assert model.objective_degree() == 1

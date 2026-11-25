import math
import os

import numpy as np
import pytest

import mathoptmodel as mom

INF = math.inf
P0033 = os.path.join(os.path.dirname(__file__), "..", "..", "data", "miplib3", "p0033.mps")

needs_highs = pytest.mark.skipif(
    mom.Solver.Highs not in mom.available_solvers(), reason="HiGHS not available in this build")


def small_milp():
    """max x + 2 y + 3 z  s.t.  x + y + z <= 10,  x - y >= -3,  x in [0, 10],
    y binary, z integer in [0, 5]. Optimal: x = 4, y = 1, z = 5, objective 21."""
    model = mom.Model()
    model.add_variables(
        lower_bounds=0.0,
        upper_bounds=[10.0, 1.0, 5.0],
        types=np.array([mom.VariableType.Continuous, mom.VariableType.Binary,
                        mom.VariableType.Integer], dtype=np.int8))
    model.set_objective(mom.ObjectiveDirection.Maximize, coefficients=[1.0, 2.0, 3.0])
    model.add_constraints(
        lower_bounds=[-INF, -3.0],
        upper_bounds=[10.0, INF],
        starts=np.array([0, 3], dtype=np.int32),
        variables=np.array([0, 1, 2, 0, 1], dtype=np.int32),
        coefficients=[1.0, 1.0, 1.0, 1.0, -1.0])
    return model


def infeasible_lp():
    """x >= 1 and x <= 0."""
    model = mom.Model()
    model.add_variables(1, lower_bounds=-INF, upper_bounds=INF)
    model.add_constraints(
        lower_bounds=[1.0, -INF], upper_bounds=[INF, 0.0],
        starts=np.array([0, 1], dtype=np.int32),
        variables=np.array([0, 0], dtype=np.int32), coefficients=[1.0, 1.0])
    return model


def test_available_solvers():
    solvers = mom.available_solvers()
    assert all(isinstance(solver, mom.Solver) for solver in solvers)
    assert len(set(solvers)) == len(solvers)


@needs_highs
def test_solve_milp():
    result = mom.solve(small_milp(), "highs", verbose=False)
    assert isinstance(result, mom.SolveResult)
    assert result.found_solution
    assert result.objective_value == pytest.approx(21)
    assert result.solution == pytest.approx([4, 1, 5])
    assert result.solution.dtype == np.float64
    assert result.has_bound
    assert result.bound == pytest.approx(21)
    assert "objective value 21" in repr(result)


@needs_highs
@pytest.mark.parametrize("solver", ["highs", "HiGHS", "HIGHS", mom.Solver.Highs])
def test_solver_names(solver):
    assert mom.solve(small_milp(), solver, verbose=False).objective_value == pytest.approx(21)


@needs_highs
def test_solve_infeasible():
    result = mom.solve(infeasible_lp(), "highs", verbose=False)
    assert not result.found_solution
    assert len(result.solution) == 0
    assert math.isnan(result.objective_value)
    assert "no solution found" in repr(result)


@needs_highs
def test_time_limit():
    result = mom.solve(small_milp(), "highs", time_limit=10, verbose=False)
    assert result.objective_value == pytest.approx(21)
    with pytest.raises(ValueError, match="time limit"):
        mom.solve(small_milp(), "highs", time_limit=-1, verbose=False)


@needs_highs
@pytest.mark.skipif(not os.path.exists(P0033), reason="repository data not available")
def test_solve_file():
    result = mom.solve(mom.read(P0033), "highs", verbose=False)
    assert result.objective_value == pytest.approx(3089)
    assert result.bound == pytest.approx(3089)


@needs_highs
def test_nonlinear_model_rejected():
    model = small_milp()
    model.set_objective(quadratic_variables_1=np.array([0], dtype=np.int32),
                        quadratic_variables_2=np.array([0], dtype=np.int32),
                        quadratic_coefficients=[1.0])
    with pytest.raises(ValueError, match="MILP"):
        mom.solve(model, "highs", verbose=False)


def test_unknown_solver():
    with pytest.raises(ValueError, match="Unknown solver"):
        mom.solve(small_milp(), "nosolver")


def test_unavailable_solver():
    unavailable = [solver for solver in mom.Solver if solver not in mom.available_solvers()]
    if not unavailable:
        pytest.skip("every solver is available in this build")
    with pytest.raises(ValueError, match="not available in this build"):
        mom.solve(small_milp(), unavailable[0])

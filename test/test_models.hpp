#pragma once

#include "mathoptsolverscmake/mathopt.hpp"

#include <ostream>
#include <string>
#include <vector>

namespace mathoptsolverscmaketest
{

/**
 * A MathOptModel bundled with its known, unique optimal solution, so a test
 * has a single source of truth for what a model builder is expected to
 * produce (no separate table to keep in sync). The optimal objective value
 * and the solution's feasibility can both be derived from the model itself
 * via MathOptModel::evaluate_objective() / check_solution().
 *
 * For an infeasible model, "feasible" is false and "optimal_solution" is
 * empty: there is nothing to evaluate, the model itself is the fixture and
 * the test just checks that the solver reports infeasibility.
 */
struct TestModel
{
    std::string name;
    mathoptsolverscmake::MathOptModel model;
    std::vector<double> optimal_solution;
    bool feasible = true;
};

std::ostream& operator<<(std::ostream& os, const TestModel& test_model);

/*
 * Each function below builds a small TestModel, exercising one structure
 * type (linear, quadratic, nonlinear expression tree, black-box) in a
 * continuous and a mixed-integer variant, plus one infeasible variant of
 * each (a continuous-relaxation-infeasible one, and one that is only
 * infeasible because of integrality).
 */

TestModel linear_continuous_model();
TestModel linear_mixed_integer_model();
TestModel linear_infeasible_model();
TestModel linear_mixed_integer_infeasible_model();

TestModel quadratic_continuous_model();
TestModel quadratic_mixed_integer_model();
TestModel quadratic_infeasible_model();
TestModel quadratic_mixed_integer_infeasible_model();

TestModel nonlinear_continuous_model();
TestModel nonlinear_mixed_integer_model();
TestModel nonlinear_infeasible_model();
TestModel nonlinear_mixed_integer_infeasible_model();

TestModel black_box_continuous_model();
TestModel black_box_mixed_integer_model();
TestModel black_box_infeasible_model();
TestModel black_box_mixed_integer_infeasible_model();

}

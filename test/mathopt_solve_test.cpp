#include "test_models.hpp"

#include "mathoptsolverscmake/mathopt_solve.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <sstream>
#include <stdexcept>

using namespace mathoptsolverscmake;
using namespace mathoptsolverscmaketest;

namespace
{

std::string name(SolverName solver)
{
    std::ostringstream os;
    os << solver;
    return os.str();
}

/** Solvers of models with integer variables (Ipopt only solves continuous ones). */
bool solves_integer_models(SolverName solver)
{
    return solver != SolverName::Ipopt;
}

void expect_optimal(
        const TestModel& test_model,
        SolverName solver)
{
    SolveParameters parameters;
    parameters.solver = solver;
    parameters.verbose = false;
    SolveOutput output = solve(test_model.model, parameters);
    ASSERT_TRUE(output.found_solution) << name(solver) << " " << test_model;
    ASSERT_EQ(output.solution.size(), (size_t)test_model.model.number_of_variables());
    double expected = test_model.model.evaluate_objective(test_model.optimal_solution);
    EXPECT_NEAR(output.objective_value, expected, 1e-5) << name(solver) << " " << test_model;
    EXPECT_NEAR(
            output.objective_value,
            test_model.model.evaluate_objective(output.solution),
            1e-5) << name(solver);
    if (output.has_bound)
        EXPECT_NEAR(output.bound, expected, 1e-4) << name(solver) << " " << test_model;
}

}

TEST(MathOptSolve, AvailableSolvers)
{
    for (SolverName solver: available_solvers())
        EXPECT_TRUE(is_available(solver)) << name(solver);
    EXPECT_FALSE(is_available(SolverName::Dlib));
}

TEST(MathOptSolve, LinearContinuous)
{
    for (SolverName solver: available_solvers())
        expect_optimal(linear_continuous_model(), solver);
}

TEST(MathOptSolve, LinearMixedInteger)
{
    for (SolverName solver: available_solvers())
        if (solves_integer_models(solver))
            expect_optimal(linear_mixed_integer_model(), solver);
}

TEST(MathOptSolve, Infeasible)
{
    for (SolverName solver: available_solvers()) {
        SolveParameters parameters;
        parameters.solver = solver;
        parameters.verbose = false;
        SolveOutput output = solve(linear_infeasible_model().model, parameters);
        EXPECT_FALSE(output.found_solution) << name(solver);
        EXPECT_TRUE(output.solution.empty()) << name(solver);
        EXPECT_TRUE(std::isnan(output.objective_value)) << name(solver);
    }
}

TEST(MathOptSolve, TimeLimit)
{
    for (SolverName solver: available_solvers()) {
        SolveParameters parameters;
        parameters.solver = solver;
        parameters.verbose = false;
        parameters.time_limit = 10;
        SolveOutput output = solve(linear_continuous_model().model, parameters);
        EXPECT_TRUE(output.found_solution) << name(solver);
    }
}

TEST(MathOptSolve, InvalidParameters)
{
    SolveParameters parameters;
    parameters.solver = SolverName::Dlib;
    EXPECT_THROW(solve(linear_continuous_model().model, parameters), std::invalid_argument);

    for (SolverName solver: {
            SolverName::Cbc, SolverName::Highs, SolverName::Xpress, SolverName::Cplex,
            SolverName::Gurobi, SolverName::Knitro, SolverName::Ipopt, SolverName::Scip}) {
        parameters.solver = solver;
        parameters.time_limit = -1;
        EXPECT_THROW(solve(linear_continuous_model().model, parameters), std::invalid_argument)
            << name(solver);
        parameters.time_limit = std::numeric_limits<double>::infinity();
        if (!is_available(solver)) {
            EXPECT_THROW(solve(linear_continuous_model().model, parameters), std::invalid_argument)
                << name(solver);
        }
    }
}

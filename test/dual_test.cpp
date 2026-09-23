#include "mathoptsolverscmake/dual.hpp"
#include "mathoptsolverscmake/mathopt_mps.hpp"
#if HIGHS_FOUND
#include "mathoptsolverscmake/mathopt_highs.hpp"
#endif

#include <gtest/gtest.h>

#include <limits>
#include <stdexcept>

using namespace mathoptsolverscmake;

namespace
{

const double inf = std::numeric_limits<double>::infinity();

void add_constraint(
        MathOptModel& model,
        const std::vector<int>& variables,
        const std::vector<double>& coefficients,
        double lower_bound,
        double upper_bound)
{
    model.constraints_starts.push_back(model.number_of_elements());
    model.elements_variables.insert(model.elements_variables.end(), variables.begin(), variables.end());
    model.elements_coefficients.insert(model.elements_coefficients.end(), coefficients.begin(), coefficients.end());
    model.constraints_lower_bounds.push_back(lower_bound);
    model.constraints_upper_bounds.push_back(upper_bound);
    model.quadratic_elements_constraints_starts.push_back(0);
    model.constraints_functions.push_back(BlackBoxFunction());
}

/**
 * An LP exercising every kind of constraint side and variable bound.
 *
 * min 2 x0 + 3 x1 - x2 + x3 + x4
 * s.t.      x0 + x1 + x2 + x4 >= 2
 *           x0 - x1           <= 3
 *      1 <= x1 + x2 + x3      <= 5
 *           x0 + x3            = 1.5
 *           x2 - x3              free
 *      0 <= x0, -1 <= x1 <= 4, x2 <= 2, x3 free, x4 = 2.5
 */
MathOptModel general_lp(ObjectiveDirection objective_direction)
{
    MathOptModel model(5);
    model.objective_direction = objective_direction;
    model.objective_coefficients = {2, 3, -1, 1, 1};
    if (objective_direction == ObjectiveDirection::Maximize)
        for (double& coefficient: model.objective_coefficients)
            coefficient = -coefficient;
    model.variables_lower_bounds = {0, -1, -inf, -inf, 2.5};
    model.variables_upper_bounds = {inf, 4, 2, inf, 2.5};
    add_constraint(model, {0, 1, 2, 4}, {1, 1, 1, 1}, 2, inf);
    add_constraint(model, {0, 1}, {1, -1}, -inf, 3);
    add_constraint(model, {1, 2, 3}, {1, 1, 1}, 1, 5);
    add_constraint(model, {0, 3}, {1, 1}, 1.5, 1.5);
    add_constraint(model, {2, 3}, {1, -1}, -inf, inf);
    return model;
}

}

TEST(MathOptDual, TextbookForm)
{
    // min c^T x s.t. A x >= b, x >= 0  =>  max b^T y s.t. A^T y <= c, y >= 0
    MathOptModel model(2);
    model.objective_direction = ObjectiveDirection::Minimize;
    model.objective_coefficients = {3, 5};
    model.variables_lower_bounds = {0, 0};
    add_constraint(model, {0, 1}, {1, 2}, 4, inf);
    add_constraint(model, {0}, {1}, 1, inf);

    MathOptModel dual_model = dual(model);
    ASSERT_TRUE(dual_model.check(1));
    EXPECT_EQ(dual_model.objective_direction, ObjectiveDirection::Maximize);
    ASSERT_EQ(dual_model.number_of_variables(), 2);
    ASSERT_EQ(dual_model.number_of_constraints(), 2);
    EXPECT_EQ(dual_model.objective_coefficients, (std::vector<double>{4, 1}));
    EXPECT_EQ(dual_model.variables_lower_bounds, (std::vector<double>{0, 0}));
    EXPECT_EQ(dual_model.variables_upper_bounds, (std::vector<double>{inf, inf}));
    EXPECT_EQ(dual_model.variable_name(0), "c0_lb");
    EXPECT_EQ(dual_model.variable_name(1), "c1_lb");
    EXPECT_EQ(dual_model.constraint_name(0), "x0");

    // Row x0: y0 + y1 <= 3, row x1: 2 y0 <= 5.
    EXPECT_EQ(dual_model.constraints_starts, (std::vector<int>{0, 2}));
    EXPECT_EQ(dual_model.elements_variables, (std::vector<int>{0, 1, 0}));
    EXPECT_EQ(dual_model.elements_coefficients, (std::vector<double>{1, 1, 2}));
    EXPECT_EQ(dual_model.constraints_lower_bounds, (std::vector<double>{-inf, -inf}));
    EXPECT_EQ(dual_model.constraints_upper_bounds, (std::vector<double>{3, 5}));

    // Optimality certificate: x = (1, 1.5), y = (2.5, 0.5) both have value 10.5.
    EXPECT_TRUE(model.check_solution({1, 1.5}));
    EXPECT_TRUE(dual_model.check_solution({2.5, 0.5}));
    EXPECT_EQ(model.evaluate_objective({1, 1.5}), 10.5);
    EXPECT_EQ(dual_model.evaluate_objective({2.5, 0.5}), 10.5);
}

TEST(MathOptDual, GeneralStructure)
{
    MathOptModel model = general_lp(ObjectiveDirection::Minimize);
    MathOptModel dual_model = dual(model);
    ASSERT_TRUE(dual_model.check(1));

    // Constraints: c0 lb, c1 ub, c2 lb+ub, c3 eq, c4 none.
    // Variables: x0 none (zero bound), x1 lb+ub, x2 ub, x3 none, x4 eq.
    std::vector<std::string> names;
    for (int variable_id = 0; variable_id < dual_model.number_of_variables(); ++variable_id)
        names.push_back(dual_model.variable_name(variable_id));
    EXPECT_EQ(names, (std::vector<std::string>{
                "c0_lb", "c1_ub", "c2_lb", "c2_ub", "c3_eq",
                "x1_lb", "x1_ub", "x2_ub", "x4_eq"}));
    EXPECT_EQ(dual_model.objective_coefficients,
            (std::vector<double>{2, 3, 1, 5, 1.5, -1, 4, 2, 2.5}));
    EXPECT_EQ(dual_model.variables_lower_bounds,
            (std::vector<double>{0, -inf, 0, -inf, -inf, 0, -inf, -inf, -inf}));
    EXPECT_EQ(dual_model.variables_upper_bounds,
            (std::vector<double>{inf, 0, inf, 0, inf, inf, 0, 0, inf}));

    // x0 >= 0 turns its dual row into <= c0, x3 free into = c3.
    EXPECT_EQ(dual_model.constraints_lower_bounds, (std::vector<double>{-inf, 3, -1, 1, 1}));
    EXPECT_EQ(dual_model.constraints_upper_bounds, (std::vector<double>{2, 3, -1, 1, 1}));
}

TEST(MathOptDual, MaximizeFlipsSigns)
{
    MathOptModel dual_model = dual(general_lp(ObjectiveDirection::Maximize));
    ASSERT_TRUE(dual_model.check(1));
    EXPECT_EQ(dual_model.objective_direction, ObjectiveDirection::Minimize);
    EXPECT_EQ(dual_model.variables_lower_bounds,
            (std::vector<double>{-inf, 0, -inf, 0, -inf, -inf, 0, 0, -inf}));
    EXPECT_EQ(dual_model.variables_upper_bounds,
            (std::vector<double>{0, inf, 0, inf, inf, 0, inf, inf, inf}));
    // x0 >= 0 now gives a >= c0 dual row.
    EXPECT_EQ(dual_model.constraints_lower_bounds[0], -2);
    EXPECT_EQ(dual_model.constraints_upper_bounds[0], inf);
}

TEST(MathOptDual, RejectsNonLp)
{
    MathOptModel model = general_lp(ObjectiveDirection::Minimize);
    model.variables_types[0] = VariableType::Integer;
    EXPECT_THROW(dual(model), std::invalid_argument);
}

#if HIGHS_FOUND

namespace
{

double solve_with_highs(MathOptModel model)
{
    model.feasibility_tolerance = 1e-6;
    Highs highs;
    reduce_printout(highs);
    load(highs, model);
    solve(highs);
    EXPECT_EQ(highs.getModelStatus(), HighsModelStatus::kOptimal);
    std::vector<double> solution = get_solution(highs);
    EXPECT_TRUE(model.check_solution(solution, 1));
    return model.evaluate_objective(solution);
}

void expect_strong_duality(const MathOptModel& model)
{
    double primal_value = solve_with_highs(model);
    MathOptModel dual_model = dual(model);
    ASSERT_TRUE(dual_model.check(1));
    EXPECT_NEAR(solve_with_highs(dual_model), primal_value, 1e-6);
    // The dual of the dual has the same optimal value.
    EXPECT_NEAR(solve_with_highs(dual(dual_model)), primal_value, 1e-6);
}

}

TEST(MathOptDual, StrongDualityMinimize)
{
    expect_strong_duality(general_lp(ObjectiveDirection::Minimize));
}

TEST(MathOptDual, StrongDualityMaximize)
{
    expect_strong_duality(general_lp(ObjectiveDirection::Maximize));
}

TEST(MathOptDual, StrongDualityP0033Relaxation)
{
    MathOptModel model = read_mps(std::string(MATHOPTSOLVERSCMAKE_DATA_DIR) + "/miplib3/p0033.mps");
    for (VariableType& type: model.variables_types)
        type = VariableType::Continuous;
    expect_strong_duality(model);
}

#endif

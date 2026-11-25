#include "mathoptsolverscmake/mathopt.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace mathoptsolverscmake;

namespace
{

const double inf = std::numeric_limits<double>::infinity();

/**
 * A model with three variables and one constraint, whose expression tree is
 * given in pre-order (see mathopt.hpp).
 */
MathOptModel tree_model(
        const std::vector<char>& operators,
        const std::vector<double>& values,
        const std::vector<int>& variables,
        const std::vector<int>& number_of_children)
{
    MathOptModel model(3, 1, 0);
    model.constraints_upper_bounds = {4};
    model.constraints_lower_bounds = {-inf};
    model.constraints_starts = {0};
    model.nonlinear_elements_constraints_starts = {0};
    model.nonlinear_elements_operators = operators;
    model.nonlinear_elements_values = values;
    model.nonlinear_elements_variables = variables;
    model.nonlinear_elements_number_of_children = number_of_children;
    return model;
}

std::string format(const MathOptModel& model)
{
    std::ostringstream os;
    model.format_constraint(os, 0);
    return os.str();
}

}

TEST(MathOptExpressions, LinearAndQuadraticTerms)
{
    MathOptModel model(3, 1, 0);
    model.constraints_lower_bounds = {1};
    model.constraints_upper_bounds = {inf};
    model.constraints_starts = {0};
    // 2 x0 - x1 + 0 x2 (the zero term isn't shown) + 3 x0^2 - 0.5 x1 x2.
    model.elements_variables = {0, 1, 2};
    model.elements_coefficients = {2, -1, 0};
    model.quadratic_elements_constraints_starts = {0};
    model.quadratic_elements_variables_1 = {0, 1};
    model.quadratic_elements_variables_2 = {0, 2};
    model.quadratic_elements_coefficients = {3, -0.5};
    EXPECT_EQ(model.format_constraint_expression(0), "2 * x0 - x1 + 3 * x0 ^ 2 - 0.5 * x1 * x2");
    EXPECT_EQ(format(model), "c0: 2 * x0 - x1 + 3 * x0 ^ 2 - 0.5 * x1 * x2 >= 1\n");
    EXPECT_EQ(model.constraint_degree(0), 2);

    model.objective_coefficients = {0, 0.1, -1};
    EXPECT_EQ(model.format_objective_expression(), "0.1 * x1 - x2");
    EXPECT_EQ(model.objective_degree(), 1);
}

TEST(MathOptExpressions, Bounds)
{
    MathOptModel model(1, 1, 0);
    model.constraints_starts = {0};
    model.elements_variables = {0};
    model.elements_coefficients = {1};
    std::vector<std::pair<std::pair<double, double>, std::string>> cases = {
        {{-inf, 4}, "c0: x0 <= 4\n"},
        {{1, inf}, "c0: x0 >= 1\n"},
        {{2, 2}, "c0: x0 = 2\n"},
        {{-1, 3}, "c0: -1 <= x0 <= 3\n"},
        {{-inf, inf}, "c0: x0\n"},
    };
    for (const auto& c: cases) {
        model.constraints_lower_bounds = {c.first.first};
        model.constraints_upper_bounds = {c.first.second};
        EXPECT_EQ(format(model), c.second);
    }
    model.elements_coefficients = {0};
    model.constraints_lower_bounds = {-inf};
    model.constraints_upper_bounds = {0};
    EXPECT_EQ(format(model), "c0: 0 <= 0\n");
}

TEST(MathOptExpressions, Precedence)
{
    // x0 * (x1 + x2)
    MathOptModel model = tree_model(
            {'*', 'v', '+', 'v', 'v'}, {0, 0, 0, 0, 0}, {-1, 0, -1, 1, 2}, {2, 0, 2, 0, 0});
    EXPECT_EQ(model.format_constraint_expression(0), "x0 * (x1 + x2)");
    EXPECT_EQ(model.constraint_degree(0), 2);

    // x0 - (x1 - x2)
    model = tree_model({'-', 'v', '-', 'v', 'v'}, {0, 0, 0, 0, 0}, {-1, 0, -1, 1, 2}, {2, 0, 2, 0, 0});
    EXPECT_EQ(model.format_constraint_expression(0), "x0 - (x1 - x2)");

    // x0 / (x1 * x2)
    model = tree_model({'/', 'v', '*', 'v', 'v'}, {0, 0, 0, 0, 0}, {-1, 0, -1, 1, 2}, {2, 0, 2, 0, 0});
    EXPECT_EQ(model.format_constraint_expression(0), "x0 / (x1 * x2)");
    EXPECT_EQ(model.constraint_degree(0), inf);

    // (-x0) ^ 2 + x1 ^ -1.5: the base and the exponent are parenthesized.
    model = tree_model(
            {'+', 'p', 'n', 'v', 'k', 'p', 'v', 'k'},
            {0, 0, 0, 0, 2, 0, 0, -1.5},
            {-1, -1, -1, 0, -1, -1, 1, -1},
            {2, 2, 1, 0, 0, 2, 0, 0});
    EXPECT_EQ(model.format_constraint_expression(0), "(-x0) ^ 2 + x1 ^ (-1.5)");
    EXPECT_EQ(model.constraint_degree(0), inf);

    // -(x0 + x1), and a sum with negative terms: x2 + -3 * x0 + -x1.
    model = tree_model({'n', '+', 'v', 'v'}, {0, 0, 0, 0}, {-1, -1, 0, 1}, {1, 2, 0, 0});
    EXPECT_EQ(model.format_constraint_expression(0), "-(x0 + x1)");
    model = tree_model(
            {'+', 'v', '*', 'k', 'v', 'n', 'v'},
            {0, 0, 0, -3, 0, 0, 0},
            {-1, 2, -1, -1, 0, -1, 1},
            {3, 0, 2, 0, 0, 1, 0});
    EXPECT_EQ(model.format_constraint_expression(0), "x2 - 3 * x0 - x1");
}

TEST(MathOptExpressions, Functions)
{
    // exp(x0) + log(x1) + max(x0, x1, 2)
    MathOptModel model = tree_model(
            {'+', 'e', 'v', 'l', 'v', 'M', 'v', 'v', 'k'},
            {0, 0, 0, 0, 0, 0, 0, 0, 2},
            {-1, -1, 0, -1, 1, -1, 0, 1, -1},
            {3, 1, 0, 1, 0, 3, 0, 0, 0});
    EXPECT_EQ(model.format_constraint_expression(0), "exp(x0) + log(x1) + max(x0, x1, 2)");
    EXPECT_EQ(model.constraint_degree(0), inf);

    // sqrt(4): a constant.
    model = tree_model({'q', 'k'}, {0, 4}, {-1, -1}, {1, 0});
    EXPECT_EQ(model.constraint_degree(0), 0);
}

TEST(MathOptExpressions, PolynomialDegrees)
{
    // (x0 * x1 + x2) ^ 3: degree 6.
    MathOptModel model = tree_model(
            {'p', '+', '*', 'v', 'v', 'v', 'k'},
            {0, 0, 0, 0, 0, 0, 3},
            {-1, -1, -1, 0, 1, 2, -1},
            {2, 2, 2, 0, 0, 0, 0});
    EXPECT_EQ(model.constraint_degree(0), 6);
    // x0 / 2: degree 1.
    model = tree_model({'/', 'v', 'k'}, {0, 0, 2}, {-1, 0, -1}, {2, 0, 0});
    EXPECT_EQ(model.constraint_degree(0), 1);
}

TEST(MathOptExpressions, DefinedVariables)
{
    // d0 = x0 * x1; d1 = d0 + x2; constraint: d1 ^ 2.
    MathOptModel model = tree_model({'p', 'd', 'k'}, {0, 0, 2}, {-1, 1, -1}, {2, 0, 0});
    model.defined_variables_starts = {0, 3};
    model.defined_variables_elements_operators = {'*', 'v', 'v', '+', 'd', 'v'};
    model.defined_variables_elements_values = {0, 0, 0, 0, 0, 0};
    model.defined_variables_elements_variables = {-1, 0, 1, -1, 0, 2};
    model.defined_variables_elements_number_of_children = {2, 0, 0, 2, 0, 0};
    EXPECT_EQ(model.format_constraint_expression(0), "d1 ^ 2");
    EXPECT_EQ(model.format_constraint_expression(0, true), "(x0 * x1 + x2) ^ 2");
    EXPECT_EQ(model.format_defined_variable_expression(1), "d0 + x2");
    EXPECT_EQ(model.format_defined_variable_expression(1, true), "x0 * x1 + x2");
    EXPECT_EQ(model.constraint_degree(0), 4);
}

TEST(MathOptExpressions, InvalidIndices)
{
    MathOptModel model(1, 1, 0);
    model.constraints_starts = {0};
    EXPECT_THROW(model.format_constraint_expression(1), std::invalid_argument);
    EXPECT_THROW(model.constraint_degree(-1), std::invalid_argument);
    EXPECT_THROW(model.format_defined_variable_expression(0), std::invalid_argument);
}

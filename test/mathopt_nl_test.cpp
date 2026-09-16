#include "mathoptsolverscmake/mathopt_nl.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <sstream>

using namespace mathoptsolverscmake;

namespace
{

const double inf = std::numeric_limits<double>::infinity();

MathOptModel read(const std::string& nl)
{
    std::istringstream in(nl);
    return read_nl(in);
}

std::string write(const MathOptModel& model)
{
    std::ostringstream out;
    write_nl(model, out);
    return out.str();
}

}

TEST(MathOptNl, LinearBoundTypesAndSenses)
{
    // minimize x0 + 2 x1
    // s.t. c0: x0 + x1 <= 5       (bound type 1: upper only)
    //      c1: x0 - x1 = 0        (bound type 4: equality)
    // x0 in [0, 10] (type 0: range), x1 free (type 3).
    MathOptModel model = read(
            "g0\n"
            "2 2 1 0 1\n"
            "0 0\n"
            "0 0\n"
            "0 0 0\n"
            "0 0\n"
            "0 0 0 0 0\n"
            "4 2\n"
            "0 0\n"
            "0 0 0 0 0\n"
            "b\n"
            "0 0 10\n"
            "3\n"
            "r\n"
            "1 5\n"
            "4 0\n"
            "C0\n"
            "n0\n"
            "C1\n"
            "n0\n"
            "O0 0\n"
            "n0\n"
            "J0 2\n"
            "0 1\n"
            "1 1\n"
            "J1 2\n"
            "0 1\n"
            "1 -1\n"
            "G0 2\n"
            "0 1\n"
            "1 2\n");

    EXPECT_EQ(model.objective_direction, ObjectiveDirection::Minimize);
    ASSERT_EQ(model.number_of_variables(), 2);
    ASSERT_EQ(model.number_of_constraints(), 2);

    EXPECT_EQ(model.variables_lower_bounds[0], 0.0);
    EXPECT_EQ(model.variables_upper_bounds[0], 10.0);
    EXPECT_EQ(model.variables_lower_bounds[1], -inf);
    EXPECT_EQ(model.variables_upper_bounds[1], inf);

    EXPECT_EQ(model.objective_coefficients[0], 1.0);
    EXPECT_EQ(model.objective_coefficients[1], 2.0);

    EXPECT_EQ(model.constraints_lower_bounds[0], -inf);
    EXPECT_EQ(model.constraints_upper_bounds[0], 5.0);
    EXPECT_EQ(model.constraints_lower_bounds[1], 0.0);
    EXPECT_EQ(model.constraints_upper_bounds[1], 0.0);

    std::vector<double> x = {1.0, 2.0};
    EXPECT_EQ(model.evaluate_objective(x), 1.0 + 2.0 * 2.0);
    EXPECT_EQ(model.evaluate_constraint(x, 0), 3.0);
    EXPECT_EQ(model.evaluate_constraint(x, 1), -1.0);
}

TEST(MathOptNl, RangeConstraintAndFixedVariable)
{
    // c0: 2 <= x0 + x1 <= 8 (bound type 0), x0 fixed at 3 (bound type 4).
    MathOptModel model = read(
            "g0\n"
            "2 1 1 0 0\n"
            "0 0\n"
            "0 0\n"
            "0 0 0\n"
            "0 0\n"
            "0 0 0 0 0\n"
            "2 0\n"
            "0 0\n"
            "0 0 0 0 0\n"
            "b\n"
            "4 3\n"
            "2 0\n"
            "r\n"
            "0 2 8\n"
            "C0\n"
            "n0\n"
            "O0 0\n"
            "n0\n"
            "J0 2\n"
            "0 1\n"
            "1 1\n");

    ASSERT_EQ(model.number_of_variables(), 2);
    EXPECT_EQ(model.variables_lower_bounds[0], 3.0);
    EXPECT_EQ(model.variables_upper_bounds[0], 3.0);
    EXPECT_EQ(model.variables_lower_bounds[1], 0.0);
    EXPECT_EQ(model.variables_upper_bounds[1], inf);
    EXPECT_EQ(model.constraints_lower_bounds[0], 2.0);
    EXPECT_EQ(model.constraints_upper_bounds[0], 8.0);
}

TEST(MathOptNl, NonlinearObjective)
{
    // minimize exp(x0) + x1
    // s.t. c0: x0 + x1 <= 5
    // x0 in [0, 10], x1 in [0, inf).
    MathOptModel model = read(
            "g0\n"
            "2 1 1 0 0\n"
            "0 1\n"
            "0 0\n"
            "0 1 0\n"
            "0 0\n"
            "0 0 0 0 0\n"
            "2 1\n"
            "0 0\n"
            "0 0 0 0 0\n"
            "b\n"
            "0 0 10\n"
            "2 0\n"
            "r\n"
            "1 5\n"
            "C0\n"
            "n0\n"
            "O0 0\n"
            "o44\n"
            "v0\n"
            "J0 2\n"
            "0 1\n"
            "1 1\n"
            "G0 1\n"
            "1 1\n");

    ASSERT_EQ(model.number_of_variables(), 2);
    EXPECT_TRUE(model.has_nonlinear());
    EXPECT_FALSE(model.has_quadratic());

    std::vector<double> x = {0.0, 3.0};
    EXPECT_NEAR(model.evaluate_objective(x), std::exp(0.0) + 3.0, 1e-12);
    x = {2.0, 3.0};
    EXPECT_NEAR(model.evaluate_objective(x), std::exp(2.0) + 3.0, 1e-12);
}

TEST(MathOptNl, LinearBinaryFromHeaderNbv)
{
    // A variable declared "binary" via the header's nbv count, but whose
    // own 'b' segment bound is [0, +infinity): read_nl() must trust the
    // header's explicit nbv/niv split for the linear portion, not infer
    // Binary-vs-Integer from bounds (that heuristic only applies to the
    // nonlinear-variable positions, which have no such explicit signal).
    MathOptModel model = read(
            "g0\n"
            "1 0 1 0 0\n"
            "0 0\n"
            "0 0\n"
            "0 0 0\n"
            "0 0\n"
            "1 0 0 0 0\n"
            "0 1\n"
            "0 0\n"
            "0 0 0 0 0\n"
            "b\n"
            "2 0\n"
            "r\n"
            "O0 0\n"
            "n0\n"
            "G0 1\n"
            "0 1\n");

    ASSERT_EQ(model.number_of_variables(), 1);
    EXPECT_EQ(model.variables_types[0], VariableType::Binary);
    EXPECT_EQ(model.variables_lower_bounds[0], 0.0);
    EXPECT_EQ(model.variables_upper_bounds[0], inf);
}

TEST(MathOptNl, BinaryFormatThrows)
{
    EXPECT_THROW(read("b1\n"), std::invalid_argument);
}

TEST(MathOptNl, CommonSubexpressionsThrow)
{
    EXPECT_THROW(
            read(
                    "g0\n"
                    "1 0 1 0 0\n"
                    "0 0\n"
                    "0 0\n"
                    "0 0 0\n"
                    "0 0\n"
                    "0 0 0 0 0\n"
                    "0 1\n"
                    "0 0\n"
                    "1 0 0 0 0\n"),
            std::invalid_argument);
}

TEST(MathOptNl, UnsupportedOpcodeThrows)
{
    // o12 is MAXLIST, not one of MathOptModel's supported operators.
    EXPECT_THROW(
            read(
                    "g0\n"
                    "1 0 1 0 0\n"
                    "0 1\n"
                    "0 0\n"
                    "0 1 0\n"
                    "0 0\n"
                    "0 0 0 0 0\n"
                    "0 0\n"
                    "0 0\n"
                    "0 0 0 0 0\n"
                    "b\n"
                    "3\n"
                    "r\n"
                    "O0 0\n"
                    "o12\n"
                    "2\n"
                    "v0\n"
                    "n1\n"),
            std::invalid_argument);
}

TEST(MathOptNl, MissingFileThrows)
{
    EXPECT_THROW(read_nl("/no/such/file.nl"), std::invalid_argument);
}

TEST(MathOptNl, WriteThenReadLinearRoundTrip)
{
    // All-continuous linear model: write_nl()'s variable partition is the
    // identity in this case (no nonlinear content to sort to the front),
    // so an exact field-by-field comparison is meaningful, not just a
    // semantic (evaluate_*) one.
    MathOptModel original(2, 2, 4);
    original.objective_direction = ObjectiveDirection::Maximize;
    original.variables_lower_bounds = {0.0, -5.0};
    original.variables_upper_bounds = {10.0, 5.0};
    original.objective_coefficients = {3.0, 2.0};
    original.constraints_starts = {0, 2};
    original.elements_variables = {0, 1, 0, 1};
    original.elements_coefficients = {1.0, 1.0, 1.0, -1.0};
    original.constraints_lower_bounds = {-inf, 0.0};
    original.constraints_upper_bounds = {4.0, 0.0};

    MathOptModel round_tripped = read(write(original));

    ASSERT_EQ(round_tripped.number_of_variables(), 2);
    ASSERT_EQ(round_tripped.number_of_constraints(), 2);
    EXPECT_EQ(round_tripped.objective_direction, original.objective_direction);
    for (int v = 0; v < 2; ++v) {
        EXPECT_EQ(round_tripped.variables_lower_bounds[v], original.variables_lower_bounds[v]);
        EXPECT_EQ(round_tripped.variables_upper_bounds[v], original.variables_upper_bounds[v]);
        EXPECT_EQ(round_tripped.objective_coefficients[v], original.objective_coefficients[v]);
    }
    for (int c = 0; c < 2; ++c) {
        EXPECT_EQ(round_tripped.constraints_lower_bounds[c], original.constraints_lower_bounds[c]);
        EXPECT_EQ(round_tripped.constraints_upper_bounds[c], original.constraints_upper_bounds[c]);
    }
}

TEST(MathOptNl, WriteThenReadNonlinearAndQuadraticRoundTrip)
{
    // minimize exp(x0) + x1 * x2
    // s.t. c0: x0 + x1 + x2 <= 10
    // .nl reorders variables (nonlinear ones first), so compare
    // semantically via evaluate_objective/evaluate_constraint instead of
    // by index.
    MathOptModel original(3, 1, 3);
    original.objective_direction = ObjectiveDirection::Minimize;
    original.variables_lower_bounds = {0.0, 0.0, 0.0};
    original.variables_upper_bounds = {5.0, 5.0, 5.0};
    original.objective_quadratic_elements_variables_1 = {1};
    original.objective_quadratic_elements_variables_2 = {2};
    original.objective_quadratic_elements_coefficients = {1.0};
    original.objective_nonlinear_elements_operators = {'e', 'v'};
    original.objective_nonlinear_elements_values = {0.0, 0.0};
    original.objective_nonlinear_elements_variables = {-1, 0};
    original.objective_nonlinear_elements_number_of_children = {1, 0};
    original.constraints_starts = {0};
    original.elements_variables = {0, 1, 2};
    original.elements_coefficients = {1.0, 1.0, 1.0};
    original.constraints_lower_bounds = {-inf};
    original.constraints_upper_bounds = {10.0};

    MathOptModel round_tripped = read(write(original));

    ASSERT_EQ(round_tripped.number_of_variables(), 3);
    ASSERT_EQ(round_tripped.number_of_constraints(), 1);
    EXPECT_TRUE(round_tripped.check_solution({1.0, 2.0, 3.0}));

    std::vector<double> values = {0.5, 1.5, 2.5};
    EXPECT_NEAR(
            round_tripped.evaluate_objective(values),
            original.evaluate_objective(values),
            1e-9);
    EXPECT_NEAR(
            round_tripped.evaluate_constraint(values, 0),
            original.evaluate_constraint(values, 0),
            1e-9);
}

TEST(MathOptNl, WriteBlackBoxThrows)
{
    MathOptModel model(1, 0, 0);
    model.objective_direction = ObjectiveDirection::Minimize;
    model.objective_function = [](const std::vector<double>& x) -> BlackBoxFunctionOutput
    {
        BlackBoxFunctionOutput output;
        output.objective_value = x[0];
        return output;
    };

    std::ostringstream out;
    EXPECT_THROW(write_nl(model, out), std::invalid_argument);
}

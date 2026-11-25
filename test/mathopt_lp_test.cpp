#include "mathoptsolverscmake/mathopt_lp.hpp"

#include <gtest/gtest.h>

#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace mathoptsolverscmake;

namespace
{

const double inf = std::numeric_limits<double>::infinity();

MathOptModel read(const std::string& lp)
{
    std::istringstream in(lp);
    return read_lp(in);
}

std::string write(const MathOptModel& model)
{
    std::ostringstream out;
    write_lp(model, out);
    return out.str();
}

}

TEST(MathOptLp, LinearRowTypesAndBounds)
{
    // maximize 3 x0 + 2 x1
    // s.t. x0 + x1 <= 4   (L)
    //      x0 + 2 x1 >= 1  (G)
    //      x0 - x1 = 0     (E)
    //      0 <= x0 <= 10
    //      x1 free
    MathOptModel model = read(
            "\\ a comment\n"
            "Maximize\n"
            " obj: 3 x0 + 2 x1\n"
            "Subject To\n"
            " c1: x0 + x1 <= 4\n"
            " c2: x0 + 2 x1 >= 1\n"
            " c3: x0 - x1 = 0\n"
            "Bounds\n"
            " x0 <= 10\n"
            " x1 free\n"
            "End\n");

    EXPECT_EQ(model.objective_direction, ObjectiveDirection::Maximize);
    ASSERT_EQ(model.number_of_variables(), 2);
    ASSERT_EQ(model.number_of_constraints(), 3);

    EXPECT_EQ(model.variable_name(0), "x0");
    EXPECT_EQ(model.variable_name(1), "x1");
    EXPECT_EQ(model.variables_lower_bounds[0], 0.0);
    EXPECT_EQ(model.variables_upper_bounds[0], 10.0);
    EXPECT_EQ(model.variables_lower_bounds[1], -inf);
    EXPECT_EQ(model.variables_upper_bounds[1], inf);

    EXPECT_EQ(model.objective_coefficients[0], 3.0);
    EXPECT_EQ(model.objective_coefficients[1], 2.0);

    EXPECT_EQ(model.constraints_lower_bounds[0], -inf);
    EXPECT_EQ(model.constraints_upper_bounds[0], 4.0);
    EXPECT_EQ(model.constraints_lower_bounds[1], 1.0);
    EXPECT_EQ(model.constraints_upper_bounds[1], inf);
    EXPECT_EQ(model.constraints_lower_bounds[2], 0.0);
    EXPECT_EQ(model.constraints_upper_bounds[2], 0.0);
    EXPECT_EQ(model.constraint_name(0), "c1");
}

TEST(MathOptLp, RangeConstraintBothDirections)
{
    MathOptModel model = read(
            "Minimize\n"
            " obj: x0\n"
            "Subject To\n"
            " c1: 2 <= x0 + x1 <= 8\n"
            " c2: 8 >= x0 - x1 >= 2\n"
            "End\n");

    ASSERT_EQ(model.number_of_constraints(), 2);
    EXPECT_EQ(model.constraints_lower_bounds[0], 2.0);
    EXPECT_EQ(model.constraints_upper_bounds[0], 8.0);
    EXPECT_EQ(model.constraints_lower_bounds[1], 2.0);
    EXPECT_EQ(model.constraints_upper_bounds[1], 8.0);
}

TEST(MathOptLp, GeneralAndBinarySections)
{
    MathOptModel model = read(
            "Minimize\n"
            " obj: x0 + x1 + x2\n"
            "Subject To\n"
            " c1: x0 + x1 + x2 <= 10\n"
            "General\n"
            " x0\n"
            "Binary\n"
            " x1\n"
            "Bounds\n"
            " x2 <= 20\n"
            "General\n"
            " x2\n"
            "End\n");

    ASSERT_EQ(model.number_of_variables(), 3);
    // General with no bounds entry: default stays [0, +infinity), unlike MPS.
    EXPECT_EQ(model.variables_types[0], VariableType::Integer);
    EXPECT_EQ(model.variables_lower_bounds[0], 0.0);
    EXPECT_EQ(model.variables_upper_bounds[0], inf);
    // Binary with no bounds entry: upper defaults to 1.
    EXPECT_EQ(model.variables_types[1], VariableType::Binary);
    EXPECT_EQ(model.variables_upper_bounds[1], 1.0);
    // General with an explicit bound: type is Integer regardless of bounds.
    EXPECT_EQ(model.variables_types[2], VariableType::Integer);
    EXPECT_EQ(model.variables_upper_bounds[2], 20.0);
}

TEST(MathOptLp, BoundForms)
{
    MathOptModel model = read(
            "Minimize\n"
            " obj: x0 + x1 + x2 + x3\n"
            "Subject To\n"
            " c1: x0 + x1 + x2 + x3 <= 100\n"
            "Bounds\n"
            " x0 = 5\n"
            " -10 <= x1 <= 10\n"
            " -infinity <= x2 <= 5\n"
            " x3 >= -3\n"
            "End\n");

    ASSERT_EQ(model.number_of_variables(), 4);
    EXPECT_EQ(model.variables_lower_bounds[0], 5.0);
    EXPECT_EQ(model.variables_upper_bounds[0], 5.0);
    EXPECT_EQ(model.variables_lower_bounds[1], -10.0);
    EXPECT_EQ(model.variables_upper_bounds[1], 10.0);
    EXPECT_EQ(model.variables_lower_bounds[2], -inf);
    EXPECT_EQ(model.variables_upper_bounds[2], 5.0);
    EXPECT_EQ(model.variables_lower_bounds[3], -3.0);
    EXPECT_EQ(model.variables_upper_bounds[3], inf);
}

TEST(MathOptLp, DuplicateVariableInExpressionIsSummed)
{
    MathOptModel model = read(
            "Minimize\n"
            " obj: x0 + 2 x0\n"
            "Subject To\n"
            " c1: x0 + x0 <= 10\n"
            "End\n");

    ASSERT_EQ(model.number_of_variables(), 1);
    EXPECT_EQ(model.objective_coefficients[0], 3.0);
    ASSERT_EQ(model.number_of_variables(0), 1);
    EXPECT_EQ(model.elements_coefficients[0], 2.0);
}

TEST(MathOptLp, UnsupportedSectionThrows)
{
    EXPECT_THROW(
            read(
                    "Minimize\n"
                    " obj: x0\n"
                    "Subject To\n"
                    " c1: x0 <= 10\n"
                    "SOS\n"
                    " s1: S1:: x0:1\n"
                    "End\n"),
            std::invalid_argument);
}

TEST(MathOptLp, QuadraticObjectiveAndConstraints)
{
    // min x + 1/2 [ 2 x^2 + 2 x y ] = x + x^2 + x y
    // s.t. c1: x + y >= 1
    //      q1: y + [ x^2 + 3 x y - y^2 ] <= 4
    MathOptModel model = read(
            "Minimize\n"
            " obj: x + [ 2 x ^ 2 + 2 x * y ] / 2\n"
            "Subject To\n"
            " c1: x + y >= 1\n"
            " q1: y + [ x^2 + 3 x * y - y ^ 2 ] <= 4\n"
            "End\n");

    ASSERT_EQ(model.number_of_variables(), 2);
    ASSERT_EQ(model.number_of_constraints(), 2);
    EXPECT_EQ(model.objective_coefficients, (std::vector<double>{1.0, 0.0}));
    EXPECT_EQ(model.objective_quadratic_elements_variables_1, (std::vector<int>{0, 0}));
    EXPECT_EQ(model.objective_quadratic_elements_variables_2, (std::vector<int>{0, 1}));
    EXPECT_EQ(model.objective_quadratic_elements_coefficients, (std::vector<double>{1.0, 1.0}));
    EXPECT_EQ(model.quadratic_elements_constraints_starts, (std::vector<int>{0, 0}));
    EXPECT_EQ(model.quadratic_elements_variables_1, (std::vector<int>{0, 0, 1}));
    EXPECT_EQ(model.quadratic_elements_variables_2, (std::vector<int>{0, 1, 1}));
    EXPECT_EQ(model.quadratic_elements_coefficients, (std::vector<double>{1.0, 3.0, -1.0}));
    EXPECT_EQ(model.number_of_variables(1), 1);
    EXPECT_EQ(model.constraints_upper_bounds[1], 4.0);
    EXPECT_TRUE(model.check());

    // At x = 1, y = 2: objective 1 + 1 + 2 = 4; q1's body 2 + 1 + 6 - 4 = 5.
    EXPECT_EQ(model.evaluate_objective({1.0, 2.0}), 4.0);
    EXPECT_EQ(model.evaluate_constraint({1.0, 2.0}, 1), 5.0);
}

TEST(MathOptLp, LinearModelHasNoQuadraticStarts)
{
    MathOptModel model = read(
            "Minimize\n"
            " obj: x\n"
            "Subject To\n"
            " c1: x >= 1\n"
            "End\n");
    EXPECT_FALSE(model.has_quadratic());
    EXPECT_TRUE(model.quadratic_elements_constraints_starts.empty());
}

TEST(MathOptLp, MalformedQuadraticTermsThrow)
{
    for (const char* objective: {
            "[ x ^ 3 ] / 2", "[ x y ] / 2", "[ x * 2 ] / 2", "[ x ^ 2 / 2", "[ x ^ 2 ] / 0", "[ x ^ 2 x ^ 2 ]"}) {
        EXPECT_THROW(read(
                "Minimize\n"
                " obj: " + std::string(objective) + "\n"
                "Subject To\n"
                " c1: x >= 1\n"
                "End\n"), std::invalid_argument) << objective;
    }
}

TEST(MathOptLp, QuadraticRoundTrip)
{
    MathOptModel original = read(
            "Maximize\n"
            " obj: 2 x - [ 3 x ^ 2 - 2 x * y + 0.5 y ^ 2 ] / 2\n"
            "Subject To\n"
            " c1: x + y <= 4\n"
            " q1: [ x ^ 2 + y ^ 2 ] <= 9\n"
            " c2: x - y >= -1\n"
            " q2: -2 <= x + [ - x * y ] <= 3\n"
            "End\n");

    std::string text = write(original);
    MathOptModel round_tripped = read(text);

    EXPECT_EQ(round_tripped.objective_coefficients, original.objective_coefficients);
    EXPECT_EQ(round_tripped.objective_quadratic_elements_variables_1, original.objective_quadratic_elements_variables_1);
    EXPECT_EQ(round_tripped.objective_quadratic_elements_variables_2, original.objective_quadratic_elements_variables_2);
    EXPECT_EQ(round_tripped.objective_quadratic_elements_coefficients, original.objective_quadratic_elements_coefficients);
    EXPECT_EQ(round_tripped.quadratic_elements_constraints_starts, original.quadratic_elements_constraints_starts);
    EXPECT_EQ(round_tripped.quadratic_elements_variables_1, original.quadratic_elements_variables_1);
    EXPECT_EQ(round_tripped.quadratic_elements_variables_2, original.quadratic_elements_variables_2);
    EXPECT_EQ(round_tripped.quadratic_elements_coefficients, original.quadratic_elements_coefficients);
    EXPECT_EQ(round_tripped.constraints_lower_bounds, original.constraints_lower_bounds);
    EXPECT_EQ(round_tripped.constraints_upper_bounds, original.constraints_upper_bounds);
    EXPECT_NE(text.find("obj: +2 x + [ -3 x ^ 2 +2 x * y -0.5 y ^ 2 ] / 2"), std::string::npos) << text;
    EXPECT_NE(text.find("q1: + [ + x ^ 2 + y ^ 2 ] <=9"), std::string::npos) << text;
}

TEST(MathOptLp, MissingFileThrows)
{
    EXPECT_THROW(read_lp("/no/such/file.lp"), std::invalid_argument);
}

TEST(MathOptLp, KeywordAsConstraintName)
{
    // A section keyword followed by ':' is a constraint name (as in HiGHS's
    // check/instances/1451.lp).
    MathOptModel model = read(
            "Minimize\n"
            "  x\n"
            "Subject to\n"
            " end: x >= 3\n"
            "Integer\n"
            "  x\n"
            "End\n");

    ASSERT_EQ(model.number_of_variables(), 1);
    EXPECT_EQ(model.variables_types[0], VariableType::Integer);
    ASSERT_EQ(model.number_of_constraints(), 1);
    EXPECT_EQ(model.constraint_name(0), "end");
    EXPECT_EQ(model.constraints_lower_bounds[0], 3.0);
    EXPECT_EQ(model.constraints_upper_bounds[0], inf);
}

TEST(MathOptLp, KeywordConstraintNamesRoundTrip)
{
    std::vector<std::string> names = {
        "end", "bounds", "st", "subject", "general", "binary", "free", "minimize", "semi", "sos"};
    MathOptModel original;
    original.variables_names = {"x"};
    original.variables_lower_bounds = {0.0};
    original.variables_upper_bounds = {inf};
    original.variables_types = {VariableType::Integer};
    original.objective_coefficients = {1.0};
    for (const std::string& name: names) {
        original.constraints_names.push_back(name);
        original.constraints_lower_bounds.push_back(1.0);
        original.constraints_upper_bounds.push_back(inf);
        original.constraints_starts.push_back((int)original.elements_variables.size());
        original.elements_variables.push_back(0);
        original.elements_coefficients.push_back(1.0);
    }

    MathOptModel round_tripped = read(write(original));

    ASSERT_EQ(round_tripped.number_of_constraints(), (int)names.size());
    for (int constraint_id = 0; constraint_id < (int)names.size(); ++constraint_id)
        EXPECT_EQ(round_tripped.constraint_name(constraint_id), names[constraint_id]);
    EXPECT_EQ(round_tripped.variables_types[0], VariableType::Integer);
}

TEST(MathOptLp, WriteThenReadRoundTrip)
{
    MathOptModel original = read(
            "Maximize\n"
            " obj: 3 x0 + 2 x1\n"
            "Subject To\n"
            " c1: x0 + x1 <= 4\n"
            " c2: 2 <= x0 - x1 <= 6\n"
            " c3: x0 - x1 = 0\n"
            "Bounds\n"
            " x0 <= 10\n"
            " -5 <= x1 <= 5\n"
            "End\n");

    MathOptModel round_tripped = read(write(original));

    ASSERT_EQ(round_tripped.number_of_variables(), original.number_of_variables());
    ASSERT_EQ(round_tripped.number_of_constraints(), original.number_of_constraints());
    EXPECT_EQ(round_tripped.objective_direction, original.objective_direction);
    for (int variable_id = 0; variable_id < original.number_of_variables(); ++variable_id) {
        EXPECT_EQ(round_tripped.variable_name(variable_id), original.variable_name(variable_id));
        EXPECT_EQ(round_tripped.variables_lower_bounds[variable_id], original.variables_lower_bounds[variable_id]);
        EXPECT_EQ(round_tripped.variables_upper_bounds[variable_id], original.variables_upper_bounds[variable_id]);
        EXPECT_EQ(round_tripped.variables_types[variable_id], original.variables_types[variable_id]);
        EXPECT_EQ(round_tripped.objective_coefficients[variable_id], original.objective_coefficients[variable_id]);
    }
    for (int constraint_id = 0; constraint_id < original.number_of_constraints(); ++constraint_id) {
        EXPECT_EQ(round_tripped.constraint_name(constraint_id), original.constraint_name(constraint_id));
        EXPECT_EQ(round_tripped.constraints_lower_bounds[constraint_id], original.constraints_lower_bounds[constraint_id]);
        EXPECT_EQ(round_tripped.constraints_upper_bounds[constraint_id], original.constraints_upper_bounds[constraint_id]);
        EXPECT_EQ(round_tripped.number_of_variables(constraint_id), original.number_of_variables(constraint_id));
    }
}

TEST(MathOptLp, WriteFreeConstraintRoundTrip)
{
    MathOptModel model(1, 1, 1);
    model.objective_direction = ObjectiveDirection::Minimize;
    model.objective_coefficients[0] = 1.0;
    model.constraints_starts[0] = 0;
    model.elements_variables[0] = 0;
    model.elements_coefficients[0] = 1.0;
    // Both bounds infinite: a Free constraint.
    ASSERT_EQ(model.constraint_sense(0), ConstraintSense::Free);

    MathOptModel round_tripped = read(write(model));

    ASSERT_EQ(round_tripped.number_of_constraints(), 1);
    EXPECT_EQ(round_tripped.constraint_sense(0), ConstraintSense::Free);
}

TEST(MathOptLp, WriteIntegerAndBinaryRoundTrip)
{
    MathOptModel model(3, 0, 0);
    model.objective_direction = ObjectiveDirection::Minimize;
    model.objective_coefficients = {1.0, 1.0, 1.0};
    model.variables_lower_bounds[0] = 0.0;
    model.variables_upper_bounds[0] = 1.0;
    model.variables_types[0] = VariableType::Binary;
    model.variables_lower_bounds[1] = 0.0;
    model.variables_upper_bounds[1] = 20.0;
    model.variables_types[1] = VariableType::Integer;
    model.variables_lower_bounds[2] = 0.0;
    model.variables_upper_bounds[2] = 1.0;
    model.variables_types[2] = VariableType::Integer;

    MathOptModel round_tripped = read(write(model));

    ASSERT_EQ(round_tripped.number_of_variables(), 3);
    EXPECT_EQ(round_tripped.variables_types[0], VariableType::Binary);
    EXPECT_EQ(round_tripped.variables_types[1], VariableType::Integer);
    EXPECT_EQ(round_tripped.variables_upper_bounds[1], 20.0);
    // Integer (not Binary) with bounds [0, 1]: section placement, not bounds,
    // determines the type, so it stays Integer (unlike the MPS reader).
    EXPECT_EQ(round_tripped.variables_types[2], VariableType::Integer);
    EXPECT_EQ(round_tripped.variables_upper_bounds[2], 1.0);
}

TEST(MathOptLp, WriteNonlinearThrows)
{
    // Objective: exp(x0).
    MathOptModel model(1, 0, 0);
    model.objective_direction = ObjectiveDirection::Minimize;
    model.objective_nonlinear_elements_operators = {'e', 'v'};
    model.objective_nonlinear_elements_values = {0.0, 0.0};
    model.objective_nonlinear_elements_variables = {-1, 0};
    model.objective_nonlinear_elements_number_of_children = {1, 0};

    std::ostringstream out;
    EXPECT_THROW(write_lp(model, out), std::invalid_argument);
}

TEST(MathOptLp, WriteInvalidNameThrows)
{
    MathOptModel model(1, 0, 0);
    model.objective_direction = ObjectiveDirection::Minimize;
    model.objective_coefficients[0] = 1.0;
    model.variables_names = {"bad-name"};

    std::ostringstream out;
    EXPECT_THROW(write_lp(model, out), std::invalid_argument);
}

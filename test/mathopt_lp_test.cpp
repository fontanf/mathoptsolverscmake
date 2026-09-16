#include "mathoptsolverscmake/mathopt_lp.hpp"

#include <gtest/gtest.h>

#include <limits>
#include <sstream>

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

TEST(MathOptLp, QuadraticBracketThrows)
{
    EXPECT_THROW(
            read(
                    "Minimize\n"
                    " obj: [ x0 * x0 ]\n"
                    "Subject To\n"
                    " c1: x0 <= 10\n"
                    "End\n"),
            std::invalid_argument);
}

TEST(MathOptLp, MissingFileThrows)
{
    EXPECT_THROW(read_lp("/no/such/file.lp"), std::invalid_argument);
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

TEST(MathOptLp, WriteNonMilpThrows)
{
    MathOptModel model(1, 0, 0);
    model.objective_direction = ObjectiveDirection::Minimize;
    model.objective_quadratic_elements_variables_1.push_back(0);
    model.objective_quadratic_elements_variables_2.push_back(0);
    model.objective_quadratic_elements_coefficients.push_back(1.0);

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

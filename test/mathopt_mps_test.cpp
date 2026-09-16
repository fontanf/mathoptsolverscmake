#include "mathoptsolverscmake/mathopt_mps.hpp"

#include <gtest/gtest.h>

#include <limits>
#include <sstream>

using namespace mathoptsolverscmake;

namespace
{

const double inf = std::numeric_limits<double>::infinity();

MathOptModel read(const std::string& mps)
{
    std::istringstream in(mps);
    return read_mps(in);
}

std::string write(const MathOptModel& model)
{
    std::ostringstream out;
    write_mps(model, out);
    return out.str();
}

}

TEST(MathOptMps, LinearRowTypesAndBounds)
{
    // maximize 3 x0 + 2 x1
    // s.t. x0 + x1 <= 4  (L)
    //      x0 + 2 x1 >= 1  (G)
    //      x0 - x1 = 0  (E)
    //      0 <= x0 <= 10
    //      x1 free
    MathOptModel model = read(
            "NAME          TEST\n"
            "OBJSENSE\n"
            " MAX\n"
            "ROWS\n"
            " N  obj\n"
            " L  c1\n"
            " G  c2\n"
            " E  c3\n"
            "COLUMNS\n"
            "    x0        obj              3   c1               1\n"
            "    x0        c2               1   c3               1\n"
            "    x1        obj              2   c1               1\n"
            "    x1        c2               2   c3              -1\n"
            "RHS\n"
            "    RHS       c1               4   c2               1\n"
            "BOUNDS\n"
            " UP BND       x0              10\n"
            " FR BND       x1\n"
            "ENDATA\n");

    EXPECT_EQ(model.objective_direction, ObjectiveDirection::Maximize);
    ASSERT_EQ(model.number_of_variables(), 2);
    ASSERT_EQ(model.number_of_constraints(), 3);

    EXPECT_EQ(model.variable_name(0), "x0");
    EXPECT_EQ(model.variable_name(1), "x1");
    EXPECT_EQ(model.variables_lower_bounds[0], 0.0);
    EXPECT_EQ(model.variables_upper_bounds[0], 10.0);
    EXPECT_EQ(model.variables_lower_bounds[1], -inf);
    EXPECT_EQ(model.variables_upper_bounds[1], inf);
    EXPECT_EQ(model.variables_types[0], VariableType::Continuous);

    EXPECT_EQ(model.objective_coefficients[0], 3.0);
    EXPECT_EQ(model.objective_coefficients[1], 2.0);

    EXPECT_EQ(model.constraint_name(0), "c1");
    EXPECT_EQ(model.constraints_lower_bounds[0], -inf);
    EXPECT_EQ(model.constraints_upper_bounds[0], 4.0);

    EXPECT_EQ(model.constraint_name(1), "c2");
    EXPECT_EQ(model.constraints_lower_bounds[1], 1.0);
    EXPECT_EQ(model.constraints_upper_bounds[1], inf);

    EXPECT_EQ(model.constraint_name(2), "c3");
    EXPECT_EQ(model.constraints_lower_bounds[2], 0.0);
    EXPECT_EQ(model.constraints_upper_bounds[2], 0.0);
}

TEST(MathOptMps, IntegerMarkerDefaultsToBinary)
{
    MathOptModel model = read(
            "ROWS\n"
            " N  obj\n"
            " L  c1\n"
            "COLUMNS\n"
            "    MARKER1   'MARKER'                 'INTORG'\n"
            "    x0        obj              1   c1               1\n"
            "    MARKER2   'MARKER'                 'INTEND'\n"
            "RHS\n"
            "    RHS       c1               5\n"
            "ENDATA\n");

    ASSERT_EQ(model.number_of_variables(), 1);
    // No explicit BOUNDS entry: an integer column defaults to [0, 1].
    EXPECT_EQ(model.variables_lower_bounds[0], 0.0);
    EXPECT_EQ(model.variables_upper_bounds[0], 1.0);
    EXPECT_EQ(model.variables_types[0], VariableType::Binary);
}

TEST(MathOptMps, IntegerMarkerWithExplicitBoundStaysInteger)
{
    MathOptModel model = read(
            "ROWS\n"
            " N  obj\n"
            " L  c1\n"
            "COLUMNS\n"
            "    MARKER1   'MARKER'                 'INTORG'\n"
            "    x0        obj              1   c1               1\n"
            "    MARKER2   'MARKER'                 'INTEND'\n"
            "RHS\n"
            "    RHS       c1               5\n"
            "BOUNDS\n"
            " UP BND       x0              20\n"
            "ENDATA\n");

    ASSERT_EQ(model.number_of_variables(), 1);
    // An explicit BOUNDS entry overrides the [0, 1] integer default.
    EXPECT_EQ(model.variables_lower_bounds[0], 0.0);
    EXPECT_EQ(model.variables_upper_bounds[0], 20.0);
    EXPECT_EQ(model.variables_types[0], VariableType::Integer);
}

TEST(MathOptMps, BoundTypes)
{
    MathOptModel model = read(
            "ROWS\n"
            " N  obj\n"
            "COLUMNS\n"
            "    x0        obj              1\n"
            "    x1        obj              1\n"
            "    x2        obj              1\n"
            "    x3        obj              1\n"
            "    x4        obj              1\n"
            "BOUNDS\n"
            " LO BND       x0              -5\n"
            " FX BND       x1               2\n"
            " MI BND       x2\n"
            " BV BND       x3\n"
            " UI BND       x4              10\n"
            "ENDATA\n");

    ASSERT_EQ(model.number_of_variables(), 5);

    // LO: only the lower bound is set, upper stays at its +inf default.
    EXPECT_EQ(model.variables_lower_bounds[0], -5.0);
    EXPECT_EQ(model.variables_upper_bounds[0], inf);

    // FX: fixed to a single value.
    EXPECT_EQ(model.variables_lower_bounds[1], 2.0);
    EXPECT_EQ(model.variables_upper_bounds[1], 2.0);

    // MI: -inf lower bound, upper stays at its +inf default.
    EXPECT_EQ(model.variables_lower_bounds[2], -inf);
    EXPECT_EQ(model.variables_upper_bounds[2], inf);

    // BV: binary, regardless of any integer marker.
    EXPECT_EQ(model.variables_lower_bounds[3], 0.0);
    EXPECT_EQ(model.variables_upper_bounds[3], 1.0);
    EXPECT_EQ(model.variables_types[3], VariableType::Binary);

    // UI: integer upper bound, and the column becomes integer.
    EXPECT_EQ(model.variables_upper_bounds[4], 10.0);
    EXPECT_EQ(model.variables_types[4], VariableType::Integer);
}

TEST(MathOptMps, Ranges)
{
    MathOptModel model = read(
            "ROWS\n"
            " N  obj\n"
            " L  c1\n"
            " G  c2\n"
            " E  c3\n"
            " E  c4\n"
            "COLUMNS\n"
            "    x0        obj              1   c1               1\n"
            "    x0        c2               1   c3               1\n"
            "    x0        c4               1\n"
            "RHS\n"
            "    RHS       c1              10   c2               2\n"
            "    RHS       c3               5   c4               5\n"
            "RANGES\n"
            "    RNG       c1               3   c2               4\n"
            "    RNG       c3               2   c4              -2\n"
            "ENDATA\n");

    ASSERT_EQ(model.number_of_constraints(), 4);

    // L row: [upper - |r|, upper].
    EXPECT_EQ(model.constraints_lower_bounds[0], 7.0);
    EXPECT_EQ(model.constraints_upper_bounds[0], 10.0);

    // G row: [lower, lower + |r|].
    EXPECT_EQ(model.constraints_lower_bounds[1], 2.0);
    EXPECT_EQ(model.constraints_upper_bounds[1], 6.0);

    // E row, r >= 0: [rhs, rhs + |r|].
    EXPECT_EQ(model.constraints_lower_bounds[2], 5.0);
    EXPECT_EQ(model.constraints_upper_bounds[2], 7.0);

    // E row, r < 0: [rhs - |r|, rhs].
    EXPECT_EQ(model.constraints_lower_bounds[3], 3.0);
    EXPECT_EQ(model.constraints_upper_bounds[3], 5.0);
}

TEST(MathOptMps, FreeRowIsDropped)
{
    MathOptModel model = read(
            "ROWS\n"
            " N  obj\n"
            " N  ignored\n"
            " L  c1\n"
            "COLUMNS\n"
            "    x0        obj              1   ignored          1\n"
            "    x0        c1               1\n"
            "RHS\n"
            "    RHS       c1               1\n"
            "ENDATA\n");

    // The second N row is a free row: dropped entirely, not just left empty.
    ASSERT_EQ(model.number_of_constraints(), 1);
    EXPECT_EQ(model.constraint_name(0), "c1");
}

TEST(MathOptMps, UnsupportedSectionThrows)
{
    EXPECT_THROW(
            read(
                    "ROWS\n"
                    " N  obj\n"
                    "COLUMNS\n"
                    "    x0        obj              1\n"
                    "SOS\n"
                    " S1 SET               sos1\n"
                    "    x0                              1\n"
                    "ENDATA\n"),
            std::invalid_argument);
}

TEST(MathOptMps, MissingFileThrows)
{
    EXPECT_THROW(read_mps("/no/such/file.mps"), std::invalid_argument);
}

TEST(MathOptMps, WriteThenReadRoundTrip)
{
    MathOptModel original = read(
            "NAME          TEST\n"
            "OBJSENSE\n"
            " MAX\n"
            "ROWS\n"
            " N  obj\n"
            " L  c1\n"
            " G  c2\n"
            " E  c3\n"
            "COLUMNS\n"
            "    x0        obj              3   c1               1\n"
            "    x0        c2               1   c3               1\n"
            "    x1        obj              2   c1               1\n"
            "    x1        c2               2   c3              -1\n"
            "RHS\n"
            "    RHS       c1               4   c2               1\n"
            "BOUNDS\n"
            " UP BND       x0              10\n"
            " FR BND       x1\n"
            "ENDATA\n");

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

TEST(MathOptMps, WriteRangeConstraintRoundTrip)
{
    MathOptModel model(1, 1, 1);
    model.objective_direction = ObjectiveDirection::Minimize;
    model.variables_lower_bounds[0] = 0.0;
    model.variables_upper_bounds[0] = 100.0;
    model.objective_coefficients[0] = 1.0;
    model.constraints_starts[0] = 0;
    model.elements_variables[0] = 0;
    model.elements_coefficients[0] = 1.0;
    model.constraints_lower_bounds[0] = 3.0;
    model.constraints_upper_bounds[0] = 8.0;
    ASSERT_EQ(model.constraint_sense(0), ConstraintSense::Range);

    MathOptModel round_tripped = read(write(model));

    ASSERT_EQ(round_tripped.number_of_constraints(), 1);
    EXPECT_EQ(round_tripped.constraints_lower_bounds[0], 3.0);
    EXPECT_EQ(round_tripped.constraints_upper_bounds[0], 8.0);
}

TEST(MathOptMps, WriteOrphanVariableRoundTrip)
{
    // x1 has a zero objective coefficient and appears in no constraint; it
    // must still come back as a variable of the round-tripped model.
    MathOptModel model(2, 1, 1);
    model.objective_direction = ObjectiveDirection::Minimize;
    model.objective_coefficients[0] = 1.0;
    model.objective_coefficients[1] = 0.0;
    model.constraints_starts[0] = 0;
    model.elements_variables[0] = 0;
    model.elements_coefficients[0] = 1.0;
    model.constraints_lower_bounds[0] = -inf;
    model.constraints_upper_bounds[0] = 10.0;

    MathOptModel round_tripped = read(write(model));

    EXPECT_EQ(round_tripped.number_of_variables(), 2);
}

TEST(MathOptMps, WriteIntegerAndBinaryRoundTrip)
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
    model.variables_lower_bounds[2] = -5.0;
    model.variables_upper_bounds[2] = 5.0;
    model.variables_types[2] = VariableType::Integer;

    MathOptModel round_tripped = read(write(model));

    ASSERT_EQ(round_tripped.number_of_variables(), 3);
    EXPECT_EQ(round_tripped.variables_types[0], VariableType::Binary);
    EXPECT_EQ(round_tripped.variables_types[1], VariableType::Integer);
    EXPECT_EQ(round_tripped.variables_upper_bounds[1], 20.0);
    EXPECT_EQ(round_tripped.variables_types[2], VariableType::Integer);
    EXPECT_EQ(round_tripped.variables_lower_bounds[2], -5.0);
    EXPECT_EQ(round_tripped.variables_upper_bounds[2], 5.0);
}

TEST(MathOptMps, WriteNonMilpThrows)
{
    MathOptModel model(1, 0, 0);
    model.objective_direction = ObjectiveDirection::Minimize;
    model.objective_quadratic_elements_variables_1.push_back(0);
    model.objective_quadratic_elements_variables_2.push_back(0);
    model.objective_quadratic_elements_coefficients.push_back(1.0);

    std::ostringstream out;
    EXPECT_THROW(write_mps(model, out), std::invalid_argument);
}

TEST(MathOptMps, WriteWhitespaceNameThrows)
{
    MathOptModel model(1, 0, 0);
    model.objective_direction = ObjectiveDirection::Minimize;
    model.objective_coefficients[0] = 1.0;
    model.variables_names = {"bad name"};

    std::ostringstream out;
    EXPECT_THROW(write_mps(model, out), std::invalid_argument);
}

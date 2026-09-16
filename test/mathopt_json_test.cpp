#include "mathoptsolverscmake/mathopt_json.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <sstream>

using namespace mathoptsolverscmake;

namespace
{

const double inf = std::numeric_limits<double>::infinity();

MathOptModel read(const std::string& json)
{
    std::istringstream in(json);
    return read_json(in);
}

std::string write(const MathOptModel& model)
{
    std::ostringstream out;
    write_json(model, out);
    return out.str();
}

}

TEST(MathOptJson, MinimalRequiredFieldsDefault)
{
    MathOptModel model = read("{\"objective_direction\": \"Minimize\"}");

    EXPECT_EQ(model.objective_direction, ObjectiveDirection::Minimize);
    EXPECT_EQ(model.number_of_variables(), 0);
    EXPECT_EQ(model.number_of_constraints(), 0);
    EXPECT_EQ(model.feasibility_tolerance, 0.0);
    EXPECT_EQ(model.integrality_tolerance, 0.0);
}

TEST(MathOptJson, VariablesTypesDefaultsToContinuous)
{
    MathOptModel model = read(
            "{\n"
            "  \"objective_direction\": \"Maximize\",\n"
            "  \"variables_lower_bounds\": [0, 0, 0],\n"
            "  \"variables_upper_bounds\": [1, 1, 1]\n"
            "}\n");

    ASSERT_EQ(model.number_of_variables(), 3);
    for (int variable_id = 0; variable_id < 3; ++variable_id)
        EXPECT_EQ(model.variables_types[variable_id], VariableType::Continuous);
}

TEST(MathOptJson, InfinityAndNaNTokens)
{
    MathOptModel model = read(
            "{\n"
            "  \"objective_direction\": \"Minimize\",\n"
            "  \"variables_lower_bounds\": [-Infinity, 0],\n"
            "  \"variables_upper_bounds\": [Infinity, NaN]\n"
            "}\n");

    ASSERT_EQ(model.number_of_variables(), 2);
    EXPECT_EQ(model.variables_lower_bounds[0], -inf);
    EXPECT_EQ(model.variables_upper_bounds[0], inf);
    EXPECT_TRUE(std::isnan(model.variables_upper_bounds[1]));
}

TEST(MathOptJson, MissingObjectiveDirectionThrows)
{
    EXPECT_THROW(read("{\"variables_lower_bounds\": [0]}"), std::invalid_argument);
}

TEST(MathOptJson, WrongFieldTypeThrows)
{
    EXPECT_THROW(
            read("{\"objective_direction\": \"Minimize\", \"variables_lower_bounds\": \"not an array\"}"),
            std::invalid_argument);
}

TEST(MathOptJson, MalformedJsonThrows)
{
    EXPECT_THROW(read("{\"objective_direction\": "), std::invalid_argument);
}

TEST(MathOptJson, MissingFileThrows)
{
    EXPECT_THROW(read_json("/no/such/file.json"), std::invalid_argument);
}

TEST(MathOptJson, WriteThenReadRoundTrip)
{
    // maximize 3 x0 + 2 x1 + x1 * x1
    // s.t. x0 + x1 <= 4  (linear)
    //      2 <= x0 - x1 <= 6  (range)
    // x0 in [0, 10], x1 integer in [-5, 5]
    MathOptModel original(2, 2, 2);
    original.objective_direction = ObjectiveDirection::Maximize;

    original.variables_lower_bounds = {0.0, -5.0};
    original.variables_upper_bounds = {10.0, 5.0};
    original.variables_types = {VariableType::Continuous, VariableType::Integer};
    original.variables_names = {"x0", "x1"};
    original.variables_initial_values = {1.0, -1.0};

    original.objective_coefficients = {3.0, 2.0};
    original.objective_quadratic_elements_variables_1 = {1};
    original.objective_quadratic_elements_variables_2 = {1};
    original.objective_quadratic_elements_coefficients = {1.0};

    original.constraints_lower_bounds = {-inf, 2.0};
    original.constraints_upper_bounds = {4.0, 6.0};
    original.constraints_names = {"c0", "c1"};

    original.constraints_starts = {0, 2};
    original.elements_variables = {0, 1, 0, 1};
    original.elements_coefficients = {1.0, 1.0, 1.0, -1.0};

    original.feasibility_tolerance = 1e-6;
    original.integrality_tolerance = 1e-5;

    MathOptModel round_tripped = read(write(original));

    ASSERT_EQ(round_tripped.number_of_variables(), original.number_of_variables());
    ASSERT_EQ(round_tripped.number_of_constraints(), original.number_of_constraints());
    EXPECT_EQ(round_tripped.objective_direction, original.objective_direction);
    EXPECT_EQ(round_tripped.variables_lower_bounds, original.variables_lower_bounds);
    EXPECT_EQ(round_tripped.variables_upper_bounds, original.variables_upper_bounds);
    EXPECT_EQ(round_tripped.variables_types, original.variables_types);
    EXPECT_EQ(round_tripped.variables_names, original.variables_names);
    EXPECT_EQ(round_tripped.variables_initial_values, original.variables_initial_values);
    EXPECT_EQ(round_tripped.objective_coefficients, original.objective_coefficients);
    EXPECT_EQ(round_tripped.objective_quadratic_elements_variables_1, original.objective_quadratic_elements_variables_1);
    EXPECT_EQ(round_tripped.objective_quadratic_elements_variables_2, original.objective_quadratic_elements_variables_2);
    EXPECT_EQ(round_tripped.objective_quadratic_elements_coefficients, original.objective_quadratic_elements_coefficients);
    EXPECT_EQ(round_tripped.constraints_lower_bounds, original.constraints_lower_bounds);
    EXPECT_EQ(round_tripped.constraints_upper_bounds, original.constraints_upper_bounds);
    EXPECT_EQ(round_tripped.constraints_names, original.constraints_names);
    EXPECT_EQ(round_tripped.constraints_starts, original.constraints_starts);
    EXPECT_EQ(round_tripped.elements_variables, original.elements_variables);
    EXPECT_EQ(round_tripped.elements_coefficients, original.elements_coefficients);
    EXPECT_EQ(round_tripped.feasibility_tolerance, original.feasibility_tolerance);
    EXPECT_EQ(round_tripped.integrality_tolerance, original.integrality_tolerance);
}

TEST(MathOptJson, NonlinearAndQuadraticConstraintRoundTrip)
{
    // A single constraint with a quadratic term (x0 * x1) and a nonlinear
    // term (x0 + 2), stored as a DFS pre-order tree: '+' (2 children),
    // 'v' (variable 0), 'k' (constant 2).
    MathOptModel original(2, 1, 0);
    original.objective_direction = ObjectiveDirection::Minimize;
    original.objective_coefficients = {0.0, 0.0};

    original.constraints_lower_bounds = {-inf};
    original.constraints_upper_bounds = {10.0};
    original.constraints_starts = {0};

    original.quadratic_elements_constraints_starts = {0};
    original.quadratic_elements_variables_1 = {0};
    original.quadratic_elements_variables_2 = {1};
    original.quadratic_elements_coefficients = {1.0};

    original.nonlinear_elements_constraints_starts = {0};
    original.nonlinear_elements_operators = {'+', 'v', 'k'};
    original.nonlinear_elements_values = {0.0, 0.0, 2.0};
    original.nonlinear_elements_variables = {-1, 0, -1};
    original.nonlinear_elements_number_of_children = {2, 0, 0};

    MathOptModel round_tripped = read(write(original));

    EXPECT_EQ(round_tripped.quadratic_elements_constraints_starts, original.quadratic_elements_constraints_starts);
    EXPECT_EQ(round_tripped.quadratic_elements_variables_1, original.quadratic_elements_variables_1);
    EXPECT_EQ(round_tripped.quadratic_elements_variables_2, original.quadratic_elements_variables_2);
    EXPECT_EQ(round_tripped.quadratic_elements_coefficients, original.quadratic_elements_coefficients);
    EXPECT_EQ(round_tripped.nonlinear_elements_constraints_starts, original.nonlinear_elements_constraints_starts);
    EXPECT_EQ(round_tripped.nonlinear_elements_operators, original.nonlinear_elements_operators);
    EXPECT_EQ(round_tripped.nonlinear_elements_values, original.nonlinear_elements_values);
    EXPECT_EQ(round_tripped.nonlinear_elements_variables, original.nonlinear_elements_variables);
    EXPECT_EQ(round_tripped.nonlinear_elements_number_of_children, original.nonlinear_elements_number_of_children);

    // Both structures are present, so the model is neither pure MILP nor solvable by MPS.
    EXPECT_TRUE(original.has_quadratic());
    EXPECT_TRUE(original.has_nonlinear());
}

TEST(MathOptJson, EmptyOptionalFieldsOmittedFromOutput)
{
    MathOptModel model(1, 0, 0);
    model.objective_direction = ObjectiveDirection::Minimize;
    model.objective_coefficients[0] = 1.0;

    std::string json = write(model);
    EXPECT_EQ(json.find("variables_names"), std::string::npos);
    EXPECT_EQ(json.find("quadratic"), std::string::npos);
    EXPECT_EQ(json.find("nonlinear"), std::string::npos);
}

TEST(MathOptJson, BlackBoxObjectiveThrows)
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
    EXPECT_THROW(write_json(model, out), std::invalid_argument);
}

#include "mathoptsolverscmake/mathopt_statistics.hpp"

#include "mathoptsolverscmake/mathopt_json.hpp"
#include "mathoptsolverscmake/mathopt_lp.hpp"
#include "mathoptsolverscmake/mathopt_mps.hpp"
#include "mathoptsolverscmake/mathopt_nl.hpp"

#include <gtest/gtest.h>

#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>

using namespace mathoptsolverscmake;

namespace
{

const double inf = std::numeric_limits<double>::infinity();

/**
 * A linear model with a continuous, a binary and an integer variable, and two
 * constraints with 3 nonzeros.
 */
MathOptModel linear_model()
{
    MathOptModel model;
    model.objective_direction = ObjectiveDirection::Maximize;
    model.variables_lower_bounds = {0, 0, 0};
    model.variables_upper_bounds = {10, 1, 5};
    model.variables_types = {VariableType::Continuous, VariableType::Binary, VariableType::Integer};
    model.objective_coefficients = {1, 2, 3};
    model.constraints_lower_bounds = {-inf, 1};
    model.constraints_upper_bounds = {10, inf};
    model.constraints_starts = {0, 2};
    model.elements_variables = {0, 1, 2};
    model.elements_coefficients = {1, 1, 1};
    return model;
}

void expect_equal(const MathOptModelStatistics& a, const MathOptModelStatistics& b, const std::string& what)
{
    EXPECT_EQ(a.number_of_variables, b.number_of_variables) << what;
    EXPECT_EQ(a.number_of_continuous_variables, b.number_of_continuous_variables) << what;
    EXPECT_EQ(a.number_of_binary_variables, b.number_of_binary_variables) << what;
    EXPECT_EQ(a.number_of_integer_variables, b.number_of_integer_variables) << what;
    EXPECT_EQ(a.number_of_constraints, b.number_of_constraints) << what;
    EXPECT_EQ(a.number_of_elements, b.number_of_elements) << what;
    EXPECT_EQ(a.number_of_nonlinear_constraints, b.number_of_nonlinear_constraints) << what;
    EXPECT_EQ(a.number_of_objectives, b.number_of_objectives) << what;
}

}

TEST(MathOptStatistics, LinearModel)
{
    MathOptModelStatistics statistics = compute_statistics(linear_model());
    EXPECT_EQ(statistics.number_of_variables, 3);
    EXPECT_EQ(statistics.number_of_continuous_variables, 1);
    EXPECT_EQ(statistics.number_of_binary_variables, 1);
    EXPECT_EQ(statistics.number_of_integer_variables, 1);
    EXPECT_EQ(statistics.number_of_constraints, 2);
    EXPECT_EQ(statistics.number_of_elements, 3);
    EXPECT_EQ(statistics.number_of_quadratic_constraints, 0);
    EXPECT_EQ(statistics.number_of_nonlinear_constraints, 0);
    EXPECT_EQ(statistics.number_of_objectives, 1);
    EXPECT_EQ(statistics.number_of_quadratic_objectives, 0);
    EXPECT_EQ(statistics.number_of_nonlinear_objectives, 0);
}

TEST(MathOptStatistics, VariableTypes)
{
    MathOptModel model;
    model.variables_lower_bounds = {0, 0, 0, 0, 1, 1};
    model.variables_upper_bounds = {1, 1, 1, 5, 4, 4};
    model.variables_types = {VariableType::Continuous, VariableType::Binary, VariableType::Integer,
        VariableType::Integer, VariableType::SemiContinuous, VariableType::SemiInteger};
    model.objective_coefficients = {0, 0, 0, 0, 0, 0};
    MathOptModelStatistics statistics = compute_statistics(model);
    EXPECT_EQ(statistics.number_of_continuous_variables, 1);
    // The Integer variable bounded to [0, 1] counts as binary.
    EXPECT_EQ(statistics.number_of_binary_variables, 2);
    EXPECT_EQ(statistics.number_of_integer_variables, 1);
    EXPECT_EQ(statistics.number_of_semi_continuous_variables, 1);
    EXPECT_EQ(statistics.number_of_semi_integer_variables, 1);
    // No objective term.
    EXPECT_EQ(statistics.number_of_objectives, 0);
}

TEST(MathOptStatistics, QuadraticAndNonlinear)
{
    MathOptModel model = linear_model();
    // Constraint 0: + x0 * x1; constraint 1: exp(x2).
    model.quadratic_elements_constraints_starts = {0, 1};
    model.quadratic_elements_variables_1 = {0};
    model.quadratic_elements_variables_2 = {1};
    model.quadratic_elements_coefficients = {1};
    model.nonlinear_elements_constraints_starts = {0, 0};
    model.nonlinear_elements_operators = {'e', 'v'};
    model.nonlinear_elements_values = {0, 0};
    model.nonlinear_elements_variables = {-1, 2};
    model.nonlinear_elements_number_of_children = {1, 0};
    model.objective_quadratic_elements_variables_1 = {0};
    model.objective_quadratic_elements_variables_2 = {0};
    model.objective_quadratic_elements_coefficients = {1};
    MathOptModelStatistics statistics = compute_statistics(model);
    EXPECT_EQ(statistics.number_of_quadratic_constraints, 1);
    EXPECT_EQ(statistics.number_of_nonlinear_constraints, 1);
    EXPECT_EQ(statistics.number_of_quadratic_objectives, 1);
    EXPECT_EQ(statistics.number_of_nonlinear_objectives, 0);
}

TEST(MathOptStatistics, ReadStatisticsOfEachFormat)
{
    MathOptModel model = linear_model();
    MathOptModelStatistics expected = compute_statistics(model);
    std::filesystem::path directory = std::filesystem::temp_directory_path() / "mathopt_statistics_test";
    std::filesystem::create_directories(directory);
    write_mps(model, (directory / "model.mps").string());
    write_lp(model, (directory / "model.LP").string());
    write_json(model, (directory / "model.json").string());
    write_nl(model, (directory / "model.nl").string());
    for (const char* name: {"model.mps", "model.LP", "model.json", "model.nl"})
        expect_equal(read_statistics((directory / name).string()), expected, name);
    EXPECT_THROW(read_statistics((directory / "model.txt").string()), std::invalid_argument);
    std::filesystem::remove_all(directory);
}

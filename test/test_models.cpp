#include "test_models.hpp"

#include "mathoptsolverscmake/mathopt_nonlinear.hpp"

#include <cmath>
#include <vector>

using namespace mathoptsolverscmake;
using mathoptsolverscmaketest::TestModel;

std::ostream& mathoptsolverscmaketest::operator<<(
        std::ostream& os,
        const TestModel& test_model)
{
    return os << test_model.name;
}

TestModel mathoptsolverscmaketest::linear_continuous_model()
{
    // maximize 3 x0 + 2 x1
    // s.t.     x0 +   x1 <= 4
    //          x0 + 3 x1 <= 6
    //          x0, x1 >= 0
    // Optimal: x0 = 4, x1 = 0, objective = 12.
    MathOptModel model(2, 2, 4);
    model.objective_direction = ObjectiveDirection::Maximize;

    model.variables_lower_bounds[0] = 0.0;
    model.variables_lower_bounds[1] = 0.0;

    model.objective_coefficients[0] = 3.0;
    model.objective_coefficients[1] = 2.0;

    model.constraints_starts[0] = 0;
    model.elements_variables[0] = 0;
    model.elements_coefficients[0] = 1.0;
    model.elements_variables[1] = 1;
    model.elements_coefficients[1] = 1.0;
    model.constraints_upper_bounds[0] = 4.0;

    model.constraints_starts[1] = 2;
    model.elements_variables[2] = 0;
    model.elements_coefficients[2] = 1.0;
    model.elements_variables[3] = 1;
    model.elements_coefficients[3] = 3.0;
    model.constraints_upper_bounds[1] = 6.0;

    return TestModel{"linear_continuous", model, {4.0, 0.0}};
}

TestModel mathoptsolverscmaketest::linear_mixed_integer_model()
{
    // 0/1 knapsack: values [3, 4, 5], weights [2, 3, 4], capacity 5.
    // Optimal: pick items 0 and 1, x = [1, 1, 0], objective = 7.
    MathOptModel model(3, 1, 3);
    model.objective_direction = ObjectiveDirection::Maximize;

    for (int variable_id = 0; variable_id < 3; ++variable_id) {
        model.variables_lower_bounds[variable_id] = 0.0;
        model.variables_upper_bounds[variable_id] = 1.0;
        model.variables_types[variable_id] = VariableType::Binary;
    }
    model.objective_coefficients[0] = 3.0;
    model.objective_coefficients[1] = 4.0;
    model.objective_coefficients[2] = 5.0;

    model.constraints_starts[0] = 0;
    model.elements_variables[0] = 0;
    model.elements_coefficients[0] = 2.0;
    model.elements_variables[1] = 1;
    model.elements_coefficients[1] = 3.0;
    model.elements_variables[2] = 2;
    model.elements_coefficients[2] = 4.0;
    model.constraints_upper_bounds[0] = 5.0;

    return TestModel{"linear_mixed_integer", model, {1.0, 1.0, 0.0}};
}

TestModel mathoptsolverscmaketest::linear_infeasible_model()
{
    // minimize x0
    // s.t.     x0 <= 1
    //          x0 >= 5
    // Infeasible: no x0 satisfies both constraints.
    MathOptModel model(1, 2, 2);
    model.objective_direction = ObjectiveDirection::Minimize;

    model.objective_coefficients[0] = 1.0;

    model.constraints_starts[0] = 0;
    model.elements_variables[0] = 0;
    model.elements_coefficients[0] = 1.0;
    model.constraints_upper_bounds[0] = 1.0;

    model.constraints_starts[1] = 1;
    model.elements_variables[1] = 0;
    model.elements_coefficients[1] = 1.0;
    model.constraints_lower_bounds[1] = 5.0;

    return TestModel{"linear_infeasible", model, {}, false};
}

TestModel mathoptsolverscmaketest::linear_mixed_integer_infeasible_model()
{
    // minimize x0, x0 integer
    // s.t.     2 x0 = 1
    // The continuous relaxation is feasible (x0 = 0.5), but no integer x0
    // satisfies 2 x0 = 1.
    MathOptModel model(1, 1, 1);
    model.objective_direction = ObjectiveDirection::Minimize;

    model.variables_lower_bounds[0] = 0.0;
    model.variables_upper_bounds[0] = 10.0;
    model.variables_types[0] = VariableType::Integer;

    model.objective_coefficients[0] = 1.0;

    model.constraints_starts[0] = 0;
    model.elements_variables[0] = 0;
    model.elements_coefficients[0] = 2.0;
    model.constraints_lower_bounds[0] = 1.0;
    model.constraints_upper_bounds[0] = 1.0;

    return TestModel{"linear_mixed_integer_infeasible", model, {}, false};
}

TestModel mathoptsolverscmaketest::quadratic_continuous_model()
{
    // minimize x0^2 - 2 x0 + x1^2 - 4 x1
    // s.t.     0 <= x0, x1 <= 10
    // Optimal: x0 = 1, x1 = 2, objective = -5.
    MathOptModel model(2, 0, 0);
    model.objective_direction = ObjectiveDirection::Minimize;

    model.variables_lower_bounds[0] = 0.0;
    model.variables_upper_bounds[0] = 10.0;
    model.variables_lower_bounds[1] = 0.0;
    model.variables_upper_bounds[1] = 10.0;

    model.objective_coefficients[0] = -2.0;
    model.objective_coefficients[1] = -4.0;
    model.objective_quadratic_elements_variables_1 = {0, 1};
    model.objective_quadratic_elements_variables_2 = {0, 1};
    model.objective_quadratic_elements_coefficients = {1.0, 1.0};

    return TestModel{"quadratic_continuous", model, {1.0, 2.0}};
}

TestModel mathoptsolverscmaketest::quadratic_mixed_integer_model()
{
    // minimize x0^2 - 2.4 x0 + x1^2 - 4 x1, x0 integer
    // s.t.     0 <= x0, x1 <= 10
    // Unconstrained real optimum of the x0 term is at 1.2; the integer
    // optimum is x0 = 1 (f(1) = -1.4 < f(2) = -0.8).
    // Optimal: x0 = 1, x1 = 2, objective = -5.4.
    MathOptModel model(2, 0, 0);
    model.objective_direction = ObjectiveDirection::Minimize;

    model.variables_lower_bounds[0] = 0.0;
    model.variables_upper_bounds[0] = 10.0;
    model.variables_types[0] = VariableType::Integer;
    model.variables_lower_bounds[1] = 0.0;
    model.variables_upper_bounds[1] = 10.0;

    model.objective_coefficients[0] = -2.4;
    model.objective_coefficients[1] = -4.0;
    model.objective_quadratic_elements_variables_1 = {0, 1};
    model.objective_quadratic_elements_variables_2 = {0, 1};
    model.objective_quadratic_elements_coefficients = {1.0, 1.0};

    return TestModel{"quadratic_mixed_integer", model, {1.0, 2.0}};
}

TestModel mathoptsolverscmaketest::quadratic_infeasible_model()
{
    // minimize x0^2 - 2 x0 + x1^2 - 4 x1
    // s.t.     x0 <= 1
    //          x0 >= 5
    // Infeasible regardless of the (otherwise unconstrained) quadratic
    // objective.
    MathOptModel model(2, 2, 2);
    model.objective_direction = ObjectiveDirection::Minimize;

    model.objective_coefficients[0] = -2.0;
    model.objective_coefficients[1] = -4.0;
    model.objective_quadratic_elements_variables_1 = {0, 1};
    model.objective_quadratic_elements_variables_2 = {0, 1};
    model.objective_quadratic_elements_coefficients = {1.0, 1.0};

    model.constraints_starts[0] = 0;
    model.elements_variables[0] = 0;
    model.elements_coefficients[0] = 1.0;
    model.constraints_upper_bounds[0] = 1.0;

    model.constraints_starts[1] = 1;
    model.elements_variables[1] = 0;
    model.elements_coefficients[1] = 1.0;
    model.constraints_lower_bounds[1] = 5.0;

    return TestModel{"quadratic_infeasible", model, {}, false};
}

TestModel mathoptsolverscmaketest::quadratic_mixed_integer_infeasible_model()
{
    // minimize x0^2 - 2.4 x0 + x1^2 - 4 x1, x0 integer
    // s.t.     2 x0 = 1
    // The continuous relaxation is feasible (x0 = 0.5), but no integer x0
    // satisfies 2 x0 = 1.
    MathOptModel model(2, 1, 1);
    model.objective_direction = ObjectiveDirection::Minimize;

    model.variables_lower_bounds[0] = 0.0;
    model.variables_upper_bounds[0] = 10.0;
    model.variables_types[0] = VariableType::Integer;
    model.variables_lower_bounds[1] = 0.0;
    model.variables_upper_bounds[1] = 10.0;

    model.objective_coefficients[0] = -2.4;
    model.objective_coefficients[1] = -4.0;
    model.objective_quadratic_elements_variables_1 = {0, 1};
    model.objective_quadratic_elements_variables_2 = {0, 1};
    model.objective_quadratic_elements_coefficients = {1.0, 1.0};

    model.constraints_starts[0] = 0;
    model.elements_variables[0] = 0;
    model.elements_coefficients[0] = 2.0;
    model.constraints_lower_bounds[0] = 1.0;
    model.constraints_upper_bounds[0] = 1.0;

    return TestModel{"quadratic_mixed_integer_infeasible", model, {}, false};
}

TestModel mathoptsolverscmaketest::nonlinear_continuous_model()
{
    // minimize (x0 - 1)^2 + exp(x1) - x1
    // s.t.     exp(x0) <= 2
    //          -10 <= x0 <= 10, -5 <= x1 <= 5
    // The unconstrained optimum x0 = 1 is cut off by the constraint, so the
    // constrained optimum sits on the boundary: x0 = ln(2).
    // Optimal: x0 = ln(2), x1 = 0, objective = (ln(2) - 1)^2 + 1.
    MathOptModel model(2, 1, 0);
    model.objective_direction = ObjectiveDirection::Minimize;

    model.variables_lower_bounds[0] = -10.0;
    model.variables_upper_bounds[0] = 10.0;
    model.variables_lower_bounds[1] = -5.0;
    model.variables_upper_bounds[1] = 5.0;

    // Objective: pow(x0 - 1, 2) + (exp(x1) - x1).
    set_objective_nonlinear_expression(
            model,
            nonlinear_sum({
                nonlinear_power(
                        nonlinear_subtract(nonlinear_variable(0), nonlinear_constant(1.0)),
                        nonlinear_constant(2.0)),
                nonlinear_subtract(nonlinear_exp(nonlinear_variable(1)), nonlinear_variable(1)),
            }));

    // Constraint: exp(x0) <= 2.
    model.constraints_upper_bounds[0] = 2.0;
    add_constraint_nonlinear_expression(model, 0, nonlinear_exp(nonlinear_variable(0)));

    double x0 = std::log(2.0);
    return TestModel{"nonlinear_continuous", model, {x0, 0.0}};
}

TestModel mathoptsolverscmaketest::nonlinear_mixed_integer_model()
{
    // minimize (x0 - 1)^2 + exp(x1) - x1, x0 integer
    // s.t.     exp(x0) <= 1.5
    //          -10 <= x0 <= 10, -5 <= x1 <= 5
    // The constraint rules out x0 = 1 (exp(1) > 1.5); among the remaining
    // integers x0 = 0 is the closest to the unconstrained optimum.
    // Optimal: x0 = 0, x1 = 0, objective = 2.
    MathOptModel model(2, 1, 0);
    model.objective_direction = ObjectiveDirection::Minimize;

    model.variables_lower_bounds[0] = -10.0;
    model.variables_upper_bounds[0] = 10.0;
    model.variables_types[0] = VariableType::Integer;
    model.variables_lower_bounds[1] = -5.0;
    model.variables_upper_bounds[1] = 5.0;

    set_objective_nonlinear_expression(
            model,
            nonlinear_sum({
                nonlinear_power(
                        nonlinear_subtract(nonlinear_variable(0), nonlinear_constant(1.0)),
                        nonlinear_constant(2.0)),
                nonlinear_subtract(nonlinear_exp(nonlinear_variable(1)), nonlinear_variable(1)),
            }));

    model.constraints_upper_bounds[0] = 1.5;
    add_constraint_nonlinear_expression(model, 0, nonlinear_exp(nonlinear_variable(0)));

    return TestModel{"nonlinear_mixed_integer", model, {0.0, 0.0}};
}

TestModel mathoptsolverscmaketest::nonlinear_infeasible_model()
{
    // minimize (x0 - 1)^2 + exp(x1) - x1
    // s.t.     exp(x0) <= 0.5
    //          0 <= x0 <= 10, -5 <= x1 <= 5
    // Infeasible: exp(x0) >= 1 for every x0 >= 0, but the constraint
    // requires exp(x0) <= 0.5.
    MathOptModel model(2, 1, 0);
    model.objective_direction = ObjectiveDirection::Minimize;

    model.variables_lower_bounds[0] = 0.0;
    model.variables_upper_bounds[0] = 10.0;
    model.variables_lower_bounds[1] = -5.0;
    model.variables_upper_bounds[1] = 5.0;

    set_objective_nonlinear_expression(
            model,
            nonlinear_sum({
                nonlinear_power(
                        nonlinear_subtract(nonlinear_variable(0), nonlinear_constant(1.0)),
                        nonlinear_constant(2.0)),
                nonlinear_subtract(nonlinear_exp(nonlinear_variable(1)), nonlinear_variable(1)),
            }));

    model.constraints_upper_bounds[0] = 0.5;
    add_constraint_nonlinear_expression(model, 0, nonlinear_exp(nonlinear_variable(0)));

    return TestModel{"nonlinear_infeasible", model, {}, false};
}

TestModel mathoptsolverscmaketest::nonlinear_mixed_integer_infeasible_model()
{
    // minimize (x0 - 1)^2 + exp(x1) - x1, x0 integer
    // s.t.     exp(x0) <= 1.5   (nonlinear constraint 0)
    //          2 x0 = 1         (linear constraint 1)
    // The continuous relaxation is feasible (x0 = 0.5 satisfies both), but
    // no integer x0 satisfies 2 x0 = 1.
    MathOptModel model(2, 2, 1);
    model.objective_direction = ObjectiveDirection::Minimize;

    model.variables_lower_bounds[0] = -10.0;
    model.variables_upper_bounds[0] = 10.0;
    model.variables_types[0] = VariableType::Integer;
    model.variables_lower_bounds[1] = -5.0;
    model.variables_upper_bounds[1] = 5.0;

    set_objective_nonlinear_expression(
            model,
            nonlinear_sum({
                nonlinear_power(
                        nonlinear_subtract(nonlinear_variable(0), nonlinear_constant(1.0)),
                        nonlinear_constant(2.0)),
                nonlinear_subtract(nonlinear_exp(nonlinear_variable(1)), nonlinear_variable(1)),
            }));

    // Constraint 0 (nonlinear): exp(x0) <= 1.5.
    model.constraints_upper_bounds[0] = 1.5;
    add_constraint_nonlinear_expression(model, 0, nonlinear_exp(nonlinear_variable(0)));
    model.nonlinear_elements_constraints_starts[1] = (int)model.nonlinear_elements_operators.size();

    // Constraint 1 (linear): 2 x0 = 1.
    model.constraints_starts[1] = 0;
    model.elements_variables[0] = 0;
    model.elements_coefficients[0] = 2.0;
    model.constraints_lower_bounds[1] = 1.0;
    model.constraints_upper_bounds[1] = 1.0;

    return TestModel{"nonlinear_mixed_integer_infeasible", model, {}, false};
}

TestModel mathoptsolverscmaketest::black_box_continuous_model()
{
    // minimize (x0 - 3)^2 + (x1 + 2)^2
    // s.t.     -10 <= x0, x1 <= 10
    // Optimal: x0 = 3, x1 = -2, objective = 0.
    MathOptModel model(2, 0, 0);
    model.objective_direction = ObjectiveDirection::Minimize;

    model.variables_lower_bounds[0] = -10.0;
    model.variables_upper_bounds[0] = 10.0;
    model.variables_lower_bounds[1] = -10.0;
    model.variables_upper_bounds[1] = 10.0;

    model.objective_function = [](const std::vector<double>& x)
    {
        BlackBoxFunctionOutput output;
        output.objective_value = (x[0] - 3.0) * (x[0] - 3.0) + (x[1] + 2.0) * (x[1] + 2.0);
        output.gradient = {2.0 * (x[0] - 3.0), 2.0 * (x[1] + 2.0)};
        return output;
    };

    return TestModel{"black_box_continuous", model, {3.0, -2.0}};
}

TestModel mathoptsolverscmaketest::black_box_infeasible_model()
{
    // minimize (x0 - 3)^2 + (x1 + 2)^2
    // s.t.     x0 <= 1
    //          x0 >= 5
    // Infeasible regardless of the (otherwise unconstrained) black-box
    // objective.
    MathOptModel model(2, 2, 2);
    model.objective_direction = ObjectiveDirection::Minimize;

    model.variables_lower_bounds[0] = -10.0;
    model.variables_upper_bounds[0] = 10.0;
    model.variables_lower_bounds[1] = -10.0;
    model.variables_upper_bounds[1] = 10.0;

    model.objective_function = [](const std::vector<double>& x)
    {
        BlackBoxFunctionOutput output;
        output.objective_value = (x[0] - 3.0) * (x[0] - 3.0) + (x[1] + 2.0) * (x[1] + 2.0);
        output.gradient = {2.0 * (x[0] - 3.0), 2.0 * (x[1] + 2.0)};
        return output;
    };

    model.constraints_starts[0] = 0;
    model.elements_variables[0] = 0;
    model.elements_coefficients[0] = 1.0;
    model.constraints_upper_bounds[0] = 1.0;

    model.constraints_starts[1] = 1;
    model.elements_variables[1] = 0;
    model.elements_coefficients[1] = 1.0;
    model.constraints_lower_bounds[1] = 5.0;

    return TestModel{"black_box_infeasible", model, {}, false};
}

TestModel mathoptsolverscmaketest::black_box_mixed_integer_model()
{
    // minimize (x0 - 3.2)^2 + (x1 + 2)^2, x0 integer
    // s.t.     -10 <= x0, x1 <= 10
    // The integer optimum for x0 is 3 (|3 - 3.2| < |4 - 3.2|).
    // Optimal: x0 = 3, x1 = -2, objective = 0.04.
    MathOptModel model(2, 0, 0);
    model.objective_direction = ObjectiveDirection::Minimize;

    model.variables_lower_bounds[0] = -10.0;
    model.variables_upper_bounds[0] = 10.0;
    model.variables_types[0] = VariableType::Integer;
    model.variables_lower_bounds[1] = -10.0;
    model.variables_upper_bounds[1] = 10.0;

    model.objective_function = [](const std::vector<double>& x)
    {
        BlackBoxFunctionOutput output;
        output.objective_value = (x[0] - 3.2) * (x[0] - 3.2) + (x[1] + 2.0) * (x[1] + 2.0);
        output.gradient = {2.0 * (x[0] - 3.2), 2.0 * (x[1] + 2.0)};
        return output;
    };

    return TestModel{"black_box_mixed_integer", model, {3.0, -2.0}};
}

TestModel mathoptsolverscmaketest::black_box_mixed_integer_infeasible_model()
{
    // minimize (x0 - 3.2)^2 + (x1 + 2)^2, x0 integer
    // s.t.     2 x0 = 1
    // The continuous relaxation is feasible (x0 = 0.5), but no integer x0
    // satisfies 2 x0 = 1.
    MathOptModel model(2, 1, 1);
    model.objective_direction = ObjectiveDirection::Minimize;

    model.variables_lower_bounds[0] = -10.0;
    model.variables_upper_bounds[0] = 10.0;
    model.variables_types[0] = VariableType::Integer;
    model.variables_lower_bounds[1] = -10.0;
    model.variables_upper_bounds[1] = 10.0;

    model.objective_function = [](const std::vector<double>& x)
    {
        BlackBoxFunctionOutput output;
        output.objective_value = (x[0] - 3.2) * (x[0] - 3.2) + (x[1] + 2.0) * (x[1] + 2.0);
        output.gradient = {2.0 * (x[0] - 3.2), 2.0 * (x[1] + 2.0)};
        return output;
    };

    model.constraints_starts[0] = 0;
    model.elements_variables[0] = 0;
    model.elements_coefficients[0] = 2.0;
    model.constraints_lower_bounds[0] = 1.0;
    model.constraints_upper_bounds[0] = 1.0;

    return TestModel{"black_box_mixed_integer_infeasible", model, {}, false};
}

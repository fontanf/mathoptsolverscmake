#include "mathoptsolverscmake/mathopt_build.hpp"
#include "mathoptsolverscmake/mathopt_json.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <sstream>

using namespace mathoptsolverscmake;

namespace
{

const double INF = std::numeric_limits<double>::infinity();

std::string to_json(const MathOptModel& model)
{
    std::ostringstream out;
    write_json(model, out);
    return out.str();
}

/** max x + 2 y + 3 z s.t. x + y + z <= 10, x - y >= -3. */
MathOptModel small_milp()
{
    MathOptModel model;
    double upper_bounds[] = {10.0, 1.0, 5.0};
    signed char types[] = {0, 1, 2};
    std::vector<std::string> names = {"x", "y", "z"};
    VariablesBatch variables;
    variables.number = 3;
    variables.lower_bound = 0.0;
    variables.upper_bounds = upper_bounds;
    variables.types = types;
    variables.names = &names;
    EXPECT_EQ(add_variables(model, variables), 0);

    double objective[] = {1.0, 2.0, 3.0};
    ObjectiveBatch objective_batch;
    objective_batch.set_direction = true;
    objective_batch.direction = ObjectiveDirection::Maximize;
    objective_batch.set_coefficients = true;
    objective_batch.coefficients = objective;
    set_objective(model, objective_batch);

    double lower_bounds[] = {-INF, -3.0};
    double constraints_upper_bounds[] = {10.0, INF};
    int starts[] = {0, 3};
    int elements_variables[] = {0, 1, 2, 0, 1};
    double coefficients[] = {1, 1, 1, 1, -1};
    ConstraintsBatch constraints;
    constraints.number = 2;
    constraints.lower_bounds = lower_bounds;
    constraints.upper_bounds = constraints_upper_bounds;
    constraints.starts = starts;
    constraints.number_of_elements = 5;
    constraints.variables = elements_variables;
    constraints.coefficients = coefficients;
    EXPECT_EQ(add_constraints(model, constraints), 0);
    return model;
}

}

TEST(MathOptBuild, SmallMilp)
{
    MathOptModel model = small_milp();
    EXPECT_TRUE(model.check(1));
    EXPECT_EQ(model.objective_direction, ObjectiveDirection::Maximize);
    EXPECT_EQ(model.variables_types[2], VariableType::Integer);
    EXPECT_EQ(model.variables_names[1], "y");
    EXPECT_EQ(model.constraints_names, std::vector<std::string>());
    EXPECT_EQ(model.constraint_end(0), 3);
    EXPECT_EQ(model.evaluate_objective({1.0, 1.0, 2.0}), 9.0);
    EXPECT_EQ(model.evaluate_constraint({1.0, 1.0, 2.0}, 1), 0.0);
}

TEST(MathOptBuild, AppendedBatchesAndNames)
{
    MathOptModel model = small_milp();
    VariablesBatch variables;
    variables.number = 2;
    EXPECT_EQ(add_variables(model, variables), 3);
    EXPECT_EQ(model.variables_names, (std::vector<std::string>{"x", "y", "z", "x3", "x4"}));
    EXPECT_EQ(model.variables_lower_bounds[4], -INF);

    std::vector<std::string> names = {"last"};
    ConstraintsBatch constraints;
    constraints.number = 1;
    constraints.upper_bound = 1.0;
    constraints.names = &names;
    EXPECT_EQ(add_constraints(model, constraints), 2);
    EXPECT_EQ(model.constraints_names, (std::vector<std::string>{"c0", "c1", "last"}));
    EXPECT_EQ(model.constraints_starts, (std::vector<int>{0, 3, 5}));
    EXPECT_EQ(model.quadratic_elements_constraints_starts, (std::vector<int>{0, 0, 0}));
    EXPECT_TRUE(model.check(1));
}

TEST(MathOptBuild, NonlinearAndDefinedVariables)
{
    MathOptModel model;
    VariablesBatch variables;
    variables.number = 2;
    variables.lower_bound = 0.0;
    add_variables(model, variables);

    // A linear constraint, then exp(x0) + x1 ^ 2 <= 4.
    ConstraintsBatch linear;
    linear.number = 1;
    add_constraints(model, linear);
    int nonlinear_starts[] = {0};
    char operators[] = {'+', 'e', 'v', 'p', 'v', 'k'};
    int number_of_children[] = {2, 1, 0, 2, 0, 0};
    int tree_variables[] = {-1, -1, 0, -1, 1, -1};
    double values[] = {0, 0, 0, 0, 0, 2};
    ConstraintsBatch nonlinear;
    nonlinear.number = 1;
    nonlinear.upper_bound = 4.0;
    nonlinear.nonlinear_starts = nonlinear_starts;
    nonlinear.nonlinear.number_of_nodes = 6;
    nonlinear.nonlinear.operators = operators;
    nonlinear.nonlinear.number_of_children = number_of_children;
    nonlinear.nonlinear.variables = tree_variables;
    nonlinear.nonlinear.values = values;
    EXPECT_EQ(add_constraints(model, nonlinear), 1);
    EXPECT_EQ(model.nonlinear_elements_constraints_starts, (std::vector<int>{0, 0}));
    EXPECT_TRUE(model.check(1));
    EXPECT_NEAR(model.evaluate_constraint({1.0, 2.0}, 1), std::exp(1.0) + 4.0, 1e-12);

    // d0 = x0 * x1; the objective is d0.
    int defined_starts[] = {0};
    char defined_operators[] = {'*', 'v', 'v'};
    int defined_children[] = {2, 0, 0};
    int defined_variables[] = {-1, 0, 1};
    DefinedVariablesBatch defined;
    defined.number = 1;
    defined.starts = defined_starts;
    defined.trees.number_of_nodes = 3;
    defined.trees.operators = defined_operators;
    defined.trees.number_of_children = defined_children;
    defined.trees.variables = defined_variables;
    EXPECT_EQ(add_defined_variables(model, defined), 0);

    char objective_operators[] = {'d'};
    int objective_children[] = {0};
    int objective_variables[] = {0};
    ObjectiveBatch objective;
    objective.set_nonlinear = true;
    objective.nonlinear.number_of_nodes = 1;
    objective.nonlinear.operators = objective_operators;
    objective.nonlinear.number_of_children = objective_children;
    objective.nonlinear.variables = objective_variables;
    set_objective(model, objective);
    EXPECT_TRUE(model.check(1));
    EXPECT_EQ(model.evaluate_objective({2.0, 3.0}), 6.0);
}

TEST(MathOptBuild, InvalidBatchesLeaveTheModelUnchanged)
{
    MathOptModel model = small_milp();
    std::string before = to_json(model);

    int starts[] = {0};
    int bad_variables[] = {3};
    double coefficients[] = {1.0};
    ConstraintsBatch constraints;
    constraints.number = 1;
    constraints.starts = starts;
    constraints.number_of_elements = 1;
    constraints.variables = bad_variables;
    constraints.coefficients = coefficients;
    EXPECT_THROW(add_constraints(model, constraints), std::invalid_argument);

    // A tree with two roots.
    int nonlinear_starts[] = {0};
    char operators[] = {'v', 'v'};
    int number_of_children[] = {0, 0};
    int tree_variables[] = {0, 1};
    ConstraintsBatch nonlinear;
    nonlinear.number = 1;
    nonlinear.nonlinear_starts = nonlinear_starts;
    nonlinear.nonlinear.number_of_nodes = 2;
    nonlinear.nonlinear.operators = operators;
    nonlinear.nonlinear.number_of_children = number_of_children;
    nonlinear.nonlinear.variables = tree_variables;
    EXPECT_THROW(add_constraints(model, nonlinear), std::invalid_argument);

    signed char types[] = {3};
    VariablesBatch variables;
    variables.number = 1;
    variables.types = types;
    EXPECT_THROW(add_variables(model, variables), std::invalid_argument);

    // The objective: a valid direction with invalid quadratic terms changes nothing.
    int variables_1[] = {0};
    int variables_2[] = {7};
    ObjectiveBatch objective;
    objective.set_direction = true;
    objective.direction = ObjectiveDirection::Minimize;
    objective.set_quadratic = true;
    objective.number_of_quadratic_elements = 1;
    objective.quadratic_variables_1 = variables_1;
    objective.quadratic_variables_2 = variables_2;
    objective.quadratic_coefficients = coefficients;
    EXPECT_THROW(set_objective(model, objective), std::invalid_argument);

    EXPECT_EQ(to_json(model), before);
    EXPECT_EQ(model.objective_direction, ObjectiveDirection::Maximize);
}

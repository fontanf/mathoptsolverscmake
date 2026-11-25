#include "test_models.hpp"
#include "semi_test_models.hpp"

#include "mathoptsolverscmake/mathopt_cbc.hpp"
#include "mathoptsolverscmake/mathopt_mps.hpp"

#include "coin/CbcEventHandler.hpp"

#include <gtest/gtest.h>

#include <sstream>

#include <limits>

using namespace mathoptsolverscmake;
using namespace mathoptsolverscmaketest;

namespace
{

/**
 * Get the same model with the opposite objective direction and the opposite
 * objective, so that it has the same optimal solution.
 */
TestModel opposite_objective_direction(const TestModel& test_model)
{
    TestModel opposite_test_model = test_model;
    opposite_test_model.name += "_opposite";
    MathOptModel& model = opposite_test_model.model;
    model.objective_direction =
        (model.objective_direction == ObjectiveDirection::Minimize)?
        ObjectiveDirection::Maximize:
        ObjectiveDirection::Minimize;
    for (double& coefficient: model.objective_coefficients)
        coefficient = -coefficient;
    return opposite_test_model;
}

/** Get 1 for a minimization problem and -1 for a maximization problem. */
double objective_sense(const MathOptModel& model)
{
    return (model.objective_direction == ObjectiveDirection::Minimize)? 1.0: -1.0;
}

class CbcSolveTest: public testing::TestWithParam<TestModel> { };

TEST_P(CbcSolveTest, MatchesKnownResult)
{
    const TestModel& test_model = GetParam();

    mathoptsolverscmake::CbcModel cbc_model;
    reduce_printout(cbc_model);
    set_time_limit(cbc_model, 30.0);
    load(cbc_model, test_model.model);
    solve(cbc_model);

    if (!test_model.feasible) {
        EXPECT_TRUE(get_solution(cbc_model).empty());
        EXPECT_EQ(
                get_solution_value(cbc_model),
                objective_sense(test_model.model) * std::numeric_limits<double>::infinity());
        return;
    }

    double expected_objective_value =
        test_model.model.evaluate_objective(test_model.optimal_solution);
    EXPECT_NEAR(
            get_solution_value(cbc_model),
            expected_objective_value,
            1e-4);

    // The model is solved to optimality, so the bound is the optimal
    // objective value.
    EXPECT_NEAR(
            get_bound(cbc_model),
            expected_objective_value,
            1e-4);

    std::vector<double> solution = get_solution(cbc_model);
    ASSERT_EQ(solution.size(), test_model.optimal_solution.size());
    for (size_t variable_id = 0; variable_id < solution.size(); ++variable_id) {
        EXPECT_NEAR(
                solution[variable_id],
                test_model.optimal_solution[variable_id],
                1e-3)
            << "variable " << variable_id;
    }
}

INSTANTIATE_TEST_SUITE_P(
        MathOptCbc,
        CbcSolveTest,
        testing::Values(
                linear_continuous_model(),
                opposite_objective_direction(linear_continuous_model()),
                linear_mixed_integer_model(),
                opposite_objective_direction(linear_mixed_integer_model()),
                linear_infeasible_model(),
                opposite_objective_direction(linear_infeasible_model()),
                linear_mixed_integer_infeasible_model(),
                opposite_objective_direction(linear_mixed_integer_infeasible_model())),
        [](const testing::TestParamInfo<TestModel>& info)
        {
            return info.param.name;
        });

/**
 * Build a two-dimensional 0/1 knapsack model which Cbc doesn't solve at the
 * root node, so that the event handler is called during the search.
 */
MathOptModel two_dimensional_knapsack_model(ObjectiveDirection objective_direction)
{
    const int number_of_items = 40;
    MathOptModel model(number_of_items, 2, 2 * number_of_items);
    model.objective_direction = objective_direction;
    double sign = (objective_direction == ObjectiveDirection::Maximize)? 1.0: -1.0;
    unsigned seed = 1;
    auto random = [&seed]()
    {
        seed = seed * 1103515245 + 12345;
        return (double)(100 + (seed >> 16) % 900);
    };
    model.constraints_starts[0] = 0;
    model.constraints_starts[1] = number_of_items;
    for (int item_id = 0; item_id < number_of_items; ++item_id) {
        double weight_1 = random();
        double weight_2 = random();
        model.variables_lower_bounds[item_id] = 0.0;
        model.variables_upper_bounds[item_id] = 1.0;
        model.variables_types[item_id] = VariableType::Binary;
        model.objective_coefficients[item_id] = sign * (weight_1 + 37 + item_id % 7);
        model.elements_variables[item_id] = item_id;
        model.elements_coefficients[item_id] = weight_1;
        model.elements_variables[number_of_items + item_id] = item_id;
        model.elements_coefficients[number_of_items + item_id] = weight_2;
    }
    model.constraints_upper_bounds[0] = 8000;
    model.constraints_upper_bounds[1] = 8500;
    return model;
}

/** Event handler which records the solution values and bounds. */
class RecordEventHandler: public CbcEventHandler
{

public:

    RecordEventHandler(
            std::vector<double>& solution_values,
            std::vector<double>& bounds):
        solution_values_(solution_values),
        bounds_(bounds) { }

    virtual CbcAction event(CbcEvent)
    {
        const ::CbcModel& cbc_model = *model_;
        if ((cbc_model.specialOptions() & 2048) != 0)
            return noAction;
        solution_values_.push_back(get_solution_value(cbc_model));
        bounds_.push_back(get_bound(cbc_model));
        return noAction;
    }

    virtual CbcEventHandler* clone() const { return new RecordEventHandler(*this); }

private:

    std::vector<double>& solution_values_;
    std::vector<double>& bounds_;

};

class CbcEventHandlerTest: public testing::TestWithParam<ObjectiveDirection> { };

TEST_P(CbcEventHandlerTest, SolutionValuesAndBoundsHaveTheRightSign)
{
    MathOptModel model = two_dimensional_knapsack_model(GetParam());
    double sense = objective_sense(model);

    std::vector<double> solution_values;
    std::vector<double> bounds;
    mathoptsolverscmake::CbcModel cbc_model;
    reduce_printout(cbc_model);
    load(cbc_model, model);
    RecordEventHandler event_handler(solution_values, bounds);
    cbc_model.model.passInEventHandler(&event_handler);
    solve(cbc_model);

    double optimal_value = get_solution_value(cbc_model);
    EXPECT_NEAR(get_bound(cbc_model), optimal_value, 1e-4);
    EXPECT_NEAR(model.evaluate_objective(get_solution(cbc_model)), optimal_value, 1e-4);

    // The search must not end at the root node, otherwise this test doesn't
    // check anything.
    ASSERT_GT(get_number_of_nodes(cbc_model), 1);
    ASSERT_FALSE(solution_values.empty());

    // During the search, solution values are not better than the optimal
    // value, and bounds are not worse.
    for (double solution_value: solution_values)
        EXPECT_GE(sense * solution_value, sense * optimal_value - 1e-4);
    for (double bound: bounds)
        EXPECT_LE(sense * bound, sense * optimal_value + 1e-4);
}

INSTANTIATE_TEST_SUITE_P(
        MathOptCbc,
        CbcEventHandlerTest,
        testing::Values(
                ObjectiveDirection::Minimize,
                ObjectiveDirection::Maximize),
        [](const testing::TestParamInfo<ObjectiveDirection>& info)
        {
            return (info.param == ObjectiveDirection::Minimize)? "Minimize": "Maximize";
        });

}

TEST(MathOptCbc, SemiVariablesThrow)
{
    std::istringstream in(mathoptsolverscmaketest::semi_mps(false));
    mathoptsolverscmake::MathOptModel model = mathoptsolverscmake::read_mps(in);
    mathoptsolverscmake::CbcModel solver_model;
    EXPECT_THROW(load(solver_model, model), std::invalid_argument);
}

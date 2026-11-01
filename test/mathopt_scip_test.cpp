#include "test_models.hpp"

#include "mathoptsolverscmake/mathopt_scip.hpp"

#include <gtest/gtest.h>

using namespace mathoptsolverscmake;
using namespace mathoptsolverscmaketest;

namespace
{

class ScipSolveTest: public testing::TestWithParam<TestModel> { };

TEST_P(ScipSolveTest, MatchesKnownResult)
{
    const TestModel& test_model = GetParam();

    ScipModel scip_model;
    reduce_printout(scip_model);
    set_time_limit(scip_model, 30.0);
    load(scip_model, test_model.model);
    solve(scip_model);

    if (!test_model.feasible) {
        EXPECT_TRUE(is_infeasible(scip_model));
        return;
    }

    // The known solution must itself be feasible for the model that
    // produced it.
    EXPECT_TRUE(test_model.model.check_solution(test_model.optimal_solution));

    double expected_objective_value =
        test_model.model.evaluate_objective(test_model.optimal_solution);
    EXPECT_NEAR(
            get_solution_value(scip_model),
            expected_objective_value,
            1e-4);

    std::vector<double> solution = get_solution(scip_model);
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
        MathOptScip,
        ScipSolveTest,
        testing::Values(
                linear_continuous_model(),
                linear_mixed_integer_model(),
                linear_infeasible_model(),
                linear_mixed_integer_infeasible_model(),
                quadratic_continuous_model(),
                quadratic_mixed_integer_model(),
                quadratic_infeasible_model(),
                quadratic_mixed_integer_infeasible_model(),
                nonlinear_continuous_model(),
                nonlinear_mixed_integer_model(),
                nonlinear_infeasible_model(),
                nonlinear_mixed_integer_infeasible_model(),
                nonlinear_operators_model(),
                nonconvex_continuous_model(),
                nonconvex_mixed_integer_model()),
        [](const testing::TestParamInfo<TestModel>& info)
        {
            return info.param.name;
        });

}

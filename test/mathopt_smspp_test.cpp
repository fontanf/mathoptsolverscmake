#include "test_models.hpp"

#include "mathoptsolverscmake/mathopt_smspp.hpp"

#include <gtest/gtest.h>

using namespace mathoptsolverscmake;
using namespace mathoptsolverscmaketest;

namespace
{

class SmsppSolveTest: public testing::TestWithParam<TestModel> { };

TEST_P(SmsppSolveTest, MatchesKnownResult)
{
    const TestModel& test_model = GetParam();

    SmsppOutput output = solve_smspp(test_model.model);

    // The known solution must itself be feasible for the model that
    // produced it.
    EXPECT_TRUE(test_model.model.check_solution(test_model.optimal_solution));

    double expected_objective_value =
        test_model.model.evaluate_objective(test_model.optimal_solution);
    EXPECT_NEAR(
            output.objective_value,
            expected_objective_value,
            1e-4);

    ASSERT_EQ(output.solution.size(), test_model.optimal_solution.size());
    for (size_t variable_id = 0; variable_id < output.solution.size(); ++variable_id) {
        EXPECT_NEAR(
                output.solution[variable_id],
                test_model.optimal_solution[variable_id],
                1e-3)
            << "variable " << variable_id;
    }
}

// SMS++'s BundleSolver only supports box-constrained continuous convex
// models whose variable lower bounds are exactly 0 or -infinity (the
// Lagrangian-multiplier convention its box-constraint handling is built
// around, see mathopt_smspp.hpp), so only black_box_convex_model() applies
// here -- unlike Knitro, none of the other test_models.hpp fixtures (which
// use arbitrary finite bounds, non-convex cases, or integrality) fit.

INSTANTIATE_TEST_SUITE_P(
        MathOptSmspp,
        SmsppSolveTest,
        testing::Values(
                black_box_convex_model()),
        [](const testing::TestParamInfo<TestModel>& info)
        {
            return info.param.name;
        });

}

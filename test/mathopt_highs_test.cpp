#include "test_models.hpp"
#include "semi_test_models.hpp"

#include "mathoptsolverscmake/mathopt_highs.hpp"
#include "mathoptsolverscmake/mathopt_lp.hpp"
#include "mathoptsolverscmake/mathopt_mps.hpp"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <thread>

using namespace mathoptsolverscmake;
using namespace mathoptsolverscmaketest;

namespace
{

double get_double_option(
        const Highs& highs_model,
        const std::string& option)
{
    double value = 0;
    EXPECT_EQ(highs_model.getOptionValue(option, value), HighsStatus::kOk);
    return value;
}

}

TEST(MathOptHighs, DefaultTolerancesKeepHighsDefaults)
{
    // The tolerances of the model are 0 by default, which HiGHS rejects:
    // loading the model must keep the default tolerances of HiGHS, without
    // trying to set them (HiGHS logs an error for each rejected value).
    std::filesystem::path log_path =
        std::filesystem::temp_directory_path() / "mathopt_highs_test_default_tolerances.log";
    std::filesystem::remove(log_path);

    Highs default_highs_model;
    double default_primal_feasibility_tolerance = get_double_option(
            default_highs_model, "primal_feasibility_tolerance");
    double default_mip_feasibility_tolerance = get_double_option(
            default_highs_model, "mip_feasibility_tolerance");

    TestModel test_model = linear_mixed_integer_model();
    ASSERT_EQ(test_model.model.feasibility_tolerance, 0.0);
    ASSERT_EQ(test_model.model.integrality_tolerance, 0.0);
    Highs highs_model;
    reduce_printout(highs_model);
    set_log_file(highs_model, log_path.string());
    load(highs_model, test_model.model);
    highs_model.setOptionValue("log_file", "");

    std::ifstream log_file(log_path);
    std::stringstream log;
    log << log_file.rdbuf();
    EXPECT_EQ(log.str().find("ERROR"), std::string::npos) << log.str();
    std::filesystem::remove(log_path);

    EXPECT_EQ(
            get_double_option(highs_model, "primal_feasibility_tolerance"),
            default_primal_feasibility_tolerance);
    EXPECT_EQ(
            get_double_option(highs_model, "mip_feasibility_tolerance"),
            default_mip_feasibility_tolerance);
}

TEST(MathOptHighs, TolerancesArePassed)
{
    TestModel test_model = linear_mixed_integer_model();
    test_model.model.feasibility_tolerance = 1e-5;
    test_model.model.integrality_tolerance = 1e-4;
    Highs highs_model;
    reduce_printout(highs_model);
    load(highs_model, test_model.model);

    EXPECT_EQ(get_double_option(highs_model, "primal_feasibility_tolerance"), 1e-5);
    EXPECT_EQ(get_double_option(highs_model, "mip_feasibility_tolerance"), 1e-4);
}

class HighsSolveTest: public testing::TestWithParam<TestModel> { };

TEST_P(HighsSolveTest, MatchesKnownResult)
{
    const TestModel& test_model = GetParam();

    Highs highs_model;
    reduce_printout(highs_model);
    set_time_limit(highs_model, 30.0);
    load(highs_model, test_model.model);
    solve(highs_model);

    std::vector<double> solution = get_solution(highs_model);
    ASSERT_EQ(solution.size(), test_model.optimal_solution.size());
    EXPECT_NEAR(
            test_model.model.evaluate_objective(solution),
            test_model.model.evaluate_objective(test_model.optimal_solution),
            1e-4);
    for (size_t variable_id = 0; variable_id < solution.size(); ++variable_id) {
        EXPECT_NEAR(
                solution[variable_id],
                test_model.optimal_solution[variable_id],
                1e-3)
            << "variable " << variable_id;
    }
}

INSTANTIATE_TEST_SUITE_P(
        MathOptHighs,
        HighsSolveTest,
        testing::Values(
                linear_continuous_model(),
                linear_mixed_integer_model()),
        [](const testing::TestParamInfo<TestModel>& info)
        {
            return info.param.name;
        });

namespace
{

/**
 * Run 'function' in a new thread, so that it starts without any HiGHS
 * scheduler, HiGHS keeping one scheduler per thread.
 */
template<typename Function>
void run_in_new_thread(Function function)
{
    std::thread thread(function);
    thread.join();
}

}

TEST(MathOptHighs, SolveLeavesNoMultiThreadedScheduler)
{
    // After a solve with several threads, a solve with a single thread in the
    // same thread must not fail (HiGHS returns an error when the 'threads'
    // option differs from the size of the scheduler of the thread).
    run_in_new_thread([]()
    {
        TestModel test_model = linear_mixed_integer_model();
        Highs highs_model;
        reduce_printout(highs_model);
        highs_model.setOptionValue("threads", 2);
        load(highs_model, test_model.model);
        solve(highs_model);
        EXPECT_EQ(highs_model.getModelStatus(), HighsModelStatus::kOptimal);
        // The solution remains available after the scheduler reset.
        std::vector<double> solution = get_solution(highs_model);
        ASSERT_EQ(solution.size(), test_model.optimal_solution.size());
        EXPECT_NEAR(
                test_model.model.evaluate_objective(solution),
                test_model.model.evaluate_objective(test_model.optimal_solution),
                1e-4);

        Highs single_thread_highs_model;
        reduce_printout(single_thread_highs_model);
        single_thread_highs_model.setOptionValue("threads", 1);
        load(single_thread_highs_model, test_model.model);
        EXPECT_EQ(single_thread_highs_model.run(), HighsStatus::kOk);
        EXPECT_EQ(single_thread_highs_model.getModelStatus(), HighsModelStatus::kOptimal);
    });
}

TEST(MathOptHighs, SolveDoesNotReuseSingleThreadScheduler)
{
    // After a solve with a single thread outside 'solve', which leaves a
    // single-thread scheduler in the thread, a solve with several threads
    // must not fail.
    run_in_new_thread([]()
    {
        TestModel test_model = linear_mixed_integer_model();
        Highs single_thread_highs_model;
        reduce_printout(single_thread_highs_model);
        single_thread_highs_model.setOptionValue("threads", 1);
        load(single_thread_highs_model, test_model.model);
        ASSERT_EQ(single_thread_highs_model.run(), HighsStatus::kOk);

        Highs highs_model;
        reduce_printout(highs_model);
        highs_model.setOptionValue("threads", 2);
        load(highs_model, test_model.model);
        solve(highs_model);
        EXPECT_EQ(highs_model.getModelStatus(), HighsModelStatus::kOptimal);
    });
}

TEST(MathOptHighs, SemiVariables)
{
    // HiGHS's semi-continuous and semi-integer models, read from MPS and LP.
    for (bool semi_integer: {false, true}) {
        for (bool lp: {false, true}) {
            std::istringstream in(lp? semi_lp(semi_integer): semi_mps(semi_integer));
            MathOptModel model = lp? read_lp(in): read_mps(in);
            Highs highs_model;
            reduce_printout(highs_model);
            load(highs_model, model);
            solve(highs_model);
            std::vector<double> solution = get_solution(highs_model);
            ASSERT_EQ(solution.size(), 4u);
            double expected = semi_integer? SEMI_INTEGER_OPTIMUM: SEMI_CONTINUOUS_OPTIMUM;
            EXPECT_NEAR(model.evaluate_objective(solution), expected, 1e-4)
                << (semi_integer? "semi-integer": "semi-continuous") << (lp? " lp": " mps");
            model.feasibility_tolerance = 1e-6;
            model.integrality_tolerance = 1e-6;
            EXPECT_TRUE(model.check_solution(solution));
        }
    }
}

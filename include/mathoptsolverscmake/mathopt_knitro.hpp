#pragma once

extern "C"
{
#include "knitro.h"
}

#include "mathoptsolverscmake/mathopt.hpp"

#include <functional>
#include <memory>
#include <vector>

namespace mathoptsolverscmake
{

struct KnitroContext
{
    KN_context* kc = nullptr;

    struct EvalCallbackStruct
    {
        std::function<int(KN_eval_request_ptr const, KN_eval_result_ptr const)> eval;
        std::function<int(KN_eval_request_ptr const, KN_eval_result_ptr const)> grad;
    };

    std::vector<double> x_;
    std::vector<std::unique_ptr<EvalCallbackStruct>> eval_callbacks_;

    KnitroContext();
    ~KnitroContext();
    KnitroContext(const KnitroContext&) = delete;
    KnitroContext& operator=(const KnitroContext&) = delete;
};

void load(
        KnitroContext& knitro,
        const MathOptModel& model);

void reduce_printout(
        KnitroContext& knitro);

void set_time_limit(
        KnitroContext& knitro,
        double time_limit);

void solve(
        KnitroContext& knitro);

double get_solution_value(
        const KnitroContext& knitro);

std::vector<double> get_solution(
        const KnitroContext& knitro);

int get_number_of_nodes(
        const KnitroContext& knitro);

/**
 * Return true if Knitro considers the current solution infeasible. Unlike
 * an LP/MIP solver, an NLP solver cannot always certify infeasibility
 * algebraically: on top of the confirmed-infeasible status band, it may
 * instead exhaust its iteration/time/node limit while stuck at an
 * infeasible point, or report the ambiguous "unbounded or infeasible"
 * status; all of these are treated as infeasible here.
 */
bool is_infeasible(
        const KnitroContext& knitro);

}

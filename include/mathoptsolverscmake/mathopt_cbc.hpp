#pragma once

#include "mathoptsolverscmake/mathopt.hpp"

#include "coin/OsiCbcSolverInterface.hpp"

#include <vector>

namespace mathoptsolverscmake
{

void load(
        CbcModel& cbc_model,
        const MathOptModel& model);

void reduce_printout(
        CbcModel& cbc_model);

void set_time_limit(
        CbcModel& cbc_model,
        double time_limit);

void solve(
        CbcModel& cbc_model);

double get_solution_value(
        const CbcModel& cbc_model);

std::vector<double> get_solution(
        const CbcModel& cbc_model);

double get_bound(
        const CbcModel& cbc_model);

int get_number_of_nodes(
        const CbcModel& cbc_model);

}

#pragma once

#include "mathoptsolverscmake/mathopt.hpp"

#include "xprs.h"

#include <string>
#include <vector>

namespace mathoptsolverscmake
{

void load(
        XPRSprob& xpress_model,
        const MathOptModel& model);

void set_time_limit(
        XPRSprob& xpress_model,
        double time_limit);

void set_log_file(
        XPRSprob& xpress_model,
        const std::string& log_file);

void write_mps(
        const XPRSprob& xpress_model,
        const std::string& mps_file);

void solve(
        XPRSprob& xpress_model);

double get_solution_value(
        const XPRSprob& xpress_model);

std::vector<double> get_solution(
        const XPRSprob& xpress_model);

double get_bound(
        const XPRSprob& xpress_model);

int get_number_of_nodes(
        const XPRSprob& xpress_model);

}

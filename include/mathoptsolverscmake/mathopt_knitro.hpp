#pragma once

#include "mathoptsolverscmake/mathopt.hpp"

#include "knitrocpp/knitro.hpp"

#include <vector>

namespace mathoptsolverscmake
{

void solve(
        const MathOptModel& model,
        knitrocpp::Context& knitro_context);

double get_solution_value(
        const knitrocpp::Context& knitro_context);

std::vector<double> get_solution(
        const knitrocpp::Context& knitro_context);

}

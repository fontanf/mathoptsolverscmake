#pragma once

#include "mathoptsolverscmake/mathopt.hpp"

#include "CBSolver.hxx"

#include <vector>

namespace mathoptsolverscmake
{

void solve(
        const MathOptModel& model,
        ConicBundle::CBSolver& solver);

double get_solution_value(
        const MathOptModel& model,
        const ConicBundle::CBSolver& solver);

std::vector<double> get_solution(
        const MathOptModel& model,
        const ConicBundle::CBSolver& solver);

}

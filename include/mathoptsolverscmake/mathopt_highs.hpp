#pragma once

#include "mathoptsolverscmake/mathopt.hpp"

#include "Highs.h"

#include <string>
#include <vector>

namespace mathoptsolverscmake
{

void load(
        Highs& highs_model,
        const MathOptModel& model);

void reduce_printout(
        Highs& highs_model);

void set_time_limit(
        Highs& highs_model,
        double time_limit);

void set_node_limit(
        Highs& highs_model,
        int node_limit);

void set_log_file(
        Highs& highs_model,
        const std::string& log_file);

void write_mps(
        Highs& highs_model,
        const std::string& mps_file);

void solve(
        Highs& highs_model);

/** Return an empty vector if HiGHS has not found any feasible solution. */
std::vector<double> get_solution(
        const Highs& highs_model);

/**
 * Return the bound on the optimal objective value: HiGHS's MIP dual bound for
 * a MIP; for an LP, the optimal value if HiGHS solved it to optimality, and
 * an infinite bound otherwise.
 */
double get_bound(
        const Highs& highs_model);

int get_number_of_nodes(
        const Highs& highs_model);

MathOptModel to_mathopt(
        Highs& highs_model);

}

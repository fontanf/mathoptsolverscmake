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

void set_solution(
        Highs& highs_model,
        const std::vector<double>& solution);

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

std::vector<double> get_solution(
        const Highs& highs_model);

double get_bound(
        const Highs& highs_model);

}

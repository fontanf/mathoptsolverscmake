#pragma once

#include "mathoptsolverscmake/mathopt.hpp"

#include "ilcplex/cplex.h"

#include <string>
#include <vector>

namespace mathoptsolverscmake
{

/**
 * CPLEX environment and problem built from a MathOptModel.
 *
 * CPLEX solves MIPs with linear and quadratic terms. Nonlinear expression
 * trees and black-box functions are not supported.
 *
 * A non-convex quadratic objective is solved to global optimality (load()
 * sets the parameter CPXPARAM_OptimalityTarget to
 * CPX_OPTIMALITYTARGET_OPTIMALGLOBAL when the objective has quadratic
 * terms). Quadratic constraints must be convex.
 *
 * Parameters can be set with CPXset*param(env, ...).
 */
struct CplexModel
{
    CPXENVptr env = nullptr;

    CPXLPptr lp = nullptr;

    CplexModel();
    ~CplexModel();
    CplexModel(const CplexModel&) = delete;
    CplexModel& operator=(const CplexModel&) = delete;
};

void load(
        CplexModel& cplex_model,
        const MathOptModel& model);

void reduce_printout(
        CplexModel& cplex_model);

void set_time_limit(
        CplexModel& cplex_model,
        double time_limit);

void set_node_limit(
        CplexModel& cplex_model,
        int node_limit);

/**
 * Write the log of CPLEX to a file.
 *
 * It can be combined with reduce_printout, which only disables the output to
 * the standard output.
 */
void set_log_file(
        CplexModel& cplex_model,
        const std::string& log_file);

void solve(
        CplexModel& cplex_model);

double get_solution_value(
        const CplexModel& cplex_model);

/** Return an empty vector if CPLEX has not found any feasible solution. */
std::vector<double> get_solution(
        const CplexModel& cplex_model);

double get_bound(
        const CplexModel& cplex_model);

int get_number_of_nodes(
        const CplexModel& cplex_model);

/**
 * Return true if CPLEX proved that the model is infeasible. The ambiguous
 * "infeasible or unbounded" status is also treated as infeasible.
 */
bool is_infeasible(
        const CplexModel& cplex_model);

}

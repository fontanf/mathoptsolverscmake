#pragma once

#include "mathoptsolverscmake/mathopt.hpp"

#include "gurobi_c.h"

#include <string>
#include <vector>

namespace mathoptsolverscmake
{

/**
 * Gurobi environment and model built from a MathOptModel.
 *
 * Gurobi solves MINLPs with linear, quadratic and nonlinear expression tree
 * terms (Gurobi >= 12) to global optimality, up to its gap and feasibility
 * tolerances. Black-box functions are not supported.
 *
 * Nonlinear expression trees are added as general constraints y = f(x):
 * - the nonlinear part f(x) of a constraint is replaced by an auxiliary
 *   variable y
 * - the nonlinear part f(x) of the objective is replaced by an auxiliary
 *   variable z: minimize c x + x'Qx + z s.t. z = f(x)
 *
 * The environment is created by the constructor but only started by load(),
 * so that parameters set before load() (for example by reduce_printout())
 * also apply to the license messages. Parameters can be set with
 * GRBset*param(env, ...) before load() and GRBset*param(GRBgetenv(model),
 * ...) after.
 */
struct GurobiModel
{
    GRBenv* env = nullptr;

    /** Model, created by load(). */
    GRBmodel* model = nullptr;

    /** Number of variables of the MathOptModel. */
    int number_of_variables = 0;

    GurobiModel();
    ~GurobiModel();
    GurobiModel(const GurobiModel&) = delete;
    GurobiModel& operator=(const GurobiModel&) = delete;
};

void load(
        GurobiModel& gurobi_model,
        const MathOptModel& model);

void reduce_printout(
        GurobiModel& gurobi_model);

void set_time_limit(
        GurobiModel& gurobi_model,
        double time_limit);

void set_node_limit(
        GurobiModel& gurobi_model,
        int node_limit);

/**
 * Write the log of Gurobi to a file.
 *
 * It can be combined with reduce_printout, which only disables the output to
 * the standard output.
 */
void set_log_file(
        GurobiModel& gurobi_model,
        const std::string& log_file);

void solve(
        GurobiModel& gurobi_model);

double get_solution_value(
        const GurobiModel& gurobi_model);

/** Return an empty vector if Gurobi has not found any feasible solution. */
std::vector<double> get_solution(
        const GurobiModel& gurobi_model);

double get_bound(
        const GurobiModel& gurobi_model);

int get_number_of_nodes(
        const GurobiModel& gurobi_model);

/**
 * Return true if Gurobi proved that the model is infeasible. The ambiguous
 * "infeasible or unbounded" status is also treated as infeasible.
 */
bool is_infeasible(
        const GurobiModel& gurobi_model);

}

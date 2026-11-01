#pragma once

#include "mathoptsolverscmake/mathopt.hpp"

#include "scip/scip.h"

#include <string>
#include <vector>

namespace mathoptsolverscmake
{

/**
 * SCIP problem built from a MathOptModel.
 *
 * SCIP solves MINLPs (with linear, quadratic and nonlinear expression tree
 * terms; black-box functions are not supported) to global optimality, up to
 * its gap and feasibility tolerances ("limits/gap", "numerics/feastol"...),
 * which can be set through the "scip" field.
 *
 * SCIP only supports a linear objective. When the objective has quadratic or
 * nonlinear terms, load() moves them into a constraint with an auxiliary
 * variable (objective_variable): minimize c x + z s.t. f(x) <= z.
 */
struct ScipModel
{
    SCIP* scip = nullptr;

    /**
     * SCIP variables, in the order of the variables of the MathOptModel.
     *
     * They are owned by the SCIP problem.
     */
    std::vector<SCIP_VAR*> variables;

    /** Auxiliary objective variable (see above), or nullptr. */
    SCIP_VAR* objective_variable = nullptr;

    ScipModel();
    ~ScipModel();
    ScipModel(const ScipModel&) = delete;
    ScipModel& operator=(const ScipModel&) = delete;
};

void load(
        ScipModel& scip_model,
        const MathOptModel& model);

void reduce_printout(
        ScipModel& scip_model);

void set_time_limit(
        ScipModel& scip_model,
        double time_limit);

void set_node_limit(
        ScipModel& scip_model,
        int node_limit);

/**
 * Write the log of SCIP to a file.
 *
 * Unlike for Cbc, it can be combined with reduce_printout, which only
 * disables the output to the standard output.
 */
void set_log_file(
        ScipModel& scip_model,
        const std::string& log_file);

void solve(
        ScipModel& scip_model);

double get_solution_value(
        const ScipModel& scip_model);

/** Return an empty vector if SCIP has not found any feasible solution. */
std::vector<double> get_solution(
        const ScipModel& scip_model);

double get_bound(
        const ScipModel& scip_model);

int get_number_of_nodes(
        const ScipModel& scip_model);

/**
 * Return true if SCIP proved that the model is infeasible. The ambiguous
 * "infeasible or unbounded" status is also treated as infeasible.
 */
bool is_infeasible(
        const ScipModel& scip_model);

}

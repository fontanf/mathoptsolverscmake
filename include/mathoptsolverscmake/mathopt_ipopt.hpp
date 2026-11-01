#pragma once

#include "mathoptsolverscmake/mathopt.hpp"

#include "IpIpoptApplication.hpp"
#include "IpTNLP.hpp"

#include <limits>
#include <string>
#include <vector>

namespace mathoptsolverscmake
{

/**
 * Ipopt application and problem built from a MathOptModel.
 *
 * Ipopt is a local solver for continuous nonlinear programs: load() rejects
 * models with integer variables, and on non-convex models the returned
 * solution is only locally optimal.
 *
 * First and second derivatives of linear, quadratic and nonlinear
 * expression tree terms are computed exactly. Black-box functions only
 * provide first derivatives, so when the model has some, the Hessian is
 * approximated (Ipopt option "hessian_approximation" set to
 * "limited-memory").
 *
 * Ipopt options can be set through "application->Options()".
 */
struct IpoptModel
{
    Ipopt::SmartPtr<Ipopt::IpoptApplication> application;

    /** Problem, set by load(). */
    Ipopt::SmartPtr<Ipopt::TNLP> problem;

    /*
     * Results, set by solve().
     */

    Ipopt::ApplicationReturnStatus status = Ipopt::Internal_Error;
    std::vector<double> solution;
    double objective_value = std::numeric_limits<double>::quiet_NaN();

    IpoptModel();
    IpoptModel(const IpoptModel&) = delete;
    IpoptModel& operator=(const IpoptModel&) = delete;
};

void load(
        IpoptModel& ipopt_model,
        const MathOptModel& model);

void reduce_printout(
        IpoptModel& ipopt_model);

void set_time_limit(
        IpoptModel& ipopt_model,
        double time_limit);

void set_log_file(
        IpoptModel& ipopt_model,
        const std::string& log_file);

void solve(
        IpoptModel& ipopt_model);

double get_solution_value(
        const IpoptModel& ipopt_model);

std::vector<double> get_solution(
        const IpoptModel& ipopt_model);

/**
 * Return true if Ipopt converged to a point of local infeasibility. Ipopt is
 * a local solver, so this doesn't prove that the model is infeasible.
 */
bool is_infeasible(
        const IpoptModel& ipopt_model);

}

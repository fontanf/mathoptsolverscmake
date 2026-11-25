#include "mathoptsolverscmake/mathopt_solve.hpp"

#ifdef CBC_FOUND
#include "mathoptsolverscmake/mathopt_cbc.hpp"
#endif
#ifdef HIGHS_FOUND
#include "mathoptsolverscmake/mathopt_highs.hpp"
#endif
#ifdef XPRESS_FOUND
#include "mathoptsolverscmake/mathopt_xpress.hpp"
#endif
#ifdef CPLEX_FOUND
#include "mathoptsolverscmake/mathopt_cplex.hpp"
#endif
#ifdef GUROBI_FOUND
#include "mathoptsolverscmake/mathopt_gurobi.hpp"
#endif
#ifdef KNITRO_FOUND
#include "mathoptsolverscmake/mathopt_knitro.hpp"
#endif
#ifdef IPOPT_FOUND
#include "mathoptsolverscmake/mathopt_ipopt.hpp"
#endif
#ifdef SCIP_FOUND
#include "mathoptsolverscmake/mathopt_scip.hpp"
#endif

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace mathoptsolverscmake;

namespace
{

const SolverName SOLVERS[] = {
    SolverName::Cbc,
    SolverName::Highs,
    SolverName::Xpress,
    SolverName::Cplex,
    SolverName::Gurobi,
    SolverName::Knitro,
    SolverName::Ipopt,
    SolverName::Scip,
};

std::string solver_name(SolverName solver)
{
    std::ostringstream name;
    name << solver;
    return name.str();
}

#ifdef XPRESS_FOUND
/** Frees an Xpress problem, also when solving throws. */
struct XpressProblem
{
    XPRSprob prob = nullptr;

    XpressProblem() { XPRScreateprob(&prob); }

    ~XpressProblem() { XPRSdestroyprob(prob); }

    XpressProblem(const XpressProblem&) = delete;
    XpressProblem& operator=(const XpressProblem&) = delete;
};
#endif

}

bool mathoptsolverscmake::is_available(SolverName solver)
{
    switch (solver) {
#ifdef CBC_FOUND
    case SolverName::Cbc: return true;
#endif
#ifdef HIGHS_FOUND
    case SolverName::Highs: return true;
#endif
#ifdef XPRESS_FOUND
    case SolverName::Xpress: return true;
#endif
#ifdef CPLEX_FOUND
    case SolverName::Cplex: return true;
#endif
#ifdef GUROBI_FOUND
    case SolverName::Gurobi: return true;
#endif
#ifdef KNITRO_FOUND
    case SolverName::Knitro: return true;
#endif
#ifdef IPOPT_FOUND
    case SolverName::Ipopt: return true;
#endif
#ifdef SCIP_FOUND
    case SolverName::Scip: return true;
#endif
    default: return false;
    }
}

std::vector<SolverName> mathoptsolverscmake::available_solvers()
{
    std::vector<SolverName> solvers;
    for (SolverName solver: SOLVERS)
        if (is_available(solver))
            solvers.push_back(solver);
    return solvers;
}

SolveOutput mathoptsolverscmake::solve(
        const MathOptModel& model,
        const SolveParameters& parameters)
{
    SolverName solver = parameters.solver;
    bool has_time_limit = !std::isinf(parameters.time_limit);
    double time_limit = parameters.time_limit;
    if (std::isnan(time_limit) || time_limit < 0) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "the time limit must be non-negative.");
    }
    bool solver_known = false;
    for (SolverName s: SOLVERS)
        if (s == solver)
            solver_known = true;
    if (!solver_known) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "solver \"" + solver_name(solver) + "\" is not supported "
                "(only Cbc, HiGHS, XPRESS, CPLEX, Gurobi, Knitro, Ipopt and SCIP are).");
    }
    if (!is_available(solver)) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "solver \"" + solver_name(solver) + "\" is not available in this build.");
    }

    SolveOutput output;
    switch (solver) {
#ifdef CBC_FOUND
    case SolverName::Cbc: {
        CbcModel cbc_model;
        if (!parameters.verbose)
            reduce_printout(cbc_model);
        if (has_time_limit)
            set_time_limit(cbc_model, time_limit);
        load(cbc_model, model);
        solve(cbc_model);
        output.solution = get_solution(cbc_model);
        output.found_solution = !output.solution.empty();
        output.bound = get_bound(cbc_model);
        output.has_bound = true;
        break;
    }
#endif
#ifdef HIGHS_FOUND
    case SolverName::Highs: {
        Highs highs;
        if (!parameters.verbose)
            reduce_printout(highs);
        if (has_time_limit)
            set_time_limit(highs, time_limit);
        load(highs, model);
        solve(highs);
        output.solution = get_solution(highs);
        output.found_solution = !output.solution.empty();
        output.bound = get_bound(highs);
        output.has_bound = true;
        break;
    }
#endif
#ifdef XPRESS_FOUND
    case SolverName::Xpress: {
        XpressProblem xpress_model;
        if (!parameters.verbose)
            reduce_printout(xpress_model.prob);
        if (has_time_limit)
            set_time_limit(xpress_model.prob, time_limit);
        load(xpress_model.prob, model);
        solve(xpress_model.prob);
        // Unlike the other solvers, XPRESS's get_solution() throws rather
        // than returning an empty solution when no solution is available.
        try {
            output.solution = get_solution(xpress_model.prob);
            output.found_solution = true;
            output.bound = get_bound(xpress_model.prob);
            output.has_bound = true;
        } catch (const std::exception&) {
            output.solution.clear();
            output.found_solution = false;
        }
        break;
    }
#endif
#ifdef CPLEX_FOUND
    case SolverName::Cplex: {
        CplexModel cplex_model;
        if (!parameters.verbose)
            reduce_printout(cplex_model);
        if (has_time_limit)
            set_time_limit(cplex_model, time_limit);
        load(cplex_model, model);
        solve(cplex_model);
        output.solution = get_solution(cplex_model);
        output.found_solution = !output.solution.empty();
        output.bound = get_bound(cplex_model);
        output.has_bound = true;
        break;
    }
#endif
#ifdef GUROBI_FOUND
    case SolverName::Gurobi: {
        GurobiModel gurobi_model;
        if (!parameters.verbose)
            reduce_printout(gurobi_model);
        if (has_time_limit)
            set_time_limit(gurobi_model, time_limit);
        load(gurobi_model, model);
        solve(gurobi_model);
        output.solution = get_solution(gurobi_model);
        output.found_solution = !output.solution.empty();
        output.bound = get_bound(gurobi_model);
        output.has_bound = true;
        break;
    }
#endif
#ifdef KNITRO_FOUND
    case SolverName::Knitro: {
        KnitroContext knitro;
        if (!parameters.verbose)
            reduce_printout(knitro);
        if (has_time_limit)
            set_time_limit(knitro, time_limit);
        load(knitro, model);
        solve(knitro);
        if (!is_infeasible(knitro)) {
            output.solution = get_solution(knitro);
            output.objective_value = get_solution_value(knitro);
            output.found_solution = true;
        }
        break;
    }
#endif
#ifdef IPOPT_FOUND
    case SolverName::Ipopt: {
        IpoptModel ipopt_model;
        if (!parameters.verbose)
            reduce_printout(ipopt_model);
        if (has_time_limit)
            set_time_limit(ipopt_model, time_limit);
        load(ipopt_model, model);
        solve(ipopt_model);
        if (!is_infeasible(ipopt_model)) {
            output.solution = get_solution(ipopt_model);
            output.found_solution = !output.solution.empty();
        }
        break;
    }
#endif
#ifdef SCIP_FOUND
    case SolverName::Scip: {
        ScipModel scip_model;
        if (!parameters.verbose)
            reduce_printout(scip_model);
        if (has_time_limit)
            set_time_limit(scip_model, time_limit);
        load(scip_model, model);
        solve(scip_model);
        output.solution = get_solution(scip_model);
        output.found_solution = !output.solution.empty();
        output.bound = get_bound(scip_model);
        output.has_bound = true;
        break;
    }
#endif
    default:
        throw std::logic_error(FUNC_SIGNATURE + ": unreachable.");
    }

    if (output.found_solution && solver != SolverName::Knitro)
        output.objective_value = model.evaluate_objective(output.solution);

    return output;
}

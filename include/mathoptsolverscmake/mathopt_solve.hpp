#pragma once

/**
 * Solve a model with any of the solvers available in the build.
 *
 * Which solvers are available depends on the MATHOPTSOLVERSCMAKE_USE_* CMake
 * options the library was built with; see available_solvers().
 */

#include "mathoptsolverscmake/mathopt.hpp"

#include <limits>
#include <vector>

namespace mathoptsolverscmake
{

struct SolveParameters
{
    /** Solver: Cbc, Highs, Xpress, Cplex, Gurobi, Knitro, Ipopt or Scip. */
    SolverName solver = SolverName::Highs;

    /** Time limit, in seconds; infinity for no limit. */
    double time_limit = std::numeric_limits<double>::infinity();

    /** If false, the solver's log is not written to the standard output. */
    bool verbose = true;
};

struct SolveOutput
{
    /** True if the solver found a feasible solution. */
    bool found_solution = false;

    /** The solution found (empty if none). */
    std::vector<double> solution;

    /** Objective value of the solution (NaN if none). */
    double objective_value = std::numeric_limits<double>::quiet_NaN();

    /** True if the solver gives a bound (MILP solvers). */
    bool has_bound = false;

    /** Bound on the optimal objective value (NaN if none). */
    double bound = std::numeric_limits<double>::quiet_NaN();
};

/** Return true if the solver can be used by solve() in this build. */
bool is_available(SolverName solver);

/** The solvers solve() can use in this build. */
std::vector<SolverName> available_solvers();

/**
 * Solve a model.
 *
 * Throws std::invalid_argument if the solver isn't available in this build
 * or doesn't support the model (for example HiGHS with a nonlinear model).
 */
SolveOutput solve(
        const MathOptModel& model,
        const SolveParameters& parameters = SolveParameters());

}

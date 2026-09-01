#pragma once

#include "mathoptsolverscmake/mathopt.hpp"

#include <vector>

namespace mathoptsolverscmake
{

struct SmsppOutput
{
    double objective_value = 0;

    std::vector<double> solution;
};

/**
 * Solve a box-constrained black-box model with SMS++'s BundleSolver.
 *
 * model.objective_function is wrapped into a custom SMS++ C05Function (an
 * oracle able to return a value and a single "diagonal" linearization, i.e.
 * a gradient, at any point), plugged as the Objective of an AbstractBlock
 * whose only Variable are model's, bounded with BoxConstraint. A BundleSolver
 * is then attached to the Block and run to convergence.
 *
 * model must be box-constrained (see MathOptModel::is_box_constrained()) and
 * have model.objective_function set. model.objective_function must define a
 * convex function when model.objective_direction is Minimize, or a concave
 * one when it is Maximize (BundleSolver only accepts one or the other; this
 * cannot be checked and is the caller's responsibility -- the primary
 * intended use is Lagrangian relaxations, which always satisfy this).
 * Every variable's lower bound must be exactly 0 or -infinity (the
 * Lagrangian-multiplier convention that BundleSolver's box-constraint
 * handling is built around); upper bounds are unrestricted.
 */
SmsppOutput solve_smspp(
        const MathOptModel& model);

}

#pragma once

#include "mathoptsolverscmake/mathopt.hpp"

#include <functional>
#include <vector>

namespace mathoptsolverscmake
{

/**
 * A composable builder for a node of a nonlinear expression tree (see the
 * "Nonlinear structures" comment on MathOptModel in mathopt.hpp for the
 * on-disk layout: operators/values/variables/parent stored as parallel
 * arrays in DFS pre-order).
 *
 * Given the parallel arrays to append to and the index of its parent
 * (-1 for a root), a NonlinearExpression appends itself and then
 * recursively its operands, so the resulting arrays are already in the
 * required pre-order without the caller having to track indices by hand.
 *
 * Build one with the nonlinear_* functions below, then install it into a
 * MathOptModel with set_objective_nonlinear_expression() or
 * add_constraint_nonlinear_expression().
 */
using NonlinearExpression = std::function<void(
        std::vector<char>& operators,
        std::vector<double>& values,
        std::vector<int>& variables,
        std::vector<int>& parents,
        int parent_id)>;

/** A single variable, by index. */
NonlinearExpression nonlinear_variable(int variable_id);
/** A constant value. */
NonlinearExpression nonlinear_constant(double value);

NonlinearExpression nonlinear_negate(NonlinearExpression operand);
NonlinearExpression nonlinear_exp(NonlinearExpression operand);
NonlinearExpression nonlinear_log(NonlinearExpression operand);
NonlinearExpression nonlinear_sqrt(NonlinearExpression operand);
NonlinearExpression nonlinear_sin(NonlinearExpression operand);
NonlinearExpression nonlinear_cos(NonlinearExpression operand);
NonlinearExpression nonlinear_tan(NonlinearExpression operand);

/** n-ary sum (n >= 2 terms). */
NonlinearExpression nonlinear_sum(std::vector<NonlinearExpression> terms);
/** n-ary product (n >= 2 factors). */
NonlinearExpression nonlinear_product(std::vector<NonlinearExpression> factors);
/** Binary: left - right. */
NonlinearExpression nonlinear_subtract(NonlinearExpression left, NonlinearExpression right);
/** Binary: left / right. */
NonlinearExpression nonlinear_divide(NonlinearExpression left, NonlinearExpression right);
/** Binary: base ^ exponent. */
NonlinearExpression nonlinear_power(NonlinearExpression base, NonlinearExpression exponent);

/** Install a built expression as the model's nonlinear objective term. */
void set_objective_nonlinear_expression(
        MathOptModel& model,
        const NonlinearExpression& expression);

/**
 * Install a built expression as a constraint's nonlinear term. Grows
 * model.nonlinear_elements_constraints_starts as needed (the MathOptModel
 * constructor does not size it, unlike its linear/quadratic counterparts).
 */
void add_constraint_nonlinear_expression(
        MathOptModel& model,
        int constraint_id,
        const NonlinearExpression& expression);

}

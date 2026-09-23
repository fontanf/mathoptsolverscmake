#pragma once

#include "mathoptsolverscmake/mathopt.hpp"

namespace mathoptsolverscmake
{

/**
 * Build the Lagrangian dual of a linear program.
 *
 * Only LPs are supported (throws std::invalid_argument unless
 * model.is_lp()).
 *
 * Writing the primal as
 *
 *     min / max  c^T x
 *     s.t.       L <= A x <= U
 *                l <=   x <= u
 *
 * the dual has:
 * - one variable y_i per finite side of each constraint i (two for a Range
 *   constraint, none for a Free one, a single free one for an Equality);
 * - one variable z_j per finite side of each variable bound j, except that
 *   bounds equal to 0 get none (their z_j would have a zero objective
 *   coefficient, so they are eliminated into the dual constraint's sides
 *   instead; e.g. x >= 0 yields A^T y <= c rather than A^T y + z = c,
 *   z >= 0);
 * - one constraint per primal variable j: (A^T y)_j + sum z_j = c_j, turned
 *   into an inequality (or a Free constraint) by the eliminated z_j above;
 * - objective sum_i (L_i or U_i) y_i + sum_j (l_j or u_j) z_j.
 *
 * The dual of a minimization problem is a maximization problem and vice
 * versa. Signs follow the usual solver convention for dual values (the
 * derivative of the optimal objective value with respect to the bound):
 * for a minimization primal, a lower-side dual variable is >= 0 and an
 * upper-side one is <= 0; for a maximization primal, it is the opposite.
 * With that convention, the primal and dual constraint matrices are exactly
 * transposed of each other, and at optimality both objective values are
 * equal (no constant offset is needed).
 *
 * Dual variables are ordered as: all the y's, in primal constraint order
 * (lower side before upper side for a Range constraint), then all the z's,
 * in primal variable order (lower bound before upper bound). They are named
 * after the primal constraint/variable, suffixed with "_lb"/"_ub" (the side
 * they relax, or "_eq" for an Equality constraint / fixed variable). Dual
 * constraint j is named after primal variable j.
 */
MathOptModel dual(const MathOptModel& model);

}

#pragma once

/**
 * Statistics of a model: its numbers of variables (by type), constraints,
 * nonzeros...
 *
 * They are computed from a model (compute_statistics()), or read from a file
 * (read_statistics()): for an .nl file, from its header only, which is fast
 * even for very large files, and works for any .nl file, including the ones
 * read_nl() doesn't support.
 */

#include "mathoptsolverscmake/mathopt.hpp"

#include <string>

namespace mathoptsolverscmake
{

struct MathOptModelStatistics
{
    int number_of_variables = 0;

    int number_of_continuous_variables = 0;

    /**
     * Binary variables: Binary ones, and, as on the MIPLIB statistics page,
     * Integer ones bounded to [0, 1]. From an .nl header, only the binary
     * variables appearing linearly: the others are counted as integer.
     */
    int number_of_binary_variables = 0;

    /** Integer variables (binary ones excluded). */
    int number_of_integer_variables = 0;

    int number_of_semi_continuous_variables = 0;

    int number_of_semi_integer_variables = 0;

    int number_of_constraints = 0;

    /**
     * Nonzeros of the constraints' linear terms. From an .nl header, the
     * nonzeros of the constraints' Jacobian, nonlinear terms included.
     */
    long number_of_elements = 0;

    /**
     * Constraints with quadratic terms but no other nonlinear terms. From an
     * .nl header, which doesn't distinguish them, they are counted as
     * nonlinear.
     */
    int number_of_quadratic_constraints = 0;

    /** Constraints with nonlinear (expression tree or black-box) terms. */
    int number_of_nonlinear_constraints = 0;

    /** Objectives with at least one term. */
    int number_of_objectives = 0;

    int number_of_quadratic_objectives = 0;

    int number_of_nonlinear_objectives = 0;

    /** Common subexpressions. */
    int number_of_defined_variables = 0;

    /** Logical constraints (.nl only; read_nl() doesn't support them). */
    int number_of_logical_constraints = 0;

    /** Complementarity constraints (.nl only; read_nl() doesn't support them). */
    int number_of_complementarity_constraints = 0;
};

/** Compute the statistics of a model. */
MathOptModelStatistics compute_statistics(const MathOptModel& model);

/**
 * Read the statistics of an (uncompressed) .mps, .lp, .json or .nl file,
 * guessing the format from its extension: from its header for an .nl file
 * (see read_nl_statistics()), from the model read otherwise.
 */
MathOptModelStatistics read_statistics(const std::string& file_path);

}

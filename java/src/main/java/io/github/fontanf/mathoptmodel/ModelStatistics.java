package io.github.fontanf.mathoptmodel;

/**
 * Statistics of a model: numbers of variables (by type), constraints, nonzeros... (see
 * MathOptModelStatistics in mathopt_statistics.hpp). From {@link Model#statistics()}, or
 * {@link Model#readStatistics(java.nio.file.Path)}, which reads only the header of an .nl file.
 *
 * @param numberOfBinaryVariables binary variables, and integer ones bounded to [0, 1]; from an
 *     .nl header, only the binary variables appearing linearly (the others count as integer)
 * @param numberOfIntegerVariables integer variables, binary ones excluded
 * @param numberOfElements nonzeros of the constraints' linear terms; from an .nl header, of the
 *     constraints' Jacobian
 * @param numberOfQuadraticConstraints constraints with quadratic terms but no other nonlinear
 *     terms (counted as nonlinear from an .nl header)
 * @param numberOfObjectives objectives with at least one term
 * @param numberOfDefinedVariables common subexpressions
 * @param numberOfLogicalConstraints logical constraints (.nl only)
 * @param numberOfComplementarityConstraints complementarity constraints (.nl only)
 */
public record ModelStatistics(
        int numberOfVariables,
        int numberOfContinuousVariables,
        int numberOfBinaryVariables,
        int numberOfIntegerVariables,
        int numberOfSemiContinuousVariables,
        int numberOfSemiIntegerVariables,
        int numberOfConstraints,
        long numberOfElements,
        int numberOfQuadraticConstraints,
        int numberOfNonlinearConstraints,
        int numberOfObjectives,
        int numberOfQuadraticObjectives,
        int numberOfNonlinearObjectives,
        int numberOfDefinedVariables,
        int numberOfLogicalConstraints,
        int numberOfComplementarityConstraints) {

    /** From the values filled by the C API, indexed by the MOM_STATISTICS_* constants. */
    static ModelStatistics fromValues(long[] v) {
        return new ModelStatistics((int) v[0], (int) v[1], (int) v[2], (int) v[3], (int) v[4], (int) v[5], (int) v[6], v[7], (int) v[8], (int) v[9], (int) v[10], (int) v[11], (int) v[12], (int) v[13], (int) v[14], (int) v[15]);
    }
}

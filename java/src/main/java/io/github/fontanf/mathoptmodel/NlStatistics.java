package io.github.fontanf.mathoptmodel;

/**
 * Statistics of an .nl file, from its header (see {@link Model#readNlStatistics(java.nio.file.Path)}).
 *
 * @param binary true for the binary encoding, false for the text one
 * @param numberOfNonlinearVariables variables appearing nonlinearly in a constraint or an objective
 * @param numberOfIntegerVariables integer variables, binary ones excluded
 * @param numberOfElements nonzeros of the constraints' Jacobian
 * @param numberOfObjectiveElements nonzeros of the objectives' gradients
 * @param numberOfDefinedVariables common subexpressions
 */
public record NlStatistics(
        boolean binary,
        int numberOfVariables,
        int numberOfConstraints,
        int numberOfObjectives,
        int numberOfNonlinearConstraints,
        int numberOfNonlinearObjectives,
        int numberOfNonlinearVariables,
        int numberOfBinaryVariables,
        int numberOfIntegerVariables,
        long numberOfElements,
        long numberOfObjectiveElements,
        int numberOfDefinedVariables,
        int numberOfLogicalConstraints,
        int numberOfComplementarityConstraints) {

    /** From the values filled by mom_read_nl_statistics, in the order of the MOM_NL_* constants. */
    static NlStatistics fromValues(long[] v) {
        return new NlStatistics(v[0] != 0, (int) v[1], (int) v[2], (int) v[3], (int) v[4], (int) v[5],
                (int) v[6], (int) v[7], (int) v[8], v[9], v[10], (int) v[11], (int) v[12], (int) v[13]);
    }
}

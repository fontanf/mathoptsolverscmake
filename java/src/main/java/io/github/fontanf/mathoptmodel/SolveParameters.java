package io.github.fontanf.mathoptmodel;

/**
 * Parameters of {@link Model#solve(Solver, SolveParameters)}: by default, no
 * time limit, and the solver's log is written to the standard output.
 */
public final class SolveParameters {

    double timeLimit = Double.POSITIVE_INFINITY;
    boolean verbose = true;

    /** Time limit, in seconds ({@code Double.POSITIVE_INFINITY}: none). */
    public SolveParameters timeLimit(double seconds) {
        if (Double.isNaN(seconds) || seconds < 0)
            throw new IllegalArgumentException("The time limit must be non-negative, not " + seconds + ".");
        this.timeLimit = seconds;
        return this;
    }

    /** If false, the solver's log is not written to the standard output. */
    public SolveParameters verbose(boolean verbose) {
        this.verbose = verbose;
        return this;
    }
}

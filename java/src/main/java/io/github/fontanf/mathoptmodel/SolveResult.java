package io.github.fontanf.mathoptmodel;

import java.util.Locale;

/**
 * The result of {@link Model#solve(Solver)}: the solution found, if any, its
 * objective value, and the bound on the optimal objective value given by MILP
 * solvers.
 */
public final class SolveResult {

    private final boolean foundSolution;
    private final double[] solution;
    private final double objectiveValue;
    private final boolean hasBound;
    private final double bound;

    SolveResult(boolean foundSolution, double[] solution, double objectiveValue, boolean hasBound, double bound) {
        this.foundSolution = foundSolution;
        this.solution = solution;
        this.objectiveValue = objectiveValue;
        this.hasBound = hasBound;
        this.bound = bound;
    }

    /** True if the solver found a feasible solution. */
    public boolean foundSolution() {
        return foundSolution;
    }

    /** A copy of the solution found, one value per variable (empty if none). */
    public double[] solution() {
        return solution.clone();
    }

    /** Objective value of the solution (NaN if none). */
    public double objectiveValue() {
        return objectiveValue;
    }

    /** True if the solver gives a bound (MILP solvers). */
    public boolean hasBound() {
        return hasBound;
    }

    /** Bound on the optimal objective value (NaN if none). */
    public double bound() {
        return bound;
    }

    @Override
    public String toString() {
        String text = foundSolution
                ? String.format(Locale.ROOT, "objective value %s", objectiveValue)
                : "no solution found";
        if (hasBound)
            text += String.format(Locale.ROOT, ", bound %s", bound);
        return "SolveResult[" + text + "]";
    }
}

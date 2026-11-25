package io.github.fontanf.mathoptmodel;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * A solver, for {@link Model#solve(Solver)} (MathOptModel's SolverName; the
 * ordinals match). Which solvers can be used depends on the native library:
 * see {@link #isAvailable()}.
 */
public enum Solver {
    CBC,
    HIGHS,
    XPRESS,
    CPLEX,
    GUROBI,
    KNITRO,
    IPOPT,
    SCIP;

    /** True if the native library can solve with this solver. */
    public boolean isAvailable() {
        try {
            return (int) Native.SOLVER_AVAILABLE.invoke(ordinal()) != 0;
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /** The solvers the native library can solve with. */
    public static List<Solver> available() {
        List<Solver> solvers = new ArrayList<>();
        for (Solver solver: values())
            if (solver.isAvailable())
                solvers.add(solver);
        return Collections.unmodifiableList(solvers);
    }
}

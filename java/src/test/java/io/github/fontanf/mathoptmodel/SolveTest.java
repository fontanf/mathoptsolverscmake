package io.github.fontanf.mathoptmodel;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assumptions.assumeFalse;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

class SolveTest {

    private static final double INF = Double.POSITIVE_INFINITY;

    private static final Path P0033 = Path.of("..", "data", "miplib3", "p0033.mps");

    /** max x + 2 y + 3 z s.t. x + y + z <= 10, x - y >= -3; optimal: x = 4, y = 1, z = 5, objective 21. */
    private static Model smallMilp() {
        Model model = new Model();
        model.addVariables(new Variables(3)
                .lowerBounds(0.0)
                .upperBounds(new double[] {10.0, 1.0, 5.0})
                .types(new VariableType[] {VariableType.CONTINUOUS, VariableType.BINARY, VariableType.INTEGER}));
        model.setObjective(new Objective()
                .direction(ObjectiveDirection.MAXIMIZE)
                .coefficients(new double[] {1.0, 2.0, 3.0}));
        model.addConstraints(new Constraints(2)
                .lowerBounds(new double[] {-INF, -3.0})
                .upperBounds(new double[] {10.0, INF})
                .linear(new int[] {0, 3}, new int[] {0, 1, 2, 0, 1}, new double[] {1, 1, 1, 1, -1}));
        return model;
    }

    private static SolveParameters quiet() {
        return new SolveParameters().verbose(false);
    }

    @Test
    void availableSolvers() {
        List<Solver> available = Solver.available();
        for (Solver solver: Solver.values())
            assertEquals(available.contains(solver), solver.isAvailable(), solver.name());
    }

    @Test
    void solveMilp() {
        assumeTrue(Solver.HIGHS.isAvailable());
        try (Model model = smallMilp()) {
            SolveResult result = model.solve(Solver.HIGHS, quiet());
            assertTrue(result.foundSolution());
            assertEquals(21, result.objectiveValue(), 1e-6);
            assertArrayEquals(new double[] {4, 1, 5}, result.solution(), 1e-6);
            assertTrue(result.hasBound());
            assertEquals(21, result.bound(), 1e-6);
            assertTrue(result.toString().contains("objective value 21"), result.toString());
        }
    }

    @Test
    void solveInfeasible() {
        assumeTrue(Solver.HIGHS.isAvailable());
        try (Model model = new Model()) {
            model.addVariables(new Variables(1));
            model.addConstraints(new Constraints(2)
                    .lowerBounds(new double[] {1.0, -INF})
                    .upperBounds(new double[] {INF, 0.0})
                    .linear(new int[] {0, 1}, new int[] {0, 0}, new double[] {1, 1}));
            SolveResult result = model.solve(Solver.HIGHS, quiet());
            assertFalse(result.foundSolution());
            assertEquals(0, result.solution().length);
            assertTrue(Double.isNaN(result.objectiveValue()));
        }
    }

    @Test
    void solveFile() throws IOException {
        assumeTrue(Solver.HIGHS.isAvailable());
        try (Model model = Model.read(P0033)) {
            SolveResult result = model.solve(Solver.HIGHS, quiet().timeLimit(60));
            assertEquals(3089, result.objectiveValue(), 1e-6);
            assertEquals(3089, result.bound(), 1e-6);
            assertTrue(model.checkSolution(result.solution()));
        }
    }

    @Test
    void invalidParameters() {
        assertThrows(IllegalArgumentException.class, () -> new SolveParameters().timeLimit(-1));
        assertThrows(IllegalArgumentException.class, () -> new SolveParameters().timeLimit(Double.NaN));
    }

    @Test
    void unavailableSolver() {
        List<Solver> unavailable = java.util.Arrays.stream(Solver.values()).filter(s -> !s.isAvailable()).toList();
        assumeFalse(unavailable.isEmpty());
        try (Model model = smallMilp()) {
            IllegalArgumentException e = assertThrows(IllegalArgumentException.class,
                    () -> model.solve(unavailable.get(0), quiet()));
            assertTrue(e.getMessage().contains("not available in this build"), e.getMessage());
        }
    }

    @Test
    void nonlinearModelRejected() {
        assumeTrue(Solver.HIGHS.isAvailable());
        try (Model model = smallMilp()) {
            model.setObjective(new Objective().quadratic(new int[] {0}, new int[] {0}, new double[] {1}));
            IllegalArgumentException e = assertThrows(IllegalArgumentException.class,
                    () -> model.solve(Solver.HIGHS, quiet()));
            assertTrue(e.getMessage().contains("MILP"), e.getMessage());
        }
    }
}

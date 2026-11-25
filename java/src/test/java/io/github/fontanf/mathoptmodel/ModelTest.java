package io.github.fontanf.mathoptmodel;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Random;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class ModelTest {

    private static final double INF = Double.POSITIVE_INFINITY;

    private static final Path P0033 = Path.of("..", "data", "miplib3", "p0033.mps");

    private static final String LP = """
            Maximize
             obj: x + 2 y + 3 z
            Subject To
             c1: x + y + z <= 10
             c2: x - y >= -3
            Bounds
             0 <= z <= 5
            General
             z
            Binary
             x
            End
            """;

    /** max x + 2 y + 3 z s.t. x + y + z <= 10, x - y >= -3, x in [0, 10], y binary, z integer in [0, 5]. */
    private static Model smallMilp() {
        Model model = new Model();
        int first = model.addVariables(new Variables(3)
                .lowerBounds(0.0)
                .upperBounds(new double[] {10.0, 1.0, 5.0})
                .types(new VariableType[] {VariableType.CONTINUOUS, VariableType.BINARY, VariableType.INTEGER})
                .names("x", "y", "z"));
        assertEquals(0, first);
        model.setObjective(new Objective()
                .direction(ObjectiveDirection.MAXIMIZE)
                .coefficients(new double[] {1.0, 2.0, 3.0}));
        first = model.addConstraints(new Constraints(2)
                .lowerBounds(new double[] {-INF, -3.0})
                .upperBounds(new double[] {10.0, INF})
                .linear(new int[] {0, 3}, new int[] {0, 1, 2, 0, 1}, new double[] {1, 1, 1, 1, -1})
                .names("c1", "c2"));
        assertEquals(0, first);
        return model;
    }

    @Test
    void emptyModel() {
        try (Model model = new Model()) {
            assertEquals(0, model.numberOfVariables());
            assertEquals(0, model.numberOfConstraints());
            assertEquals(ObjectiveDirection.MINIMIZE, model.objectiveDirection());
            assertTrue(model.check());
            assertEquals("Model[0 variables, 0 constraints]", model.toString());
        }
    }

    @Test
    void readMps() throws IOException {
        try (Model model = Model.read(P0033)) {
            assertEquals(33, model.numberOfVariables());
            assertEquals(16, model.numberOfConstraints());
            assertTrue(model.check());
            assertTrue(model.isMilp());
            assertFalse(model.isLp());
            assertEquals(33, model.variablesNames().size());
            assertEquals(16, model.constraintsClasses().length);
        }
    }

    @Test
    void writeFileMatchesString(@TempDir Path directory) throws IOException {
        // Uncompressed files are written by the native library directly,
        // compressed ones through writeString: both must give the same contents.
        try (Model model = Model.read(LP.getBytes(), Format.LP)) {
            for (Format format : Format.values()) {
                Path path = directory.resolve("model é." + format.extension());
                model.write(path);
                assertEquals(model.writeString(format), Files.readString(path), format.toString());
            }
            assertThrows(IOException.class, () -> model.write(directory.resolve("missing").resolve("model.lp")));
        }
    }

    @Test
    void readGzipAndWriteRoundTrip(@TempDir Path directory) throws IOException {
        try (Model model = Model.read(LP.getBytes(), Format.LP)) {
            Path path = directory.resolve("model.mps.gz");
            model.write(path);
            try (Model read = Model.read(path)) {
                assertEquals(model.writeString(Format.LP), read.writeString(Format.LP));
            }
        }
    }

    @Test
    void readNlNames(@TempDir Path directory) throws IOException {
        try (Model model = Model.read(LP.getBytes(), Format.LP)) {
            Path path = directory.resolve("model.nl");
            model.write(path);
            try (Model read = Model.read(path)) {
                assertEquals(List.of(), read.variablesNames());
            }
            Files.writeString(directory.resolve("model.col"), "a\nb\nc\n");
            Files.writeString(directory.resolve("model.row"), "first\nsecond\nobj\n");
            try (Model read = Model.read(path)) {
                assertEquals(List.of("a", "b", "c"), read.variablesNames());
                assertEquals(List.of("first", "second"), read.constraintsNames());
            }
            // A names file that doesn't match the model is ignored.
            Files.writeString(directory.resolve("model.col"), "a\nb\n");
            try (Model read = Model.read(path)) {
                assertEquals(List.of(), read.variablesNames());
                assertEquals(List.of("first", "second"), read.constraintsNames());
            }
        }
    }

    @Test
    void statistics(@TempDir Path directory) throws IOException {
        try (Model model = Model.read(LP.getBytes(), Format.LP)) {
            ModelStatistics statistics = model.statistics();
            assertEquals(3, statistics.numberOfVariables());
            assertEquals(2, statistics.numberOfConstraints());
            // .nl files' statistics come from their header, the others' from the model.
            for (String name : List.of("model.mps", "model.lp.gz", "model.json", "model.nl", "model.nl.gz")) {
                Path path = directory.resolve(name);
                model.write(path);
                ModelStatistics read = Model.readStatistics(path);
                assertEquals(model.numberOfVariables(), read.numberOfVariables(), name);
                assertEquals(model.numberOfConstraints(), read.numberOfConstraints(), name);
                assertEquals(model.numberOfElements(), read.numberOfElements(), name);
                assertEquals(1, read.numberOfContinuousVariables(), name);
                assertEquals(1, read.numberOfBinaryVariables(), name);
                assertEquals(1, read.numberOfIntegerVariables(), name);
                assertEquals(1, read.numberOfObjectives(), name);
            }
        }
    }

    @Test
    void readErrors() {
        assertThrows(MathOptModelException.class, () -> Model.read("not an nl file".getBytes(), Format.NL));
        assertThrows(IllegalArgumentException.class, () -> Model.read(Path.of("model.txt")));
    }

    @Test
    void buildMilpMatchesLpReader() {
        try (Model built = smallMilp(); Model read = Model.read(LP.getBytes(), Format.LP)) {
            assertTrue(built.check());
            assertEquals(ObjectiveDirection.MAXIMIZE, built.objectiveDirection());
            assertEquals(List.of("x", "y", "z"), built.variablesNames());
            assertEquals(List.of("c1", "c2"), built.constraintsNames());
            assertArrayEquals(new double[] {0, 0, 0}, built.variablesLowerBounds());
            assertArrayEquals(new double[] {10, 1, 5}, built.variablesUpperBounds());
            assertArrayEquals(new int[] {0, 3}, built.constraintsStarts());
            assertArrayEquals(new int[] {0, 1, 2, 0, 1}, built.elementsVariables());
            assertArrayEquals(new double[] {1, 1, 1, 1, -1}, built.elementsCoefficients());
            assertEquals(9.0, built.evaluateObjective(new double[] {1, 1, 2}));
            assertTrue(built.checkSolution(new double[] {1, 1, 2}));
            assertFalse(built.checkSolution(new double[] {0, 1, 5.5}));
            // The LP file's model, written back, matches the built one but
            // for the types (the LP above has x binary and y continuous).
            assertEquals(ConstraintSense.LESS_THAN_OR_EQUAL_TO, built.constraintSense(0));
            assertEquals(ConstraintSense.GREATER_THAN_OR_EQUAL_TO, built.constraintSense(1));
            assertEquals(read.numberOfConstraints(), built.numberOfConstraints());
            try (Model reread = Model.read(built.writeString(Format.LP).getBytes(), Format.LP)) {
                assertArrayEquals(built.variablesTypes(), reread.variablesTypes());
                assertArrayEquals(built.objectiveCoefficients(), reread.objectiveCoefficients());
                assertArrayEquals(built.constraintsLowerBounds(), reread.constraintsLowerBounds());
                assertArrayEquals(built.constraintsUpperBounds(), reread.constraintsUpperBounds());
            }
        }
    }

    @Test
    void defaultsAndScalars() {
        try (Model model = new Model()) {
            assertEquals(0, model.addVariables(new Variables(3)));
            assertArrayEquals(new double[] {-INF, -INF, -INF}, model.variablesLowerBounds());
            assertArrayEquals(new double[] {INF, INF, INF}, model.variablesUpperBounds());
            assertEquals(3, model.addVariables(new Variables(2)
                    .lowerBounds(0).upperBounds(1).types(VariableType.BINARY).objectiveCoefficients(4)));
            assertArrayEquals(new byte[] {0, 0, 0, 1, 1}, model.variablesTypesCodes());
            assertArrayEquals(new double[] {0, 0, 0, 4, 4}, model.objectiveCoefficients());
            assertEquals(0, model.addConstraints(new Constraints(2).upperBounds(1.0)));
            assertArrayEquals(new int[] {0, 0}, model.constraintsStarts());
            assertTrue(model.check());
        }
    }

    @Test
    void batchesAreAppended() {
        try (Model model = new Model()) {
            model.addVariables(new Variables(2).lowerBounds(0));
            model.addConstraints(new Constraints(1).upperBounds(1)
                    .linear(new int[] {0}, new int[] {0, 1}, new double[] {1, 2}));
            assertEquals(2, model.addVariables(new Variables(1).lowerBounds(0)));
            assertEquals(1, model.addConstraints(new Constraints(2).upperBounds(new double[] {5, 6})
                    .linear(new int[] {0, 1}, new int[] {2, 0, 2}, new double[] {3, 4, 5})));
            assertArrayEquals(new int[] {0, 2, 3}, model.constraintsStarts());
            assertTrue(model.check());
            double[] solution = {1, 1, 1};
            assertEquals(3.0, model.evaluateConstraint(solution, 0));
            assertEquals(3.0, model.evaluateConstraint(solution, 1));
            assertEquals(9.0, model.evaluateConstraint(solution, 2));
        }
    }

    @Test
    void namesAreCompleted() {
        try (Model model = new Model()) {
            model.addVariables(new Variables(2));
            model.addVariables(new Variables(1).names("zéta"));
            model.addVariables(new Variables(1));
            assertEquals(List.of("x0", "x1", "zéta", "x3"), model.variablesNames());
            assertTrue(model.check());
        }
    }

    @Test
    void quadratic() {
        try (Model model = new Model()) {
            model.addVariables(new Variables(2).lowerBounds(0));
            model.setObjective(new Objective().quadratic(new int[] {0, 0}, new int[] {0, 1}, new double[] {1, 1}));
            model.addConstraints(new Constraints(2)
                    .lowerBounds(new double[] {-INF, 1}).upperBounds(new double[] {4, INF})
                    .linear(new int[] {0, 0}, new int[] {0, 1}, new double[] {1, 1})
                    .quadratic(new int[] {0, 2}, new int[] {0, 1}, new int[] {0, 1}, new double[] {1, 1}));
            assertTrue(model.hasQuadratic());
            assertTrue(model.check());
            assertEquals(3.0, model.evaluateObjective(new double[] {1, 2}));
            assertEquals(5.0, model.evaluateConstraint(new double[] {1, 2}, 0));
            assertEquals(3.0, model.evaluateConstraint(new double[] {1, 2}, 1));
            assertEquals(ConstraintClass.GENERAL_QUADRATIC, model.constraintClass(0));
            try (Model read = Model.read(model.writeString(Format.JSON).getBytes(), Format.JSON)) {
                assertEquals(5.0, read.evaluateConstraint(new double[] {1, 2}, 0));
            }
        }
    }

    @Test
    void nonlinearAndDefinedVariables() {
        try (Model model = new Model()) {
            model.addVariables(new Variables(2).lowerBounds(0));
            // exp(x0) + x1 ^ 2 <= 4, then a constraint without tree, then log(x0) <= 0.
            model.addConstraints(new Constraints(3)
                    .upperBounds(new double[] {4, 1, 0})
                    .nonlinear(new int[] {0, 6, 6}, new ExpressionTrees("+evpvk" + "lv",
                            new int[] {2, 1, 0, 2, 0, 0, 1, 0},
                            new int[] {-1, -1, 0, -1, 1, -1, -1, 0},
                            new double[] {0, 0, 0, 0, 0, 2, 0, 0})));
            assertTrue(model.hasNonlinear());
            assertTrue(model.check());
            double[] solution = {1, 2};
            assertEquals(Math.E + 4, model.evaluateConstraint(solution, 0), 1e-12);
            assertEquals(0.0, model.evaluateConstraint(solution, 1));
            assertEquals(0.0, model.evaluateConstraint(solution, 2));
            assertArrayEquals("+evpvklv".getBytes(), model.nonlinearElementsOperators());

            // d0 = x0 * x1, d1 = d0 + x0; minimize d1.
            assertEquals(0, model.addDefinedVariables(new DefinedVariables(new int[] {0, 3},
                    new ExpressionTrees("*vv+dv", new int[] {2, 0, 0, 2, 0, 0},
                            new int[] {-1, 0, 1, -1, 0, 0}, null))));
            model.setObjective(new Objective().nonlinear(
                    new ExpressionTrees("d", new int[] {0}, new int[] {1}, null)));
            assertTrue(model.check());
            assertEquals(8.0, model.evaluateObjective(new double[] {2, 3}));
            try (Model expanded = model.expandDefinedVariables()) {
                assertFalse(expanded.hasDefinedVariables());
                assertEquals(8.0, expanded.evaluateObjective(new double[] {2, 3}));
            }
            try (Model read = Model.read(model.writeString(Format.NL).getBytes(), Format.NL)) {
                assertEquals(8.0, read.evaluateObjective(new double[] {2, 3}));
            }
            model.setObjective(new Objective().nonlinear(new ExpressionTrees(new byte[0], new int[0], null, null)));
            assertEquals(0.0, model.evaluateObjective(new double[] {2, 3}));
        }
    }

    @Test
    void invalidBatchesLeaveTheModelUnchanged() {
        try (Model model = smallMilp()) {
            String before = model.writeString(Format.JSON);
            assertThrowsMessage("not a variable index", () -> model.addConstraints(new Constraints(1)
                    .linear(new int[] {0}, new int[] {3}, new double[] {1})));
            assertThrowsMessage("starts must begin at 0", () -> model.addConstraints(new Constraints(1)
                    .linear(new int[] {1}, new int[] {0}, new double[] {1})));
            assertThrowsMessage("missing 1 node", () -> model.addConstraints(new Constraints(1)
                    .nonlinear(new int[] {0}, new ExpressionTrees("+v", new int[] {2, 0}, new int[] {0, 0}, null))));
            assertThrowsMessage("single root", () -> model.addConstraints(new Constraints(1)
                    .nonlinear(new int[] {0}, new ExpressionTrees("vv", new int[] {0, 0}, new int[] {0, 0}, null))));
            assertThrowsMessage("unknown operator", () -> model.addConstraints(new Constraints(1)
                    .nonlinear(new int[] {0}, new ExpressionTrees("Z", new int[] {0}, null, null))));
            assertThrowsMessage("references variable -1", () -> model.addConstraints(new Constraints(1)
                    .nonlinear(new int[] {0}, new ExpressionTrees("v", new int[] {0}, null, null))));
            assertThrowsMessage("not a VariableType", () -> model.addVariables(new Variables(1)
                    .types(new byte[] {5})));
            assertThrowsMessage("only reference lower indices", () -> model.addDefinedVariables(
                    new DefinedVariables(new int[] {0}, new ExpressionTrees("d", new int[] {0}, new int[] {0}, null))));
            assertThrowsMessage("not a variable index", () -> model.setObjective(new Objective()
                    .direction(ObjectiveDirection.MINIMIZE)
                    .quadratic(new int[] {0}, new int[] {7}, new double[] {1})));
            // Array lengths are checked in Java.
            assertThrowsMessage("has 2 entries", () -> new Constraints(1).upperBounds(new double[] {1, 2}));
            assertThrowsMessage("has 1 entries", () -> model.setObjective(new Objective().coefficients(new double[] {1})));
            assertThrowsMessage("has 1 entries", () -> model.evaluateObjective(new double[] {1}));
            assertEquals(before, model.writeString(Format.JSON));
            assertEquals(ObjectiveDirection.MAXIMIZE, model.objectiveDirection());
            assertTrue(model.check());
        }
    }

    @Test
    void closedModel() {
        Model model = new Model();
        model.close();
        model.close();
        assertThrows(IllegalStateException.class, model::numberOfVariables);
        assertEquals("Model[closed]", model.toString());
    }

    @Test
    void copyAndDual() {
        try (Model model = Model.read(LP.getBytes(), Format.LP); Model copy = model.copy()) {
            copy.addVariables(new Variables(1));
            assertEquals(3, model.numberOfVariables());
            assertEquals(4, copy.numberOfVariables());
        }
        try (Model lp = new Model()) {
            lp.addVariables(new Variables(2).lowerBounds(0).objectiveCoefficients(new double[] {1, 1}));
            lp.addConstraints(new Constraints(1).lowerBounds(1).linear(new int[] {0}, new int[] {0, 1}, new double[] {1, 2}));
            try (Model dual = lp.dual()) {
                assertEquals(1, dual.numberOfVariables());
            }
        }
    }

    @Test
    void largeBatch() {
        int n = 200_000;
        int m = 100_000;
        int k = 5;
        Random random = new Random(0);
        int[] starts = new int[m];
        int[] variables = new int[m * k];
        double[] coefficients = new double[m * k];
        for (int i = 0; i < m; ++i) {
            starts[i] = i * k;
            for (int j = 0; j < k; ++j) {
                variables[i * k + j] = (i * 7 + j * 40_000) % n;
                coefficients[i * k + j] = random.nextDouble();
            }
        }
        try (Model model = new Model()) {
            model.addVariables(new Variables(n).lowerBounds(0).upperBounds(1).types(VariableType.BINARY));
            model.addConstraints(new Constraints(m).upperBounds(1).linear(starts, variables, coefficients));
            assertEquals(m * k, model.numberOfElements());
            assertTrue(model.check());
            assertArrayEquals(variables, model.elementsVariables());
        }
    }

    private static void assertThrowsMessage(String expected, org.junit.jupiter.api.function.Executable executable) {
        IllegalArgumentException e = assertThrows(IllegalArgumentException.class, executable);
        assertTrue(e.getMessage().contains(expected), () -> "'" + e.getMessage() + "' doesn't contain '" + expected + "'");
    }

    static {
        // Fail early, with a clear message, if the test data is missing.
        if (!Files.exists(P0033))
            throw new IllegalStateException("Test data not found: " + P0033.toAbsolutePath());
    }
}

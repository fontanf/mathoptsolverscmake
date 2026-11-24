package io.github.fontanf.mathoptmodel;

import java.util.Arrays;
import java.util.List;
import java.util.Objects;

/**
 * A batch of constraints {@code lowerBound <= body <= upperBound}, for
 * {@link Model#addConstraints(Constraints)}.
 *
 * <p>The body is the sum of a linear part, a quadratic part and an expression
 * tree, each optional. The terms of all the constraints of the batch are
 * concatenated (CSR layout), and a {@code starts} array gives the offset of
 * each constraint's first term (one entry per constraint, the first one 0, a
 * constraint without terms having the same start as the next one).
 *
 * <p>Bounds are a single value, for every constraint, or an array with one
 * entry per constraint; they default to {@code -inf} and {@code +inf}. Arrays
 * are referenced, not copied.
 */
public final class Constraints {

    final int number;
    double[] lowerBounds;
    double lowerBound = Double.NEGATIVE_INFINITY;
    double[] upperBounds;
    double upperBound = Double.POSITIVE_INFINITY;
    int[] starts;
    int[] variables;
    double[] coefficients;
    int[] quadraticStarts;
    int[] quadraticVariables1;
    int[] quadraticVariables2;
    double[] quadraticCoefficients;
    int[] nonlinearStarts;
    ExpressionTrees nonlinear;
    List<String> names;

    /** A batch of {@code number} constraints. */
    public Constraints(int number) {
        if (number < 0)
            throw new IllegalArgumentException("number is negative.");
        this.number = number;
    }

    public int number() {
        return number;
    }

    public Constraints lowerBounds(double lowerBound) {
        this.lowerBounds = null;
        this.lowerBound = lowerBound;
        return this;
    }

    public Constraints lowerBounds(double[] lowerBounds) {
        Batches.checkLength(lowerBounds.length, number, "lowerBounds", "the batch");
        this.lowerBounds = lowerBounds;
        return this;
    }

    public Constraints upperBounds(double upperBound) {
        this.upperBounds = null;
        this.upperBound = upperBound;
        return this;
    }

    public Constraints upperBounds(double[] upperBounds) {
        Batches.checkLength(upperBounds.length, number, "upperBounds", "the batch");
        this.upperBounds = upperBounds;
        return this;
    }

    /** Linear terms: the constraint {@code i}'s terms are {@code variables[starts[i]..]} / {@code coefficients[starts[i]..]}. */
    public Constraints linear(int[] starts, int[] variables, double[] coefficients) {
        Batches.checkLength(starts.length, number, "starts", "the batch");
        Batches.checkLength(coefficients.length, variables.length, "coefficients", "variables");
        this.starts = starts;
        this.variables = variables;
        this.coefficients = coefficients;
        return this;
    }

    /** Quadratic terms {@code coefficient * x[variable1] * x[variable2]}, laid out like the linear ones. */
    public Constraints quadratic(int[] starts, int[] variables1, int[] variables2, double[] coefficients) {
        Batches.checkLength(starts.length, number, "quadraticStarts", "the batch");
        Batches.checkLength(variables2.length, variables1.length, "quadraticVariables2", "quadraticVariables1");
        Batches.checkLength(coefficients.length, variables1.length, "quadraticCoefficients", "quadraticVariables1");
        this.quadraticStarts = starts;
        this.quadraticVariables1 = variables1;
        this.quadraticVariables2 = variables2;
        this.quadraticCoefficients = coefficients;
        return this;
    }

    /** Expression trees: constraint {@code i}'s tree starts at node {@code starts[i]} (none if empty). */
    public Constraints nonlinear(int[] starts, ExpressionTrees trees) {
        Batches.checkLength(starts.length, number, "nonlinearStarts", "the batch");
        this.nonlinearStarts = starts;
        this.nonlinear = Objects.requireNonNull(trees, "trees");
        return this;
    }

    /**
     * Names. Once a batch has names, the model's unnamed constraints are named
     * "c" followed by their index.
     */
    public Constraints names(List<String> names) {
        Batches.checkLength(names.size(), number, "names", "the batch");
        this.names = names;
        return this;
    }

    public Constraints names(String... names) {
        return names(Arrays.asList(names));
    }
}

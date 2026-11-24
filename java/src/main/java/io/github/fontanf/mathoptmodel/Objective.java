package io.github.fontanf.mathoptmodel;

import java.util.Objects;

/**
 * An objective update, for {@link Model#setObjective(Objective)}: each part
 * set replaces the model's current one, the others are kept. Arrays are
 * referenced, not copied.
 */
public final class Objective {

    ObjectiveDirection direction;
    boolean setCoefficients;
    double[] coefficients;
    double coefficient;
    boolean setQuadratic;
    int[] quadraticVariables1;
    int[] quadraticVariables2;
    double[] quadraticCoefficients;
    ExpressionTrees nonlinear;

    public Objective direction(ObjectiveDirection direction) {
        this.direction = Objects.requireNonNull(direction, "direction");
        return this;
    }

    /** The same linear coefficient for every variable. */
    public Objective coefficients(double coefficient) {
        this.setCoefficients = true;
        this.coefficients = null;
        this.coefficient = coefficient;
        return this;
    }

    /** The linear coefficients, one per variable of the model. */
    public Objective coefficients(double[] coefficients) {
        this.setCoefficients = true;
        this.coefficients = Objects.requireNonNull(coefficients, "coefficients");
        return this;
    }

    /** Quadratic terms {@code coefficient * x[variable1] * x[variable2]} (empty arrays: none). */
    public Objective quadratic(int[] variables1, int[] variables2, double[] coefficients) {
        Batches.checkLength(variables2.length, variables1.length, "quadraticVariables2", "quadraticVariables1");
        Batches.checkLength(coefficients.length, variables1.length, "quadraticCoefficients", "quadraticVariables1");
        this.setQuadratic = true;
        this.quadraticVariables1 = variables1;
        this.quadraticVariables2 = variables2;
        this.quadraticCoefficients = coefficients;
        return this;
    }

    /** An expression tree (no nodes: none). */
    public Objective nonlinear(ExpressionTrees tree) {
        this.nonlinear = Objects.requireNonNull(tree, "tree");
        return this;
    }
}

package io.github.fontanf.mathoptmodel;

import java.util.Arrays;
import java.util.List;
import java.util.Objects;

/**
 * A batch of variables, for {@link Model#addVariables(Variables)}.
 *
 * <p>Each attribute is a single value, for every variable of the batch, or an
 * array with one entry per variable. By default, variables are continuous,
 * unbounded ({@code -inf} to {@code +inf}) and have a zero objective
 * coefficient. Arrays are referenced, not copied.
 */
public final class Variables {

    final int number;
    double[] lowerBounds;
    double lowerBound = Double.NEGATIVE_INFINITY;
    double[] upperBounds;
    double upperBound = Double.POSITIVE_INFINITY;
    byte[] types;
    VariableType type = VariableType.CONTINUOUS;
    double[] objectiveCoefficients;
    double objectiveCoefficient = 0.0;
    double[] initialValues;
    List<String> names;

    /** A batch of {@code number} variables. */
    public Variables(int number) {
        if (number < 0)
            throw new IllegalArgumentException("number is negative.");
        this.number = number;
    }

    public int number() {
        return number;
    }

    public Variables lowerBounds(double lowerBound) {
        this.lowerBounds = null;
        this.lowerBound = lowerBound;
        return this;
    }

    public Variables lowerBounds(double[] lowerBounds) {
        Batches.checkLength(lowerBounds.length, number, "lowerBounds", "the batch");
        this.lowerBounds = lowerBounds;
        return this;
    }

    public Variables upperBounds(double upperBound) {
        this.upperBounds = null;
        this.upperBound = upperBound;
        return this;
    }

    public Variables upperBounds(double[] upperBounds) {
        Batches.checkLength(upperBounds.length, number, "upperBounds", "the batch");
        this.upperBounds = upperBounds;
        return this;
    }

    public Variables types(VariableType type) {
        this.types = null;
        this.type = Objects.requireNonNull(type, "type");
        return this;
    }

    public Variables types(VariableType[] types) {
        Batches.checkLength(types.length, number, "types", "the batch");
        byte[] codes = new byte[types.length];
        for (int i = 0; i < types.length; ++i)
            codes[i] = Objects.requireNonNull(types[i], "types[" + i + "]").code();
        this.types = codes;
        return this;
    }

    /** Types as codes ({@link VariableType#code()}). */
    public Variables types(byte[] types) {
        Batches.checkLength(types.length, number, "types", "the batch");
        this.types = types;
        return this;
    }

    public Variables objectiveCoefficients(double objectiveCoefficient) {
        this.objectiveCoefficients = null;
        this.objectiveCoefficient = objectiveCoefficient;
        return this;
    }

    public Variables objectiveCoefficients(double[] objectiveCoefficients) {
        Batches.checkLength(objectiveCoefficients.length, number, "objectiveCoefficients", "the batch");
        this.objectiveCoefficients = objectiveCoefficients;
        return this;
    }

    public Variables initialValues(double[] initialValues) {
        Batches.checkLength(initialValues.length, number, "initialValues", "the batch");
        this.initialValues = initialValues;
        return this;
    }

    /**
     * Names. Once a batch has names, the model's unnamed variables are named
     * "x" followed by their index.
     */
    public Variables names(List<String> names) {
        Batches.checkLength(names.size(), number, "names", "the batch");
        this.names = names;
        return this;
    }

    public Variables names(String... names) {
        return names(Arrays.asList(names));
    }
}

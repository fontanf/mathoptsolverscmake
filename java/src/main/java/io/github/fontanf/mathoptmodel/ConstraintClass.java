package io.github.fontanf.mathoptmodel;

/**
 * Constraint class (MathOptModel's ConstraintClass; the ordinals match), based
 * on <a href="https://miplib.zib.de/statistics.html">MIPLIB's classification</a>.
 */
public enum ConstraintClass {
    /** Linear constraint with no variables. */
    EMPTY,
    /** Linear constraint with no finite side. */
    FREE,
    /** Linear constraint with a single variable. */
    SINGLETON,
    /** Linear constraint of the type ax+by=c. */
    AGGREGATION,
    /** Linear constraint of the type ax−ay≤b, where x and y must have the same type. */
    PRECEDENCE,
    /** Linear constraint of the form ax+by≤c, x∈{0,1}. */
    VARIABLE_BOUND,
    /** Linear constraint of the form ∑xi=1, xi∈{0,1}. */
    SET_PARTITIONING,
    /** Linear constraint of the form ∑xi≤1, xi∈{0,1}. */
    SET_PACKING,
    /** Linear constraint of the form ∑xi≥1, xi∈{0,1}. */
    SET_COVERING,
    /** Linear constraint of the form ∑xi=k, xi∈{0,1}, k≥2. */
    CARDINALITY,
    /** Linear constraint of the form ∑xi≤b, xi∈{0,1}, b∈N≥2. */
    INVARIANT_KNAPSACK,
    /** Linear constraint of the form ∑aixi=b, xi∈{0,1}, b∈N≥2. */
    EQUATION_KNAPSACK,
    /** Linear constraint of the form ∑aixi+ax≤a, x,xi∈{0,1}, a∈N≥2. */
    BINPACKING,
    /** Linear constraint of the form ∑akxk≤b, xi∈{0,1}, b∈N≥2. */
    KNAPSACK,
    /** Linear constraint of the form ∑akxk≤b, xi∈Z, b∈N. */
    INTEGER_KNAPSACK,
    /** Linear constraint of the form ∑akxk+∑pjsj{≤,=}b, xi∈{0,1}, sj continuous. */
    MIXED_BINARY,
    /** Linear constraint with no special structure. */
    GENERAL_LINEAR,
    /** Quadratic constraint with no special structure. */
    GENERAL_QUADRATIC,
    /** Black-box constraint with no special structure. */
    GENERAL_BLACK_BOX,
    /** Nonlinear constraint expressed as an expression tree. */
    GENERAL_NONLINEAR;

    static ConstraintClass fromCode(int code) {
        return values()[code];
    }
}

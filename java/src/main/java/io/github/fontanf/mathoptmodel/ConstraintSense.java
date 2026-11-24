package io.github.fontanf.mathoptmodel;

/** Constraint sense, from its bounds (MathOptModel's ConstraintSense; the ordinals match). */
public enum ConstraintSense {
    LESS_THAN_OR_EQUAL_TO,
    GREATER_THAN_OR_EQUAL_TO,
    EQUALITY,
    RANGE,
    FREE;

    static ConstraintSense fromCode(int code) {
        return values()[code];
    }
}

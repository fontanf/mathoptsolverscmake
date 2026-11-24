package io.github.fontanf.mathoptmodel;

/** Objective direction (MathOptModel's ObjectiveDirection; the ordinals match). */
public enum ObjectiveDirection {
    MINIMIZE,
    MAXIMIZE;

    static ObjectiveDirection fromCode(int code) {
        return values()[code];
    }
}

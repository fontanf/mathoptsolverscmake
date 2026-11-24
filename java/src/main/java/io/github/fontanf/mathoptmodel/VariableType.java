package io.github.fontanf.mathoptmodel;

/** Variable type (MathOptModel's VariableType; the ordinals match). */
public enum VariableType {
    CONTINUOUS,
    BINARY,
    INTEGER;

    /** The code of this type in byte arrays of variable types. */
    public byte code() {
        return (byte) ordinal();
    }

    /** The type of a code of a byte array of variable types. */
    public static VariableType fromCode(int code) {
        return values()[code];
    }
}

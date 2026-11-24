package io.github.fontanf.mathoptmodel;

/** An error reported by the native library, other than an invalid argument. */
public class MathOptModelException extends RuntimeException {

    public MathOptModelException(String message) {
        super(message);
    }

    public MathOptModelException(String message, Throwable cause) {
        super(message, cause);
    }
}

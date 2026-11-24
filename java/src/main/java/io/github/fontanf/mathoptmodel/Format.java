package io.github.fontanf.mathoptmodel;

import java.nio.file.Path;
import java.util.Locale;
import java.util.Optional;

/** Model file formats. */
public enum Format {
    /** Free-format MPS. */
    MPS("mps"),
    /** CPLEX-style LP (linear and mixed-integer linear models only). */
    LP("lp"),
    /** MathOptSolversCMake JSON. */
    JSON("json"),
    /** AMPL .nl (text when written; text or binary when read). */
    NL("nl");

    private final String extension;

    Format(String extension) {
        this.extension = extension;
    }

    /** The file extension of this format, without the dot. */
    public String extension() {
        return extension;
    }

    int code() {
        return ordinal();
    }

    /** The format of a file from its extension, ignoring a trailing .gz. */
    public static Optional<Format> fromPath(Path path) {
        String name = path.getFileName().toString().toLowerCase(Locale.ROOT);
        if (name.endsWith(".gz"))
            name = name.substring(0, name.length() - 3);
        for (Format format : values()) {
            if (name.endsWith("." + format.extension))
                return Optional.of(format);
        }
        return Optional.empty();
    }
}

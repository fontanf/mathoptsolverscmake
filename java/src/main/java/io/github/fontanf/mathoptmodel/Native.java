package io.github.fontanf.mathoptmodel;

import java.io.IOException;
import java.io.InputStream;
import java.io.UncheckedIOException;
import java.lang.foreign.Arena;
import java.lang.foreign.FunctionDescriptor;
import java.lang.foreign.Linker;
import java.lang.foreign.MemoryLayout;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.SymbolLookup;
import java.lang.invoke.MethodHandle;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.Locale;

import static java.lang.foreign.ValueLayout.ADDRESS;
import static java.lang.foreign.ValueLayout.JAVA_BYTE;
import static java.lang.foreign.ValueLayout.JAVA_DOUBLE;
import static java.lang.foreign.ValueLayout.JAVA_INT;
import static java.lang.foreign.ValueLayout.JAVA_LONG;

/**
 * Downcall handles to the C API (mathoptmodel_c.h) of the native
 * "mathoptmodel" library.
 *
 * <p>The library is loaded from, in this order: the path given by the
 * {@code mathoptmodel.library.path} system property; the copy bundled in the
 * jar for the current platform; the system's library search path.
 *
 * <p>Functions that only copy arrays (building, array getters, evaluations)
 * are linked as "critical" functions with heap access, so that Java arrays are
 * passed to them directly, without an intermediate off-heap copy.
 */
final class Native {

    private Native() {
    }

    /* Status codes. */
    static final int OK = 0;
    static final int INVALID_ARGUMENT = 1;

    /* Model arrays (see mathoptmodel_c.h). */
    static final int VARIABLES_LOWER_BOUNDS = 0;
    static final int VARIABLES_UPPER_BOUNDS = 1;
    static final int VARIABLES_TYPES = 2;
    static final int VARIABLES_INITIAL_VALUES = 3;
    static final int CONSTRAINTS_LOWER_BOUNDS = 4;
    static final int CONSTRAINTS_UPPER_BOUNDS = 5;
    static final int CONSTRAINTS_CLASSES = 6;
    static final int OBJECTIVE_COEFFICIENTS = 7;
    static final int CONSTRAINTS_STARTS = 8;
    static final int ELEMENTS_VARIABLES = 9;
    static final int ELEMENTS_COEFFICIENTS = 10;
    static final int OBJECTIVE_QUADRATIC_ELEMENTS_VARIABLES_1 = 11;
    static final int OBJECTIVE_QUADRATIC_ELEMENTS_VARIABLES_2 = 12;
    static final int OBJECTIVE_QUADRATIC_ELEMENTS_COEFFICIENTS = 13;
    static final int QUADRATIC_ELEMENTS_CONSTRAINTS_STARTS = 14;
    static final int QUADRATIC_ELEMENTS_VARIABLES_1 = 15;
    static final int QUADRATIC_ELEMENTS_VARIABLES_2 = 16;
    static final int QUADRATIC_ELEMENTS_COEFFICIENTS = 17;
    static final int OBJECTIVE_NONLINEAR_ELEMENTS_OPERATORS = 18;
    static final int OBJECTIVE_NONLINEAR_ELEMENTS_VALUES = 19;
    static final int OBJECTIVE_NONLINEAR_ELEMENTS_VARIABLES = 20;
    static final int OBJECTIVE_NONLINEAR_ELEMENTS_NUMBER_OF_CHILDREN = 21;
    static final int NONLINEAR_ELEMENTS_CONSTRAINTS_STARTS = 22;
    static final int NONLINEAR_ELEMENTS_OPERATORS = 23;
    static final int NONLINEAR_ELEMENTS_VALUES = 24;
    static final int NONLINEAR_ELEMENTS_VARIABLES = 25;
    static final int NONLINEAR_ELEMENTS_NUMBER_OF_CHILDREN = 26;
    static final int DEFINED_VARIABLES_STARTS = 27;
    static final int DEFINED_VARIABLES_ELEMENTS_OPERATORS = 28;
    static final int DEFINED_VARIABLES_ELEMENTS_VALUES = 29;
    static final int DEFINED_VARIABLES_ELEMENTS_VARIABLES = 30;
    static final int DEFINED_VARIABLES_ELEMENTS_NUMBER_OF_CHILDREN = 31;

    /* Name lists. */
    static final int VARIABLES_NAMES = 0;
    static final int CONSTRAINTS_NAMES = 1;

    /** Size of the error message buffers. */
    static final long ERROR_SIZE = 1024;

    private static final String LIBRARY = "mathoptmodel";
    private static final Linker LINKER = Linker.nativeLinker();
    private static final SymbolLookup LOOKUP = loadLibrary();
    private static final Linker.Option CRITICAL = Linker.Option.critical(true);

    /* Life cycle. */
    static final MethodHandle MODEL_CREATE = handle("mom_model_create", false, ADDRESS);
    static final MethodHandle MODEL_FREE = handleVoid("mom_model_free", ADDRESS);
    static final MethodHandle MODEL_COPY = handle("mom_model_copy", false, ADDRESS, ADDRESS);
    static final MethodHandle FREE = handleVoid("mom_free", ADDRESS);

    /* Reading and writing. */
    static final MethodHandle READ = handle("mom_read", false, JAVA_INT,
            ADDRESS, JAVA_LONG, JAVA_INT, ADDRESS, ADDRESS, JAVA_LONG);
    static final MethodHandle READ_NL_STATISTICS = handle("mom_read_nl_statistics", false, JAVA_INT,
            ADDRESS, JAVA_LONG, ADDRESS, ADDRESS, JAVA_LONG);
    static final int NL_STATISTICS_SIZE = 14;
    static final MethodHandle READ_NL_NAMES = handle("mom_read_nl_names", false, JAVA_INT,
            ADDRESS, ADDRESS, JAVA_LONG, ADDRESS, JAVA_LONG, ADDRESS, JAVA_LONG);
    static final MethodHandle WRITE = handle("mom_write", false, JAVA_INT,
            ADDRESS, JAVA_INT, ADDRESS, ADDRESS, ADDRESS, JAVA_LONG);
    static final MethodHandle WRITE_FILE = handle("mom_write_file", false, JAVA_INT,
            ADDRESS, JAVA_INT, ADDRESS, ADDRESS, JAVA_LONG);
    static final MethodHandle FORMAT = handle("mom_format", false, JAVA_INT,
            ADDRESS, JAVA_INT, ADDRESS, ADDRESS, ADDRESS, JAVA_LONG);
    static final MethodHandle FORMAT_CONSTRAINT = handle("mom_format_constraint", false, JAVA_INT,
            ADDRESS, JAVA_INT, ADDRESS, ADDRESS, ADDRESS, JAVA_LONG);

    /* Inspection. */
    static final MethodHandle NUMBER_OF_VARIABLES = handle("mom_number_of_variables", true, JAVA_INT, ADDRESS);
    static final MethodHandle NUMBER_OF_CONSTRAINTS = handle("mom_number_of_constraints", true, JAVA_INT, ADDRESS);
    static final MethodHandle NUMBER_OF_ELEMENTS = handle("mom_number_of_elements", true, JAVA_INT, ADDRESS);
    static final MethodHandle NUMBER_OF_DEFINED_VARIABLES = handle("mom_number_of_defined_variables", true,
            JAVA_INT, ADDRESS);
    static final MethodHandle OBJECTIVE_DIRECTION = handle("mom_objective_direction", true, JAVA_INT, ADDRESS);
    static final MethodHandle ARRAY_SIZE = handle("mom_array_size", true, JAVA_LONG, ADDRESS, JAVA_INT);
    static final MethodHandle COPY_ARRAY = handle("mom_copy_array", true, JAVA_INT,
            ADDRESS, JAVA_INT, ADDRESS, ADDRESS, JAVA_LONG);
    static final MethodHandle NAMES_SIZE = handle("mom_names_size", true, JAVA_LONG,
            ADDRESS, JAVA_INT, ADDRESS);
    static final MethodHandle COPY_NAMES = handle("mom_copy_names", true, JAVA_INT,
            ADDRESS, JAVA_INT, ADDRESS, ADDRESS, ADDRESS, JAVA_LONG);
    static final MethodHandle CONSTRAINT_SENSE = handle("mom_constraint_sense", true, JAVA_INT,
            ADDRESS, JAVA_INT);
    static final MethodHandle CONSTRAINT_CLASS = handle("mom_constraint_class", true, JAVA_INT,
            ADDRESS, JAVA_INT);
    static final MethodHandle CHECK = handle("mom_check", false, JAVA_INT, ADDRESS, JAVA_INT);
    static final MethodHandle EVALUATE_OBJECTIVE = handle("mom_evaluate_objective", true, JAVA_INT,
            ADDRESS, ADDRESS, ADDRESS, ADDRESS, JAVA_LONG);
    static final MethodHandle EVALUATE_CONSTRAINT = handle("mom_evaluate_constraint", true, JAVA_INT,
            ADDRESS, ADDRESS, JAVA_INT, ADDRESS, ADDRESS, JAVA_LONG);
    static final MethodHandle CHECK_SOLUTION = handle("mom_check_solution", true, JAVA_INT,
            ADDRESS, ADDRESS, JAVA_INT, ADDRESS, ADDRESS, JAVA_LONG);

    /* Building. */
    static final MethodHandle ADD_VARIABLES = handle("mom_add_variables", true, JAVA_INT,
            ADDRESS, JAVA_LONG,
            ADDRESS, JAVA_DOUBLE,
            ADDRESS, JAVA_DOUBLE,
            ADDRESS, JAVA_BYTE,
            ADDRESS, JAVA_DOUBLE,
            ADDRESS,
            ADDRESS, ADDRESS,
            ADDRESS, ADDRESS, JAVA_LONG);
    static final MethodHandle ADD_CONSTRAINTS = handle("mom_add_constraints", true, JAVA_INT,
            ADDRESS, JAVA_LONG,
            ADDRESS, JAVA_DOUBLE,
            ADDRESS, JAVA_DOUBLE,
            ADDRESS, JAVA_LONG, ADDRESS, ADDRESS,
            ADDRESS, JAVA_LONG, ADDRESS, ADDRESS, ADDRESS,
            ADDRESS, JAVA_LONG, ADDRESS, ADDRESS, ADDRESS, ADDRESS,
            ADDRESS, ADDRESS,
            ADDRESS, ADDRESS, JAVA_LONG);
    static final MethodHandle ADD_DEFINED_VARIABLES = handle("mom_add_defined_variables", true, JAVA_INT,
            ADDRESS, JAVA_LONG, ADDRESS,
            JAVA_LONG, ADDRESS, ADDRESS, ADDRESS, ADDRESS,
            ADDRESS, ADDRESS, JAVA_LONG);
    static final MethodHandle SET_OBJECTIVE = handle("mom_set_objective", true, JAVA_INT,
            ADDRESS,
            JAVA_INT, JAVA_INT,
            JAVA_INT, ADDRESS, JAVA_DOUBLE,
            JAVA_INT, JAVA_LONG, ADDRESS, ADDRESS, ADDRESS,
            JAVA_INT, JAVA_LONG, ADDRESS, ADDRESS, ADDRESS, ADDRESS,
            ADDRESS, JAVA_LONG);

    /* Transformations. */
    static final MethodHandle DUAL = handle("mom_dual", false, JAVA_INT,
            ADDRESS, ADDRESS, ADDRESS, JAVA_LONG);
    static final MethodHandle EXPAND_DEFINED_VARIABLES = handle("mom_expand_defined_variables", false, JAVA_INT,
            ADDRESS, ADDRESS, ADDRESS, JAVA_LONG);

    private static MethodHandle handle(String name, boolean critical, MemoryLayout result, MemoryLayout... arguments) {
        return link(name, critical, FunctionDescriptor.of(result, arguments));
    }

    private static MethodHandle handleVoid(String name, MemoryLayout... arguments) {
        return link(name, false, FunctionDescriptor.ofVoid(arguments));
    }

    private static MethodHandle link(String name, boolean critical, FunctionDescriptor descriptor) {
        MemorySegment symbol = LOOKUP.find(name).orElseThrow(
                () -> new UnsatisfiedLinkError("Function " + name + " not found in the " + LIBRARY + " library."));
        return critical
                ? LINKER.downcallHandle(symbol, descriptor, CRITICAL)
                : LINKER.downcallHandle(symbol, descriptor);
    }

    /** The current platform, as named in the jar's native directory: "linux-x86-64", "macos-aarch64"... */
    static String platform() {
        String os = System.getProperty("os.name").toLowerCase(Locale.ROOT);
        String arch = System.getProperty("os.arch").toLowerCase(Locale.ROOT);
        String osName = os.contains("win") ? "windows" : (os.contains("mac") || os.contains("darwin")) ? "macos" : "linux";
        String archName = switch (arch) {
            case "amd64", "x86_64", "x86-64" -> "x86-64";
            case "aarch64", "arm64" -> "aarch64";
            default -> arch;
        };
        return osName + "-" + archName;
    }

    private static SymbolLookup loadLibrary() {
        String path = System.getProperty("mathoptmodel.library.path");
        if (path != null && !path.isEmpty())
            return SymbolLookup.libraryLookup(Path.of(path), Arena.global());

        String fileName = System.mapLibraryName(LIBRARY);
        String resource = "native/" + platform() + "/" + fileName;
        try (InputStream in = Native.class.getResourceAsStream(resource)) {
            if (in != null) {
                Path directory = Files.createTempDirectory(LIBRARY);
                Path file = directory.resolve(fileName);
                Files.copy(in, file, StandardCopyOption.REPLACE_EXISTING);
                file.toFile().deleteOnExit();
                directory.toFile().deleteOnExit();
                return SymbolLookup.libraryLookup(file, Arena.global());
            }
        } catch (IOException e) {
            throw new UncheckedIOException("Unable to extract the bundled " + LIBRARY + " library.", e);
        }

        try {
            return SymbolLookup.libraryLookup(fileName, Arena.global());
        } catch (IllegalArgumentException e) {
            UnsatisfiedLinkError error = new UnsatisfiedLinkError(
                    "The " + LIBRARY + " native library isn't bundled for platform " + platform()
                    + ", nor found on the library path; set the mathoptmodel.library.path system property"
                    + " to the path of " + fileName + ".");
            error.initCause(e);
            throw error;
        }
    }

    /** A Java array as a segment (null: MemorySegment.NULL), for critical functions. */
    static MemorySegment array(double[] array) {
        return array == null ? MemorySegment.NULL : MemorySegment.ofArray(array);
    }

    static MemorySegment array(int[] array) {
        return array == null ? MemorySegment.NULL : MemorySegment.ofArray(array);
    }

    static MemorySegment array(byte[] array) {
        return array == null ? MemorySegment.NULL : MemorySegment.ofArray(array);
    }

    static MemorySegment array(long[] array) {
        return array == null ? MemorySegment.NULL : MemorySegment.ofArray(array);
    }

    /** Throw the exception matching a C API status, if it isn't OK. */
    static void check(int status, MemorySegment errorMessage) {
        if (status == OK)
            return;
        String message = errorMessage.getString(0);
        if (status == INVALID_ARGUMENT)
            throw new IllegalArgumentException(message);
        throw new MathOptModelException(message);
    }

    /** Rethrow what a downcall handle threw. */
    static RuntimeException rethrow(Throwable throwable) {
        if (throwable instanceof RuntimeException e)
            return e;
        if (throwable instanceof Error e)
            throw e;
        return new MathOptModelException(throwable.toString(), throwable);
    }
}

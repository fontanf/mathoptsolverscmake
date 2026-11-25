package io.github.fontanf.mathoptmodel;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.lang.foreign.Arena;
import java.lang.foreign.MemorySegment;
import java.lang.ref.Cleaner;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Locale;
import java.util.Objects;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import static java.lang.foreign.ValueLayout.ADDRESS;
import static java.lang.foreign.ValueLayout.JAVA_BYTE;
import static java.lang.foreign.ValueLayout.JAVA_DOUBLE;
import static java.lang.foreign.ValueLayout.JAVA_INT;
import static java.lang.foreign.ValueLayout.JAVA_LONG;

/**
 * An optimization model: MathOptSolversCMake's {@code MathOptModel}, held by
 * the native library.
 *
 * <p>Models are read from files ({@link #read(Path)}) or built in batches of
 * arrays ({@link #addVariables(Variables)}, {@link #addConstraints(Constraints)},
 * {@link #addDefinedVariables(DefinedVariables)}, {@link #setObjective(Objective)}).
 * Each batch is validated before the model is modified: an invalid batch
 * throws {@link IllegalArgumentException} and leaves the model unchanged. The
 * getters return copies of the model's arrays; see {@code mathopt.hpp} for
 * the meaning of each of them.
 *
 * <p>A model holds native memory: close it (try-with-resources) to release it
 * early; otherwise it is released once the model is garbage collected. A model
 * must not be used from several threads at once if one of them modifies it.
 */
public final class Model implements AutoCloseable {

    private static final Cleaner CLEANER = Cleaner.create();

    private final MemorySegment handle;
    private final Cleaner.Cleanable cleanable;
    private boolean closed;

    /** An empty model. */
    public Model() {
        this(create());
    }

    private Model(MemorySegment handle) {
        this.handle = handle;
        this.cleanable = CLEANER.register(this, new Release(handle));
    }

    private static MemorySegment create() {
        try {
            MemorySegment handle = (MemorySegment) Native.MODEL_CREATE.invoke();
            if (handle.equals(MemorySegment.NULL))
                throw new OutOfMemoryError("Unable to allocate a model.");
            return handle;
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /** Frees the native model; must not reference the Model. */
    private record Release(MemorySegment handle) implements Runnable {
        @Override
        public void run() {
            try {
                Native.MODEL_FREE.invoke(handle);
            } catch (Throwable t) {
                throw Native.rethrow(t);
            }
        }
    }

    /** Release the native model. The model can't be used afterwards. */
    @Override
    public void close() {
        if (!closed) {
            closed = true;
            cleanable.clean();
        }
    }

    private MemorySegment handle() {
        if (closed)
            throw new IllegalStateException("The model is closed.");
        return handle;
    }

    /** Wrap a native model returned through a {@code mom_model**}. */
    private static Model wrap(MemorySegment pointer) {
        return new Model(pointer.get(ADDRESS, 0));
    }

    /** A copy of this model. */
    public Model copy() {
        try {
            MemorySegment copy = (MemorySegment) Native.MODEL_COPY.invoke(handle());
            if (copy.equals(MemorySegment.NULL))
                throw new OutOfMemoryError("Unable to copy the model.");
            return new Model(copy);
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /*
     * Reading and writing.
     */

    /** Read a model file, optionally .gz-compressed, in the format given by its extension. */
    public static Model read(Path path) throws IOException {
        return read(path, Format.fromPath(path).orElseThrow(
                () -> new IllegalArgumentException("Unrecognized model file format: " + path)));
    }

    /**
     * Read a model file, optionally .gz-compressed (from its extension), in the given format.
     * For an .nl file, the names come from the .col and .row files next to it, if any
     * (see {@link #readNlNames(byte[], byte[])}).
     */
    public static Model read(Path path, Format format) throws IOException {
        byte[] data;
        try (InputStream in = open(path)) {
            data = in.readAllBytes();
        }
        Model model = read(data, format);
        if (format == Format.NL) {
            byte[] col = readIfExists(nlNamesPath(path, ".col"));
            byte[] row = readIfExists(nlNamesPath(path, ".row"));
            if (col != null || row != null)
                model.readNlNames(col, row);
        }
        return model;
    }

    /**
     * Set the names of a model read from an .nl file from the contents of the .col and .row
     * files AMPL writes next to it (null: no such file). A file whose number of names doesn't
     * match the model is ignored.
     */
    public void readNlNames(byte[] col, byte[] row) {
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment colData = col == null ? MemorySegment.NULL : arena.allocateFrom(JAVA_BYTE, col);
            MemorySegment rowData = row == null ? MemorySegment.NULL : arena.allocateFrom(JAVA_BYTE, row);
            MemorySegment error = arena.allocate(Native.ERROR_SIZE);
            int status = (int) Native.READ_NL_NAMES.invoke(handle(),
                    colData, col == null ? 0L : (long) col.length,
                    rowData, row == null ? 0L : (long) row.length,
                    error, Native.ERROR_SIZE);
            Native.check(status, error);
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /** The .col / .row file of an .nl file: "model.nl.gz" gives "model.col". */
    private static Path nlNamesPath(Path path, String extension) {
        String name = path.getFileName().toString();
        String lower = name.toLowerCase(Locale.ROOT);
        if (lower.endsWith(".gz")) {
            name = name.substring(0, name.length() - 3);
            lower = lower.substring(0, lower.length() - 3);
        }
        if (lower.endsWith(".nl"))
            name = name.substring(0, name.length() - 3);
        return path.resolveSibling(name + extension);
    }

    private static byte[] readIfExists(Path path) throws IOException {
        return Files.isRegularFile(path) ? Files.readAllBytes(path) : null;
    }

    /**
     * Read the statistics of a model file, optionally .gz-compressed, in the format given by its
     * extension. For an .nl file, they come from its header only, which is fast even for very
     * large files, and works for any .nl file, including the ones {@link #read(Path)} doesn't
     * support; for the other formats, from the model read.
     */
    public static ModelStatistics readStatistics(Path path) throws IOException {
        Format format = Format.fromPath(path).orElseThrow(
                () -> new IllegalArgumentException("Unrecognized model file format: " + path));
        if (format != Format.NL) {
            try (Model model = read(path, format)) {
                return model.statistics();
            }
        }
        // The header is the first 10 lines of the file.
        ByteArrayOutputStream header = new ByteArrayOutputStream();
        try (InputStream in = open(path)) {
            int lines = 0;
            int c;
            while (lines < 10 && (c = in.read()) != -1) {
                header.write(c);
                if (c == '\n')
                    ++lines;
            }
        }
        byte[] data = header.toByteArray();
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment input = arena.allocateFrom(JAVA_BYTE, data);
            MemorySegment values = arena.allocate(JAVA_LONG, Native.STATISTICS_SIZE);
            MemorySegment error = arena.allocate(Native.ERROR_SIZE);
            int status = (int) Native.READ_NL_STATISTICS.invoke(input, (long) data.length, values, error,
                    Native.ERROR_SIZE);
            if (status != Native.OK)
                throw new MathOptModelException(error.getString(0));
            return ModelStatistics.fromValues(values.toArray(JAVA_LONG));
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /** The model's statistics: numbers of variables (by type), constraints, nonzeros... */
    public ModelStatistics statistics() {
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment values = arena.allocate(JAVA_LONG, Native.STATISTICS_SIZE);
            MemorySegment error = arena.allocate(Native.ERROR_SIZE);
            int status = (int) Native.COMPUTE_STATISTICS.invoke(handle(), values, error, Native.ERROR_SIZE);
            Native.check(status, error);
            return ModelStatistics.fromValues(values.toArray(JAVA_LONG));
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /**
     * Read a model from the contents of a file in the given format.
     *
     * @throws MathOptModelException if the contents can't be read
     */
    public static Model read(byte[] data, Format format) {
        Objects.requireNonNull(format, "format");
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment input = arena.allocateFrom(JAVA_BYTE, data);
            MemorySegment model = arena.allocate(ADDRESS);
            MemorySegment error = arena.allocate(Native.ERROR_SIZE);
            int status = (int) Native.READ.invoke(input, (long) data.length, format.code(), model, error,
                    Native.ERROR_SIZE);
            if (status != Native.OK)
                throw new MathOptModelException(error.getString(0));
            return wrap(model);
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /** Write the model to a file, .gz-compressed if its name says so, in the format given by its extension. */
    public void write(Path path) throws IOException {
        write(path, Format.fromPath(path).orElseThrow(
                () -> new IllegalArgumentException("Unrecognized model file format: " + path)));
    }

    /** Write the model to a file, .gz-compressed if its name says so, in the given format. */
    public void write(Path path, Format format) throws IOException {
        // An uncompressed file is written by the native library directly
        // (on Windows, only for ASCII paths, which it opens with the
        // narrow-string API).
        String pathString = path.toString();
        boolean windows = System.getProperty("os.name").toLowerCase(Locale.ROOT).contains("win");
        if (!isGzip(path) && (!windows || StandardCharsets.US_ASCII.newEncoder().canEncode(pathString))) {
            writeFile(pathString, format);
            return;
        }
        byte[] data = writeBytes(format);
        try (OutputStream out = isGzip(path)
                ? new GZIPOutputStream(Files.newOutputStream(path))
                : Files.newOutputStream(path)) {
            out.write(data);
        }
    }

    private void writeFile(String path, Format format) throws IOException {
        Objects.requireNonNull(format, "format");
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment nativePath = arena.allocateFrom(path, StandardCharsets.UTF_8);
            MemorySegment error = arena.allocate(Native.ERROR_SIZE);
            int status = (int) Native.WRITE_FILE.invoke(handle(), format.code(), nativePath, error,
                    Native.ERROR_SIZE);
            if (status == Native.INVALID_ARGUMENT && error.getString(0).contains("unable to open file"))
                throw new IOException(error.getString(0));
            Native.check(status, error);
        } catch (IOException | RuntimeException | Error e) {
            throw e;
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /** The model written in the given format. */
    public String writeString(Format format) {
        return new String(writeBytes(format), StandardCharsets.UTF_8);
    }

    private byte[] writeBytes(Format format) {
        Objects.requireNonNull(format, "format");
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment data = arena.allocate(ADDRESS);
            MemorySegment size = arena.allocate(JAVA_LONG);
            MemorySegment error = arena.allocate(Native.ERROR_SIZE);
            int status = (int) Native.WRITE.invoke(handle(), format.code(), data, size, error, Native.ERROR_SIZE);
            Native.check(status, error);
            return takeBytes(data, size);
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /** A description of the model (MathOptModel::format). */
    public String format(int verbosityLevel) {
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment data = arena.allocate(ADDRESS);
            MemorySegment size = arena.allocate(JAVA_LONG);
            MemorySegment error = arena.allocate(Native.ERROR_SIZE);
            int status = (int) Native.FORMAT.invoke(handle(), verbosityLevel, data, size, error, Native.ERROR_SIZE);
            Native.check(status, error);
            return new String(takeBytes(data, size), StandardCharsets.UTF_8);
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /** A description of a constraint. */
    public String formatConstraint(int constraintId) {
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment data = arena.allocate(ADDRESS);
            MemorySegment size = arena.allocate(JAVA_LONG);
            MemorySegment error = arena.allocate(Native.ERROR_SIZE);
            int status = (int) Native.FORMAT_CONSTRAINT.invoke(handle(), constraintId, data, size, error,
                    Native.ERROR_SIZE);
            Native.check(status, error);
            return new String(takeBytes(data, size), StandardCharsets.UTF_8);
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /**
     * The expression of a constraint, without its name and bounds, in infix notation. Defined
     * variables appear as {@code d<id>}, or, with {@code expandDefinedVariables}, as their
     * expressions.
     */
    public String formatConstraintExpression(int constraintId, boolean expandDefinedVariables) {
        return formatExpression(Native.EXPRESSION_CONSTRAINT, constraintId, expandDefinedVariables);
    }

    /** The expression of the objective, in infix notation. */
    public String formatObjectiveExpression(boolean expandDefinedVariables) {
        return formatExpression(Native.EXPRESSION_OBJECTIVE, 0, expandDefinedVariables);
    }

    /** The expression of a defined variable, in infix notation. */
    public String formatDefinedVariableExpression(int definedVariableId, boolean expandDefinedVariables) {
        return formatExpression(Native.EXPRESSION_DEFINED_VARIABLE, definedVariableId, expandDefinedVariables);
    }

    private String formatExpression(int kind, int id, boolean expandDefinedVariables) {
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment data = arena.allocate(ADDRESS);
            MemorySegment size = arena.allocate(JAVA_LONG);
            MemorySegment error = arena.allocate(Native.ERROR_SIZE);
            int status = (int) Native.FORMAT_EXPRESSION.invoke(handle(), kind, id, expandDefinedVariables ? 1 : 0,
                    data, size, error, Native.ERROR_SIZE);
            Native.check(status, error);
            return new String(takeBytes(data, size), StandardCharsets.UTF_8);
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /** The degree of a constraint's expression as a polynomial ({@code Double.POSITIVE_INFINITY} if it isn't one). */
    public double constraintDegree(int constraintId) {
        return degree(Native.EXPRESSION_CONSTRAINT, constraintId);
    }

    /** The degree of the objective as a polynomial ({@code Double.POSITIVE_INFINITY} if it isn't one). */
    public double objectiveDegree() {
        return degree(Native.EXPRESSION_OBJECTIVE, 0);
    }

    private double degree(int kind, int id) {
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment degree = arena.allocate(JAVA_DOUBLE);
            MemorySegment error = arena.allocate(Native.ERROR_SIZE);
            int status = (int) Native.DEGREE.invoke(handle(), kind, id, degree, error, Native.ERROR_SIZE);
            Native.check(status, error);
            return degree.get(JAVA_DOUBLE, 0);
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /** Copy and free a buffer returned by the library through {@code char** data, int64_t* size}. */
    private static byte[] takeBytes(MemorySegment dataPointer, MemorySegment sizePointer) throws Throwable {
        long size = sizePointer.get(JAVA_LONG, 0);
        MemorySegment data = dataPointer.get(ADDRESS, 0).reinterpret(size);
        try {
            if (size > Integer.MAX_VALUE - 8)
                throw new MathOptModelException("The output is too large for a Java array.");
            return data.toArray(JAVA_BYTE);
        } finally {
            Native.FREE.invoke(data);
        }
    }

    private static boolean isGzip(Path path) {
        return path.getFileName().toString().toLowerCase(Locale.ROOT).endsWith(".gz");
    }

    private static InputStream open(Path path) throws IOException {
        InputStream in = Files.newInputStream(path);
        return isGzip(path) ? new GZIPInputStream(in) : in;
    }

    /*
     * Sizes and objective direction.
     */

    public int numberOfVariables() {
        return intCall(Native.NUMBER_OF_VARIABLES);
    }

    public int numberOfConstraints() {
        return intCall(Native.NUMBER_OF_CONSTRAINTS);
    }

    public int numberOfElements() {
        return intCall(Native.NUMBER_OF_ELEMENTS);
    }

    public int numberOfDefinedVariables() {
        return intCall(Native.NUMBER_OF_DEFINED_VARIABLES);
    }

    public ObjectiveDirection objectiveDirection() {
        return ObjectiveDirection.fromCode(intCall(Native.OBJECTIVE_DIRECTION));
    }

    public void setObjectiveDirection(ObjectiveDirection direction) {
        setObjective(new Objective().direction(direction));
    }

    private int intCall(java.lang.invoke.MethodHandle function) {
        try {
            return (int) function.invoke(handle());
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /*
     * Arrays (copies).
     */

    public double[] variablesLowerBounds() {
        return doubles(Native.VARIABLES_LOWER_BOUNDS);
    }

    public double[] variablesUpperBounds() {
        return doubles(Native.VARIABLES_UPPER_BOUNDS);
    }

    public VariableType[] variablesTypes() {
        byte[] codes = variablesTypesCodes();
        VariableType[] types = new VariableType[codes.length];
        for (int i = 0; i < codes.length; ++i)
            types[i] = VariableType.fromCode(codes[i]);
        return types;
    }

    /** Variable types as codes ({@link VariableType#code()}). */
    public byte[] variablesTypesCodes() {
        return bytes(Native.VARIABLES_TYPES);
    }

    /** Initial values (empty if the model has none). */
    public double[] variablesInitialValues() {
        return doubles(Native.VARIABLES_INITIAL_VALUES);
    }

    /** Variable names (empty if the model has none). */
    public List<String> variablesNames() {
        return names(Native.VARIABLES_NAMES);
    }

    public double[] constraintsLowerBounds() {
        return doubles(Native.CONSTRAINTS_LOWER_BOUNDS);
    }

    public double[] constraintsUpperBounds() {
        return doubles(Native.CONSTRAINTS_UPPER_BOUNDS);
    }

    /** Constraint names (empty if the model has none). */
    public List<String> constraintsNames() {
        return names(Native.CONSTRAINTS_NAMES);
    }

    public ConstraintClass[] constraintsClasses() {
        byte[] codes = bytes(Native.CONSTRAINTS_CLASSES);
        ConstraintClass[] classes = new ConstraintClass[codes.length];
        for (int i = 0; i < codes.length; ++i)
            classes[i] = ConstraintClass.fromCode(codes[i]);
        return classes;
    }

    public double[] objectiveCoefficients() {
        return doubles(Native.OBJECTIVE_COEFFICIENTS);
    }

    public int[] constraintsStarts() {
        return ints(Native.CONSTRAINTS_STARTS);
    }

    public int[] elementsVariables() {
        return ints(Native.ELEMENTS_VARIABLES);
    }

    public double[] elementsCoefficients() {
        return doubles(Native.ELEMENTS_COEFFICIENTS);
    }

    public int[] objectiveQuadraticElementsVariables1() {
        return ints(Native.OBJECTIVE_QUADRATIC_ELEMENTS_VARIABLES_1);
    }

    public int[] objectiveQuadraticElementsVariables2() {
        return ints(Native.OBJECTIVE_QUADRATIC_ELEMENTS_VARIABLES_2);
    }

    public double[] objectiveQuadraticElementsCoefficients() {
        return doubles(Native.OBJECTIVE_QUADRATIC_ELEMENTS_COEFFICIENTS);
    }

    public int[] quadraticElementsConstraintsStarts() {
        return ints(Native.QUADRATIC_ELEMENTS_CONSTRAINTS_STARTS);
    }

    public int[] quadraticElementsVariables1() {
        return ints(Native.QUADRATIC_ELEMENTS_VARIABLES_1);
    }

    public int[] quadraticElementsVariables2() {
        return ints(Native.QUADRATIC_ELEMENTS_VARIABLES_2);
    }

    public double[] quadraticElementsCoefficients() {
        return doubles(Native.QUADRATIC_ELEMENTS_COEFFICIENTS);
    }

    public byte[] objectiveNonlinearElementsOperators() {
        return bytes(Native.OBJECTIVE_NONLINEAR_ELEMENTS_OPERATORS);
    }

    public double[] objectiveNonlinearElementsValues() {
        return doubles(Native.OBJECTIVE_NONLINEAR_ELEMENTS_VALUES);
    }

    public int[] objectiveNonlinearElementsVariables() {
        return ints(Native.OBJECTIVE_NONLINEAR_ELEMENTS_VARIABLES);
    }

    public int[] objectiveNonlinearElementsNumberOfChildren() {
        return ints(Native.OBJECTIVE_NONLINEAR_ELEMENTS_NUMBER_OF_CHILDREN);
    }

    /** Offset of each constraint's expression tree (empty if the model has no constraint trees). */
    public int[] nonlinearElementsConstraintsStarts() {
        return ints(Native.NONLINEAR_ELEMENTS_CONSTRAINTS_STARTS);
    }

    public byte[] nonlinearElementsOperators() {
        return bytes(Native.NONLINEAR_ELEMENTS_OPERATORS);
    }

    public double[] nonlinearElementsValues() {
        return doubles(Native.NONLINEAR_ELEMENTS_VALUES);
    }

    public int[] nonlinearElementsVariables() {
        return ints(Native.NONLINEAR_ELEMENTS_VARIABLES);
    }

    public int[] nonlinearElementsNumberOfChildren() {
        return ints(Native.NONLINEAR_ELEMENTS_NUMBER_OF_CHILDREN);
    }

    public int[] definedVariablesStarts() {
        return ints(Native.DEFINED_VARIABLES_STARTS);
    }

    public byte[] definedVariablesElementsOperators() {
        return bytes(Native.DEFINED_VARIABLES_ELEMENTS_OPERATORS);
    }

    public double[] definedVariablesElementsValues() {
        return doubles(Native.DEFINED_VARIABLES_ELEMENTS_VALUES);
    }

    public int[] definedVariablesElementsVariables() {
        return ints(Native.DEFINED_VARIABLES_ELEMENTS_VARIABLES);
    }

    public int[] definedVariablesElementsNumberOfChildren() {
        return ints(Native.DEFINED_VARIABLES_ELEMENTS_NUMBER_OF_CHILDREN);
    }

    private int arraySize(int array) throws Throwable {
        long size = (long) Native.ARRAY_SIZE.invoke(handle(), array);
        if (size < 0 || size > Integer.MAX_VALUE - 8)
            throw new MathOptModelException("Invalid array size " + size + ".");
        return (int) size;
    }

    private void copyArray(int array, MemorySegment out) throws Throwable {
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment error = arena.allocate(Native.ERROR_SIZE);
            int status = (int) Native.COPY_ARRAY.invoke(handle(), array, out, error, Native.ERROR_SIZE);
            Native.check(status, error);
        }
    }

    private double[] doubles(int array) {
        try {
            double[] res = new double[arraySize(array)];
            copyArray(array, MemorySegment.ofArray(res));
            return res;
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    private int[] ints(int array) {
        try {
            int[] res = new int[arraySize(array)];
            copyArray(array, MemorySegment.ofArray(res));
            return res;
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    private byte[] bytes(int array) {
        try {
            byte[] res = new byte[arraySize(array)];
            copyArray(array, MemorySegment.ofArray(res));
            return res;
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    private List<String> names(int list) {
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment totalSize = arena.allocate(JAVA_LONG);
            long count = (long) Native.NAMES_SIZE.invoke(handle(), list, totalSize);
            long total = totalSize.get(JAVA_LONG, 0);
            if (count == 0)
                return Collections.emptyList();
            if (count > Integer.MAX_VALUE - 8 || total > Integer.MAX_VALUE - 8)
                throw new MathOptModelException("Too many names for Java arrays.");
            byte[] data = new byte[(int) total];
            long[] offsets = new long[(int) count + 1];
            MemorySegment error = arena.allocate(Native.ERROR_SIZE);
            int status = (int) Native.COPY_NAMES.invoke(handle(), list,
                    MemorySegment.ofArray(data), MemorySegment.ofArray(offsets), error, Native.ERROR_SIZE);
            Native.check(status, error);
            List<String> names = new ArrayList<>((int) count);
            for (int i = 0; i < count; ++i) {
                names.add(new String(data, (int) offsets[i], (int) (offsets[i + 1] - offsets[i]),
                        StandardCharsets.UTF_8));
            }
            return Collections.unmodifiableList(names);
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /*
     * Model kind.
     */

    public boolean hasNonContinuousVariables() {
        for (byte code : variablesTypesCodes()) {
            if (code != VariableType.CONTINUOUS.code())
                return true;
        }
        return false;
    }

    public boolean hasQuadratic() {
        return sizeOf(Native.OBJECTIVE_QUADRATIC_ELEMENTS_VARIABLES_1) > 0
                || sizeOf(Native.QUADRATIC_ELEMENTS_VARIABLES_1) > 0;
    }

    public boolean hasNonlinear() {
        return sizeOf(Native.OBJECTIVE_NONLINEAR_ELEMENTS_OPERATORS) > 0
                || sizeOf(Native.NONLINEAR_ELEMENTS_OPERATORS) > 0;
    }

    public boolean hasDefinedVariables() {
        return numberOfDefinedVariables() > 0;
    }

    /** True for a mixed-integer linear program (no quadratic or nonlinear terms). */
    public boolean isMilp() {
        return !hasQuadratic() && !hasNonlinear();
    }

    /** True for a linear program (a MILP with only continuous variables). */
    public boolean isLp() {
        return isMilp() && !hasNonContinuousVariables();
    }

    private int sizeOf(int array) {
        try {
            return arraySize(array);
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /*
     * Constraints, evaluations and checks.
     */

    public ConstraintSense constraintSense(int constraintId) {
        checkConstraint(constraintId);
        try {
            return ConstraintSense.fromCode((int) Native.CONSTRAINT_SENSE.invoke(handle(), constraintId));
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    public ConstraintClass constraintClass(int constraintId) {
        checkConstraint(constraintId);
        try {
            return ConstraintClass.fromCode((int) Native.CONSTRAINT_CLASS.invoke(handle(), constraintId));
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    private void checkConstraint(int constraintId) {
        int number = numberOfConstraints();
        if (constraintId < 0 || constraintId >= number)
            throw new IndexOutOfBoundsException("Constraint " + constraintId + " out of range [0, " + number + ").");
    }

    /** True if the model's arrays are consistent. */
    public boolean check() {
        return check(0);
    }

    /** True if the model's arrays are consistent; the inconsistencies are printed on the standard output if {@code verbosityLevel > 0}. */
    public boolean check(int verbosityLevel) {
        try {
            return (int) Native.CHECK.invoke(handle(), verbosityLevel) != 0;
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    public double evaluateObjective(double[] solution) {
        checkSolutionLength(solution);
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment value = arena.allocate(JAVA_DOUBLE);
            MemorySegment error = arena.allocate(Native.ERROR_SIZE);
            int status = (int) Native.EVALUATE_OBJECTIVE.invoke(handle(), MemorySegment.ofArray(solution), value,
                    error, Native.ERROR_SIZE);
            Native.check(status, error);
            return value.get(JAVA_DOUBLE, 0);
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /** The value of a constraint's body. */
    public double evaluateConstraint(double[] solution, int constraintId) {
        checkSolutionLength(solution);
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment value = arena.allocate(JAVA_DOUBLE);
            MemorySegment error = arena.allocate(Native.ERROR_SIZE);
            int status = (int) Native.EVALUATE_CONSTRAINT.invoke(handle(), MemorySegment.ofArray(solution),
                    constraintId, value, error, Native.ERROR_SIZE);
            Native.check(status, error);
            return value.get(JAVA_DOUBLE, 0);
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /** True if a solution is feasible. */
    public boolean checkSolution(double[] solution) {
        checkSolutionLength(solution);
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment feasible = arena.allocate(JAVA_INT);
            MemorySegment error = arena.allocate(Native.ERROR_SIZE);
            int status = (int) Native.CHECK_SOLUTION.invoke(handle(), MemorySegment.ofArray(solution), 0,
                    feasible, error, Native.ERROR_SIZE);
            Native.check(status, error);
            return feasible.get(JAVA_INT, 0) != 0;
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    private void checkSolutionLength(double[] solution) {
        Batches.checkLength(solution.length, numberOfVariables(), "solution", "the model");
    }

    /*
     * Building.
     */

    /** Append a batch of variables; return the index of the first one. */
    public int addVariables(Variables batch) {
        Batches.EncodedNames names = Batches.EncodedNames.of(batch.names);
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment first = arena.allocate(JAVA_INT);
            MemorySegment error = arena.allocate(Native.ERROR_SIZE);
            int status = (int) Native.ADD_VARIABLES.invoke(handle(), (long) batch.number,
                    Native.array(batch.lowerBounds), batch.lowerBound,
                    Native.array(batch.upperBounds), batch.upperBound,
                    Native.array(batch.types), batch.type.code(),
                    Native.array(batch.objectiveCoefficients), batch.objectiveCoefficient,
                    Native.array(batch.initialValues),
                    Native.array(names.data()), Native.array(names.offsets()),
                    first, error, Native.ERROR_SIZE);
            Native.check(status, error);
            return first.get(JAVA_INT, 0);
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /** Append a batch of constraints; return the index of the first one. */
    public int addConstraints(Constraints batch) {
        Batches.EncodedNames names = Batches.EncodedNames.of(batch.names);
        ExpressionTrees trees = batch.nonlinear;
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment first = arena.allocate(JAVA_INT);
            MemorySegment error = arena.allocate(Native.ERROR_SIZE);
            int status = (int) Native.ADD_CONSTRAINTS.invoke(handle(), (long) batch.number,
                    Native.array(batch.lowerBounds), batch.lowerBound,
                    Native.array(batch.upperBounds), batch.upperBound,
                    Native.array(batch.starts),
                    (long) (batch.variables == null ? 0 : batch.variables.length),
                    Native.array(batch.variables),
                    Native.array(batch.coefficients),
                    Native.array(batch.quadraticStarts),
                    (long) (batch.quadraticVariables1 == null ? 0 : batch.quadraticVariables1.length),
                    Native.array(batch.quadraticVariables1),
                    Native.array(batch.quadraticVariables2),
                    Native.array(batch.quadraticCoefficients),
                    Native.array(batch.nonlinearStarts),
                    (long) (trees == null ? 0 : trees.numberOfNodes()),
                    Native.array(trees == null ? null : trees.operators),
                    Native.array(trees == null ? null : trees.values),
                    Native.array(trees == null ? null : trees.variables),
                    Native.array(trees == null ? null : trees.numberOfChildren),
                    Native.array(names.data()), Native.array(names.offsets()),
                    first, error, Native.ERROR_SIZE);
            Native.check(status, error);
            return first.get(JAVA_INT, 0);
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /** Append a batch of defined variables; return the index of the first one. */
    public int addDefinedVariables(DefinedVariables batch) {
        ExpressionTrees trees = batch.trees;
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment first = arena.allocate(JAVA_INT);
            MemorySegment error = arena.allocate(Native.ERROR_SIZE);
            int status = (int) Native.ADD_DEFINED_VARIABLES.invoke(handle(), (long) batch.starts.length,
                    Native.array(batch.starts),
                    (long) trees.numberOfNodes(),
                    Native.array(trees.operators),
                    Native.array(trees.values),
                    Native.array(trees.variables),
                    Native.array(trees.numberOfChildren),
                    first, error, Native.ERROR_SIZE);
            Native.check(status, error);
            return first.get(JAVA_INT, 0);
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /** Update the objective: each part set in {@code objective} replaces the current one. */
    public void setObjective(Objective objective) {
        if (objective.coefficients != null)
            Batches.checkLength(objective.coefficients.length, numberOfVariables(), "coefficients", "the model");
        ExpressionTrees tree = objective.nonlinear;
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment error = arena.allocate(Native.ERROR_SIZE);
            int status = (int) Native.SET_OBJECTIVE.invoke(handle(),
                    objective.direction != null ? 1 : 0,
                    objective.direction != null ? objective.direction.ordinal() : 0,
                    objective.setCoefficients ? 1 : 0,
                    Native.array(objective.coefficients), objective.coefficient,
                    objective.setQuadratic ? 1 : 0,
                    (long) (objective.quadraticVariables1 == null ? 0 : objective.quadraticVariables1.length),
                    Native.array(objective.quadraticVariables1),
                    Native.array(objective.quadraticVariables2),
                    Native.array(objective.quadraticCoefficients),
                    tree != null ? 1 : 0,
                    (long) (tree == null ? 0 : tree.numberOfNodes()),
                    Native.array(tree == null ? null : tree.operators),
                    Native.array(tree == null ? null : tree.values),
                    Native.array(tree == null ? null : tree.variables),
                    Native.array(tree == null ? null : tree.numberOfChildren),
                    error, Native.ERROR_SIZE);
            Native.check(status, error);
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    /*
     * Transformations.
     */

    /** The dual of a linear program. */
    public Model dual() {
        return transform(Native.DUAL);
    }

    /** An equivalent model whose expression trees reference no defined variable. */
    public Model expandDefinedVariables() {
        return transform(Native.EXPAND_DEFINED_VARIABLES);
    }

    private Model transform(java.lang.invoke.MethodHandle function) {
        try (Arena arena = Arena.ofConfined()) {
            MemorySegment output = arena.allocate(ADDRESS);
            MemorySegment error = arena.allocate(Native.ERROR_SIZE);
            int status = (int) function.invoke(handle(), output, error, Native.ERROR_SIZE);
            Native.check(status, error);
            return wrap(output);
        } catch (Throwable t) {
            throw Native.rethrow(t);
        }
    }

    @Override
    public String toString() {
        if (closed)
            return "Model[closed]";
        return "Model[" + numberOfVariables() + " variables, " + numberOfConstraints() + " constraints]";
    }
}

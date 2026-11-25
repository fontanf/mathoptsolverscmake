package bench;

import com.google.ortools.Loader;
import com.google.ortools.linearsolver.MPConstraint;
import com.google.ortools.linearsolver.MPConstraintProto;
import com.google.ortools.linearsolver.MPModelProto;
import com.google.ortools.linearsolver.MPObjective;
import com.google.ortools.linearsolver.MPSolver;
import com.google.ortools.linearsolver.MPVariable;
import com.google.ortools.linearsolver.MPVariableProto;
import com.google.ortools.modelbuilder.LinearExpr;
import com.google.ortools.modelbuilder.ModelBuilder;
import com.google.ortools.modelbuilder.ModelBuilderHelper;
import com.google.ortools.modelbuilder.Variable;
import io.github.fontanf.mathoptmodel.Constraints;
import io.github.fontanf.mathoptmodel.Model;
import io.github.fontanf.mathoptmodel.ObjectiveDirection;
import io.github.fontanf.mathoptmodel.VariableType;
import io.github.fontanf.mathoptmodel.Variables;

import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Time building (and writing as LP) the benchmark model with a Java tool.
 *
 * <pre>java bench.Bench TOOL N REPETITIONS OUTPUT_DIRECTORY</pre>
 *
 * <p>The model data is generated exactly as in data.py. Prints one JSON line
 * with the best build and write times over the repetitions; the first
 * repetitions also warm the JIT up.
 */
public final class Bench {

    static final int K = 5;
    static final double UPPER_BOUND = 10.0;

    /** The benchmark data (see data.py). */
    record Data(int n, double[] obj, int[] cols, double[] coefs, double[] rhs) {

        static long hash(long stream, long index) {
            long z = (stream << 40) + index + 0x9E3779B97F4A7C15L;
            z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
            z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
            return z ^ (z >>> 31);
        }

        static int mod(long value, long modulus) {
            return (int) Long.remainderUnsigned(value, modulus);
        }

        static Data generate(int n) {
            int m = n;
            int[] cols = new int[m * K];
            double[] coefs = new double[m * K];
            double[] rhs = new double[m];
            double[] obj = new double[n];
            for (int i = 0; i < m; ++i) {
                long base = mod(hash(1, i), n);
                for (int k = 0; k < K; ++k) {
                    cols[i * K + k] = (int) ((base + (long) k * (n / K)) % n);
                    coefs[i * K + k] = 1 + mod(hash(2, (long) i * K + k), 9);
                }
                rhs[i] = 10 + mod(hash(3, i), 90);
            }
            for (int j = 0; j < n; ++j)
                obj[j] = 1 + mod(hash(4, j), 9);
            return new Data(n, obj, cols, coefs, rhs);
        }

        int numberOfIntegers() {
            return n / 2;
        }
    }

    /** A tool: builds the model, then writes it as LP. */
    interface Tool {
        Object build(Data data);

        void write(Object model, Path path) throws IOException;

        default void release(Object model) {
        }
    }

    /** mathoptmodel, all the constraints in one batch. */
    static final class MathOptModelBatch implements Tool {
        @Override
        public Object build(Data data) {
            int n = data.n();
            byte[] types = new byte[n];
            for (int j = 0; j < data.numberOfIntegers(); ++j)
                types[j] = VariableType.INTEGER.code();
            int[] starts = new int[n];
            for (int i = 0; i < n; ++i)
                starts[i] = i * K;
            Model model = new Model();
            model.addVariables(new Variables(n).lowerBounds(0).upperBounds(UPPER_BOUND).types(types)
                    .objectiveCoefficients(data.obj()));
            model.setObjectiveDirection(ObjectiveDirection.MAXIMIZE);
            model.addConstraints(new Constraints(n).upperBounds(data.rhs())
                    .linear(starts, data.cols(), data.coefs()));
            return model;
        }

        @Override
        public void write(Object model, Path path) throws IOException {
            ((Model) model).write(path);
        }

        @Override
        public void release(Object model) {
            ((Model) model).close();
        }
    }

    /** mathoptmodel, one call per constraint. */
    static final class MathOptModelRow implements Tool {
        @Override
        public Object build(Data data) {
            int n = data.n();
            int integers = data.numberOfIntegers();
            Model model = new Model();
            model.addVariables(new Variables(integers).lowerBounds(0).upperBounds(UPPER_BOUND)
                    .types(VariableType.INTEGER));
            model.addVariables(new Variables(n - integers).lowerBounds(0).upperBounds(UPPER_BOUND));
            model.setObjective(new io.github.fontanf.mathoptmodel.Objective()
                    .direction(ObjectiveDirection.MAXIMIZE).coefficients(data.obj()));
            int[] starts = {0};
            int[] rowCols = new int[K];
            double[] rowCoefs = new double[K];
            for (int i = 0; i < n; ++i) {
                System.arraycopy(data.cols(), i * K, rowCols, 0, K);
                System.arraycopy(data.coefs(), i * K, rowCoefs, 0, K);
                model.addConstraints(new Constraints(1).upperBounds(data.rhs()[i])
                        .linear(starts, rowCols, rowCoefs));
            }
            return model;
        }

        @Override
        public void write(Object model, Path path) throws IOException {
            ((Model) model).write(path);
        }

        @Override
        public void release(Object model) {
            ((Model) model).close();
        }
    }

    /** OR-Tools' solver-independent ModelBuilder. */
    static final class OrToolsModelBuilder implements Tool {
        @Override
        public Object build(Data data) {
            int n = data.n();
            int integers = data.numberOfIntegers();
            ModelBuilder model = new ModelBuilder();
            Variable[] x = new Variable[n];
            for (int j = 0; j < n; ++j)
                x[j] = model.newVar(0, UPPER_BOUND, j < integers, "");
            Variable[] rowVariables = new Variable[K];
            double[] rowCoefs = new double[K];
            for (int i = 0; i < n; ++i) {
                for (int k = 0; k < K; ++k) {
                    rowVariables[k] = x[data.cols()[i * K + k]];
                    rowCoefs[k] = data.coefs()[i * K + k];
                }
                model.addLessOrEqual(LinearExpr.weightedSum(rowVariables, rowCoefs), data.rhs()[i]);
            }
            model.maximize(LinearExpr.weightedSum(x, data.obj()));
            return model;
        }

        @Override
        public void write(Object model, Path path) throws IOException {
            Files.writeString(path, ((ModelBuilder) model).exportToLpString(false));
        }
    }

    /**
     * OR-Tools' ModelBuilder, built through its low-level helper (ModelBuilder.getHelper()): one
     * native call per variable attribute and per term, without Java expression objects.
     */
    static final class OrToolsModelBuilderHelper implements Tool {
        @Override
        public Object build(Data data) {
            int n = data.n();
            int integers = data.numberOfIntegers();
            ModelBuilder model = new ModelBuilder();
            ModelBuilderHelper helper = model.getHelper();
            for (int j = 0; j < n; ++j) {
                int variable = helper.addVar();
                helper.setVarLowerBound(variable, 0);
                helper.setVarUpperBound(variable, UPPER_BOUND);
                helper.setVarIntegrality(variable, j < integers);
                helper.setVarObjectiveCoefficient(variable, data.obj()[j]);
            }
            helper.setMaximize(true);
            for (int i = 0; i < n; ++i) {
                int constraint = helper.addLinearConstraint();
                helper.setConstraintLowerBound(constraint, Double.NEGATIVE_INFINITY);
                helper.setConstraintUpperBound(constraint, data.rhs()[i]);
                for (int k = 0; k < K; ++k)
                    helper.addConstraintTerm(constraint, data.cols()[i * K + k], data.coefs()[i * K + k]);
            }
            return model;
        }

        @Override
        public void write(Object model, Path path) throws IOException {
            Files.writeString(path, ((ModelBuilder) model).exportToLpString(false));
        }
    }

    /** OR-Tools' MPSolver (linear_solver), with the SCIP backend for the integer variables. */
    static class OrToolsMpSolver implements Tool {
        @Override
        public Object build(Data data) {
            int n = data.n();
            int integers = data.numberOfIntegers();
            MPSolver solver = MPSolver.createSolver("SCIP");
            MPVariable[] x = new MPVariable[n];
            for (int j = 0; j < n; ++j)
                x[j] = solver.makeVar(0, UPPER_BOUND, j < integers, "");
            for (int i = 0; i < n; ++i) {
                MPConstraint constraint = solver.makeConstraint(Double.NEGATIVE_INFINITY, data.rhs()[i], "");
                for (int k = 0; k < K; ++k)
                    constraint.setCoefficient(x[data.cols()[i * K + k]], data.coefs()[i * K + k]);
            }
            MPObjective objective = solver.objective();
            for (int j = 0; j < n; ++j)
                objective.setCoefficient(x[j], data.obj()[j]);
            objective.setMaximization();
            return solver;
        }

        @Override
        public void write(Object model, Path path) throws IOException {
            Files.writeString(path, ((MPSolver) model).exportModelAsLpFormat());
        }

        @Override
        public void release(Object model) {
            ((MPSolver) model).delete();
        }
    }

    /** OR-Tools' MPSolver, loaded from an MPModelProto built beforehand (in the timing). */
    static final class OrToolsMpSolverProto extends OrToolsMpSolver {
        @Override
        public Object build(Data data) {
            int n = data.n();
            int integers = data.numberOfIntegers();
            MPModelProto.Builder proto = MPModelProto.newBuilder().setMaximize(true);
            for (int j = 0; j < n; ++j) {
                proto.addVariable(MPVariableProto.newBuilder().setLowerBound(0).setUpperBound(UPPER_BOUND)
                        .setIsInteger(j < integers).setObjectiveCoefficient(data.obj()[j]));
            }
            for (int i = 0; i < n; ++i) {
                MPConstraintProto.Builder constraint = MPConstraintProto.newBuilder()
                        .setLowerBound(Double.NEGATIVE_INFINITY).setUpperBound(data.rhs()[i]);
                for (int k = 0; k < K; ++k) {
                    constraint.addVarIndex(data.cols()[i * K + k]);
                    constraint.addCoefficient(data.coefs()[i * K + k]);
                }
                proto.addConstraint(constraint);
            }
            MPSolver solver = MPSolver.createSolver("SCIP");
            String error = solver.loadModelFromProto(proto.build());
            if (!error.isEmpty())
                throw new IllegalStateException(error);
            return solver;
        }
    }

    static Tool tool(String name) {
        return switch (name) {
            case "mathoptmodel-batch" -> new MathOptModelBatch();
            case "mathoptmodel-row" -> new MathOptModelRow();
            case "ortools-modelbuilder" -> new OrToolsModelBuilder();
            case "ortools-modelbuilder-helper" -> new OrToolsModelBuilderHelper();
            case "ortools-mpsolver" -> new OrToolsMpSolver();
            case "ortools-mpsolver-proto" -> new OrToolsMpSolverProto();
            default -> throw new IllegalArgumentException("Unknown tool " + name);
        };
    }

    public static void main(String[] args) throws IOException {
        String name = args[0];
        int n = Integer.parseInt(args[1]);
        int repetitions = Integer.parseInt(args[2]);
        Path output = Path.of(args[3]);
        if (name.startsWith("ortools"))
            Loader.loadNativeLibraries();
        Tool tool = tool(name);
        Data data = Data.generate(n);
        Path path = output.resolve(name + "-" + n + ".lp");
        double bestBuild = Double.POSITIVE_INFINITY;
        double bestWrite = Double.POSITIVE_INFINITY;
        for (int repetition = 0; repetition < repetitions; ++repetition) {
            System.gc();
            long start = System.nanoTime();
            Object model = tool.build(data);
            bestBuild = Math.min(bestBuild, (System.nanoTime() - start) / 1e9);
            start = System.nanoTime();
            tool.write(model, path);
            bestWrite = Math.min(bestWrite, (System.nanoTime() - start) / 1e9);
            tool.release(model);
        }
        // Peak resident memory of the process (Linux), heap and native included.
        double peakMemory = Double.NaN;
        for (String line : Files.readAllLines(Path.of("/proc/self/status"))) {
            if (line.startsWith("VmHWM:"))
                peakMemory = Long.parseLong(line.replaceAll("[^0-9]", "")) / 1024.0;
        }
        System.out.printf(java.util.Locale.ROOT,
                "{\"language\": \"java\", \"tool\": \"%s\", \"n\": %d, \"build\": %.6f, \"write\": %.6f, "
                        + "\"peak_memory_mb\": %.1f, \"lp\": \"%s\", \"jvm\": \"%s\"}%n",
                name, n, bestBuild, bestWrite, peakMemory, path, ManagementFactory.getRuntimeMXBean().getVmVersion());
    }
}

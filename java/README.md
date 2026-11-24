# mathoptmodel (Java)

Java bindings for the solver-independent part of
[MathOptSolversCMake](https://github.com/fontanf/mathoptsolverscmake): its
`MathOptModel` optimization-model structure and its MPS, LP, JSON and AMPL
`.nl` readers and writers. No solver is bundled.

They call the library's C API (`include/mathoptsolverscmake/mathoptmodel_c.h`,
built as the `mathoptmodel` native library) through the Foreign Function &
Memory API, and require **Java 22 or newer**. The jar bundles the native
library for Linux (x86-64, aarch64; glibc 2.28+), macOS (x86-64, arm64; 11+)
and Windows (x86-64).

Native access must be enabled, or Java prints a warning (and future versions
will refuse it):

```shell
java --enable-native-access=ALL-UNNAMED ...
```

## Reading and inspecting a model

```java
import io.github.fontanf.mathoptmodel.*;

try (Model model = Model.read(Path.of("p0033.mps.gz"))) {  // format from the extension; .gz handled
    System.out.println(model.numberOfVariables() + " " + model.numberOfConstraints());
    int[] starts = model.constraintsStarts();          // getters return copies
    int[] variables = model.elementsVariables();
    double[] coefficients = model.elementsCoefficients();
    ConstraintClass[] classes = model.constraintsClasses();
    model.write(Path.of("p0033.lp"));
}
```

The getters mirror `MathOptModel`'s fields in
[`mathopt.hpp`](../include/mathoptsolverscmake/mathopt.hpp), which documents
each of them, including the expression-tree layout of nonlinear objectives
and constraints.

## Building a model

Models are built in batches of arrays laid out like `MathOptModel`'s own
fields. Building a model then costs a few array copies, rather than one
native call per variable, constraint or term:

```java
try (Model model = new Model()) {
    // max x + 2 y + 3 z, x in [0, 10] continuous, y binary, z in [0, 5] integer.
    model.addVariables(new Variables(3)
            .lowerBounds(0.0)                          // a single value applies to the whole batch
            .upperBounds(new double[] {10, 1, 5})
            .types(new VariableType[] {VariableType.CONTINUOUS, VariableType.BINARY, VariableType.INTEGER})
            .names("x", "y", "z"));
    model.setObjective(new Objective()
            .direction(ObjectiveDirection.MAXIMIZE)
            .coefficients(new double[] {1, 2, 3}));

    // x + y + z <= 10 and x - y >= -3, in CSR form: the terms of all the rows
    // concatenated, and the offset of each row's first term.
    model.addConstraints(new Constraints(2)
            .lowerBounds(new double[] {Double.NEGATIVE_INFINITY, -3})
            .upperBounds(new double[] {10, Double.POSITIVE_INFINITY})
            .linear(new int[] {0, 3}, new int[] {0, 1, 2, 0, 1}, new double[] {1, 1, 1, 1, -1}));

    // exp(x) + y^2 <= 4, as a pre-order expression tree.
    model.addConstraints(new Constraints(1)
            .upperBounds(4)
            .nonlinear(new int[] {0}, new ExpressionTrees("+evpvk",
                    new int[] {2, 1, 0, 2, 0, 0},           // number of children
                    new int[] {-1, -1, 0, -1, 1, -1},       // variables
                    new double[] {0, 0, 0, 0, 0, 2})));     // constants

    model.write(Path.of("model.nl"));
}
```

- `addVariables` and `addConstraints` return the index of the first entry of
  the batch.
- Quadratic constraint terms (`Constraints.quadratic`) are laid out like the
  linear ones.
- Common subexpressions go through `addDefinedVariables`.
- `setObjective` replaces only the parts it is given: the direction, the
  linear coefficients, the quadratic terms and the expression tree.
- Each batch is validated before the model is modified. An invalid batch
  throws `IllegalArgumentException` and leaves the model unchanged.
- The batch classes reference the arrays they are given, without copying
  them.

A `Model` holds native memory: close it (try-with-resources) to release it
early; otherwise it is released once the model is garbage collected. A model
must not be used from several threads at once if one of them modifies it.

## Building from source

Build the native library, then run the tests against it:

```shell
cmake -S . -B build_java -DCMAKE_BUILD_TYPE=Release -DMATHOPTSOLVERSCMAKE_BUILD_C_API=ON \
    -DMATHOPTSOLVERSCMAKE_BUILD_TESTS=OFF -DMATHOPTSOLVERSCMAKE_BUILD_APP=OFF
cmake --build build_java
cd java
mvn test -Dmathoptmodel.library.path=$PWD/../build_java/src/libmathoptmodel.so
```

To bundle native libraries in the jar, put them in
`java/natives/<platform>/` (`linux-x86-64`, `linux-aarch64`, `macos-x86-64`,
`macos-aarch64`, `windows-x86-64`) before `mvn package`; this is what
`.github/workflows/java.yml` does. Pushing a `v*` tag publishes the jar to
Maven Central, with the version of the tag (which must match the Python
package's `python/mathoptmodel/_version.py`). At run time, the library is loaded from the
`mathoptmodel.library.path` system property if set, otherwise from the jar,
otherwise from the system's library path.

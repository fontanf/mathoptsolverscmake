# Model-building benchmark

Times building a model, then writing it as an LP file, with mathoptmodel and
other modeling tools:

- Python: mathoptmodel, Pyomo and pyoptinterface;
- Java: mathoptmodel and OR-Tools.

The model is a random sparse MILP with `n` variables (half integer, half
continuous, all in `[0, 10]`), `n` `<=` constraints with 5 nonzeros each, and
a linear objective to maximize; see `data.py`. Python and Java generate
exactly the same model, without sharing data files. Only linear models are
compared, since OR-Tools has no nonlinear modeling.

Each tool builds the model the way it's meant to be used:

| Tool | How the model is built |
|---|---|
| `mathoptmodel-batch` | one `add_variables` / `add_constraints` call with all the arrays (CSR) |
| `mathoptmodel-row` | one `add_constraints` call per constraint |
| `pyomo` | `ConcreteModel`, an indexed `Var`, a rule-based `Constraint` with `sum(...)` expressions |
| `pyomo-linear-expression` | the same, with `LinearExpression(linear_coefs=..., linear_vars=...)` instead of `sum(...)` |
| `pyoptinterface-row` | `add_variable` per variable, `add_linear_constraint(ScalarAffineFunction(...))` per constraint, into HiGHS |
| `pyoptinterface-matrix` | `add_m_variables` and `add_m_linear_constraints` with a SciPy CSR matrix, into HiGHS |
| `ortools-modelbuilder` | OR-Tools' solver-independent `ModelBuilder`, one `addLessOrEqual(LinearExpr.weightedSum(...))` per constraint |
| `ortools-modelbuilder-helper` | the same model through `ModelBuilder.getHelper()`: `addVar`, `addLinearConstraint`, `addConstraintTerm`... one native call per attribute and per term |
| `ortools-mpsolver` | OR-Tools' `MPSolver` (SCIP backend, not solved), `makeConstraint` + `setCoefficient` per term |
| `ortools-mpsolver-proto` | an `MPModelProto` built in Java, then `MPSolver.loadModelFromProto` |

The `pyomo`, `ortools-modelbuilder` and `ortools-mpsolver` rows use each
tool's usual API; the other variants of these tools are the fastest ways we
found to build the same model with them.

Timings exclude generating the data, and, for the tools that add constraints
one at a time, converting it to the tool's input types (lists, rows) too.
Variables and constraints have no names. Each size runs in its own process
under a memory cap (`systemd-run`, 4 GB by default), and reports the best time
over 5, 3 or 2 repetitions (10k, 100k, 1M); in Java, the first repetitions
also warm the JIT up. The LP files written for the check size (10k by
default) are read back with mathoptmodel and compared with the generated
data.

## Running

Python tools, in a virtual environment:

```shell
python3 -m venv venv
venv/bin/pip install mathoptmodel pyomo pyoptinterface highsbox numpy scipy
```

(`highsbox` provides the HiGHS library pyoptinterface loads.)

Java tools, with Java 22+ and Maven:

```shell
cd java
mvn compile dependency:build-classpath -Dmdep.outputFile=target/classpath.txt
cd ..
```

Then:

```shell
venv/bin/python run.py                                   # everything
venv/bin/python run.py --sizes 10000 100000 --languages python
venv/bin/python run.py --tools pyomo java:ortools-mpsolver --memory 8G
```

The table is printed and saved to `results/results.md` (and the raw numbers
to `results/results.json`). A single run can also be launched directly, for
example `venv/bin/python bench_python.py pyomo 100000 3 /tmp` or
`java -cp ... bench.Bench ortools-modelbuilder 100000 3 /tmp`.

# MathOptSolversCMake

This library includes:
* CMake wrappers for mathematical optimization solvers
* A mathematical programming modeler that supports:
  * Continuous, integer, semi-continuous and semi-integer variables (the last two only with HiGHS, CPLEX and Gurobi)
  * Linear structures
  * Quadratic structures
  * Nonlinear structures
  * Black-box functions

The goal of the modeler are:
* Minimize the modeler's overhead
* Run multiple solvers while writing the model's code once and ensuring that the model passed to each solver is the same
* Keep access to all the direct API features of the solvers
* Minimize the quantity of code to integrate a new solver
* Provide some features to help model debugging

They are not designed to be as user-friendly as possible.
And switching solver requires a bit more lines of code than changing a string.

Supported solvers:
* HiGHS (MILP) https://highs.dev/
* Cbc (MILP) https://github.com/coin-or/Cbc
* FICO Xpress (MILP) https://www.fico.com/en/products/fico-xpress-optimization
* IBM ILOG CPLEX (MIQP; non-convex quadratic objectives solved to global optimality) https://www.ibm.com/products/ilog-cplex-optimization-studio
* Gurobi (MINLP, solved to global optimality; no black-box functions) https://www.gurobi.com/
* Artelys Knitro (all) https://www.artelys.com/solvers/knitro/
* SCIP (MINLP, solved to global optimality; no black-box functions) https://www.scipopt.org/
* Ipopt (continuous NLP, local optimality) https://github.com/coin-or/Ipopt
* Dlib (box-constrained) https://dlib.net/
* ConicBundle (box-constrained) https://www-user.tu-chemnitz.de/~helmberg/ConicBundle/

Examples:
* MILP:
  * [Multiple-choice knapsack](https://github.com/fontanf/multiplechoiceknapsacksolver/blob/main/src/algorithms/milp.cpp)
  * [Set covering](https://github.com/fontanf/setcoveringsolver/blob/master/src/algorithms/milp.cpp)
  * [Generalized assignment](https://github.com/fontanf/generalizedassignmentsolver/blob/master/src/algorithms/milp.cpp)
  * [Clique](https://github.com/fontanf/stablesolver/blob/master/src/clique/algorithms/milp.cpp), [stable](https://github.com/fontanf/stablesolver/blob/master/src/stable/algorithms/milp.cpp)
  * [Knapsack with conflicts](https://github.com/fontanf/knapsackwithconflictssolver/blob/main/src/algorithms/milp.cpp)
  * [Graph coloring](https://github.com/fontanf/coloringsolver/blob/master/src/algorithms/milp.cpp)
  * Shop scheduling, [positional model](https://github.com/fontanf/shopschedulingsolver/blob/main/src/algorithms/milp_positional.cpp), [disjunctive model](https://github.com/fontanf/shopschedulingsolver/blob/main/src/algorithms/milp_disjunctive.cpp)
* Box-constrained (Lagrangian relaxations):
  * [Generalized assignment](https://github.com/fontanf/generalizedassignmentsolver/blob/master/src/algorithms/lagrangian_relaxation.cpp)
  * [Knapsack with conflicts](https://github.com/fontanf/knapsackwithconflictssolver/blob/main/src/algorithms/lagrangian_relaxation.cpp)

CMake integration example:
```cmake
# Fetch fontanf/mathoptsolverscmake.
set(MATHOPTSOLVERSCMAKE_USE_CLP ON)
FetchContent_Declare(
    mathoptsolverscmake
    GIT_REPOSITORY https://github.com/fontanf/mathoptsolverscmake.git
    GIT_TAG ...)
    #SOURCE_DIR "${PROJECT_SOURCE_DIR}/../mathoptsolverscmake/")
FetchContent_MakeAvailable(mathoptsolverscmake)

...

target_link_libraries(MyProject_my_target PUBLIC
    MathOptSolversCMake::clp)
```

## CLI app

An optional CLI app reads a model file and solves it directly, no code required.
It supports the MPS, LP, JSON and .nl formats (guessed from the file extension, or given with `--format`),
and the Cbc, HiGHS, XPRESS, CPLEX, Gurobi, Knitro, Ipopt and SCIP solvers (`--solver`), whichever were enabled at configure time.
Enable it with `MATHOPTSOLVERSCMAKE_BUILD_APP`, together with the solver(s) it should be able to use:
```shell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DMATHOPTSOLVERSCMAKE_USE_HIGHS=ON
cmake --build build --config Release --parallel
cmake --install build --config Release --prefix install
```

For example, solving the small MIPLIB3 instance [`data/miplib3/p0033.mps`](data/miplib3/p0033.mps)
(33 binary variables, optimal objective value 3089) with HiGHS:
```shell
./install/bin/mathoptsolverscmake_solve --input data/miplib3/p0033.mps --solver highs
```
```
===========================
          MathOpt          
===========================

Model
-----
Number of variables:    33
Number of constraints:  16
Number of elements:     98
Objective:              Minimize
Has non-continuous:     1
Has quadratic:          0
Has nonlinear:          0
Has black-box:          0
Is LP:                  0
Is MILP:                1
Is box-constrained:     0

Running HiGHS 1.16.0 (git hash: n/a): Copyright (c) 2026 under MIT licence terms
Includes third-party software components, see THIRD_PARTY_NOTICES.md for full details
MIP has 16 rows; 33 cols; 98 nonzeros; 33 integer variables (33 binary)
Coefficient ranges:
  Matrix  [1e+00, 4e+02]
  Cost    [5e+01, 5e+02]
  Bound   [1e+00, 1e+00]
  RHS     [1e+00, 3e+03]
Presolving model
16 rows, 33 cols, 98 nonzeros 0s
13 rows, 26 cols, 67 nonzeros 0s
12 rows, 24 cols, 64 nonzeros 0s
Presolve reductions: rows 12(-4); columns 24(-9); nonzeros 64(-34) 
Objective function is integral with scale 1

Solving MIP model with:
   12 rows
   24 cols (20 binary, 4 integer, 0 implied int., 0 continuous, 0 domain fixed)
   64 nonzeros
   Thread count: 6 available cores, 6 logical processors, using up to 6 workers across 6 threads.

Src: B => Branching; C => Central rounding; F => Feasibility pump; H => Heuristic;
     I => Shifting; J => Feasibility jump; L => Sub-MIP; P => Empty MIP; R => Randomized rounding;
     S => Solve LP; T => Evaluate node; U => Unbounded; X => User solution; Y => HiGHS solution;
     Z => ZI Round; l => Trivial lower; p => Trivial point; u => Trivial upper; z => Trivial zero

        Nodes      |    B&B Tree     |            Objective Bounds              |  Dynamic Constraints |       Work      
Src  Proc. InQueue |  Leaves   Expl. | BestBound       BestSol              Gap |   Cuts   InLp Confl. | LpIters     Time

 J       0       0         0   0.00%   -inf            3705               Large        0      0      0         0     0.0s
 R       0       0         0   0.00%   2838.546739     3347              15.19%        0      0      0        10     0.0s

4.2% inactive integer columns, restarting
Model after restart has 12 rows, 23 cols (19 bin., 4 int., 0 impl., 0 cont., 0 dom.fix.), and 61 nonzeros

         0       0         0   0.00%   3073.628205     3347               8.17%        9      0      0        55     0.0s
         0       0         0   0.00%   3073.628205     3347               8.17%        9      8      7        69     0.0s
 L       0       0         0   0.00%   3089            3089               0.00%      165     16      7        85     0.0s
         1       0         1 100.00%   3089            3089               0.00%      165     16      7        85     0.0s

Solving report
  Status            Optimal
  Primal bound      3089
  Dual bound        3089
  Gap               0% (tolerance: 0.01%)
  P-D integral      0.01
  Solution status   feasible
                    3089 (objective)
  Timing            0.03
                    0.00 (Presolve)
                        MIP    time [calls] = 0.00 [1]
                        subMIP time [calls] = 0.00 [3]
                    0.03 (Solve)
                        MIP    time [calls] = 0.02 [1]
                        subMIP time [calls] = 0.01 [1]
                    0.00 (Postsolve)
                        MIP    time [calls] = 0.00 [1]
                        subMIP time [calls] = 0.00 [3]
  Max sub-MIP depth 2
  Nodes             1
  Repair LPs        0
  LP iterations     85
                    0 (strong br.)
                    54 (separation)
                    7 (heuristics)

Final statistics
----------------
Objective value:  3089
Bound:            3089
Feasible:         1
```

For a linear program, `--dual` solves its dual (built with `mathoptsolverscmake::dual()`) instead of the model itself.
The reported objective value is the dual's, equal to the primal's at optimality, and `--output` writes the dual solution.

## Python package

The solver-independent part of the library — `MathOptModel` and the MPS, LP, JSON and `.nl` readers/writers — is also available as a Python package, `mathoptmodel`, with the model's arrays exposed as read-only NumPy views, and models buildable from Python in batches of NumPy arrays:
```shell
pip install mathoptmodel
```
```python
import mathoptmodel as mom

model = mom.read("data/miplib3/p0033.mps")  # also .gz / .bz2
print(model.number_of_variables, model.constraint_class(0))
```
See [`python/README.md`](python/README.md) for details. The wheels are built by `.github/workflows/python.yml` and published to PyPI when a `v*` tag is pushed.

## Java bindings

The same part of the library is available from Java 22+, through the C API `include/mathoptsolverscmake/mathoptmodel_c.h` (the `mathoptmodel` shared library, built with `-DMATHOPTSOLVERSCMAKE_BUILD_C_API=ON`) and the Foreign Function & Memory API:
```java
try (Model model = Model.read(Path.of("data/miplib3/p0033.mps"))) {
    System.out.println(model.numberOfVariables() + " " + model.constraintClass(0));
}
```
See [`java/README.md`](java/README.md) for details. `.github/workflows/java.yml` builds the native library for each platform, tests the bindings and packages a jar bundling all of them.

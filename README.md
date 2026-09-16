# MathOptSolversCMake

This library includes:
* CMake wrappers for mathematical optimization solvers
* A mathematical programming modeler that supports:
  * Coninous and integer variables
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
* Artelys Knitro (all) https://www.artelys.com/solvers/knitro/
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
and the Cbc, HiGHS, XPRESS and Knitro solvers (`--solver`), whichever were enabled at configure time.
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

Running HiGHS 1.12.0 (git hash: n/a): Copyright (c) 2025 HiGHS under MIT licence terms
MIP has 16 rows; 33 cols; 98 nonzeros; 33 integer variables (33 binary)
Coefficient ranges:
  Matrix  [1e+00, 4e+02]
  Cost    [5e+01, 5e+02]
  Bound   [1e+00, 1e+00]
  RHS     [1e+00, 3e+03]
Presolving model
15 rows, 32 cols, 97 nonzeros  0s
14 rows, 26 cols, 68 nonzeros  0s
Presolve reductions: rows 14(-2); columns 26(-7); nonzeros 68(-30) 
Objective function is integral with scale 1

Solving MIP model with:
   14 rows
   26 cols (22 binary, 4 integer, 0 implied int., 0 continuous, 0 domain fixed)
   68 nonzeros

Src: B => Branching; C => Central rounding; F => Feasibility pump; H => Heuristic;
     I => Shifting; J => Feasibility jump; L => Sub-MIP; P => Empty MIP; R => Randomized rounding;
     S => Solve LP; T => Evaluate node; U => Unbounded; X => User solution; Y => HiGHS solution;
     Z => ZI Round; l => Trivial lower; p => Trivial point; u => Trivial upper; z => Trivial zero

        Nodes      |    B&B Tree     |            Objective Bounds              |  Dynamic Constraints |       Work      
Src  Proc. InQueue |  Leaves   Expl. | BestBound       BestSol              Gap |   Cuts   InLp Confl. | LpIters     Time

 J       0       0         0   0.00%   -inf            3596               Large        0      0      0         0     0.0s
 R       0       0         0   0.00%   2838.546739     3347              15.19%        0      0      0        10     0.0s
 L       0       0         0   0.00%   3081.347826     3164               2.61%      209     18     25        42     0.0s
 H       0       0         0   0.00%   3086            3089               0.10%      237     19     25        46     0.0s
         1       0         1 100.00%   3089            3089               0.00%      238     19     28        56     0.0s

Solving report
  Status            Optimal
  Primal bound      3089
  Dual bound        3089
  Gap               0% (tolerance: 0.01%)
  P-D integral      0.00220620091008
  Solution status   feasible
                    3089 (objective)
                    0 (bound viol.)
                    5.59996493621e-13 (int. viol.)
                    0 (row viol.)
  Timing            0.02
  Max sub-MIP depth 1
  Nodes             1
  Repair LPs        0
  LP iterations     56
                    0 (strong br.)
                    34 (separation)
                    12 (heuristics)

Final statistics
----------------
Objective value:  3089
Bound:            3089
Feasible:         1
```

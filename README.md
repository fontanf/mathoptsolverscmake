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

cmake_minimum_required(VERSION 3.15.0)

project(MathOptSolversCMake LANGUAGES CXX)

# Solver options.
option(MATHOPTSOLVERSCMAKE_USE_CLP "Use Clp" OFF)
option(MATHOPTSOLVERSCMAKE_USE_CBC "Use Cbc" OFF)
option(MATHOPTSOLVERSCMAKE_USE_HIGHS "Use HiGHS" OFF)
option(MATHOPTSOLVERSCMAKE_USE_XPRESS "Use XPRESS" OFF)
option(MATHOPTSOLVERSCMAKE_USE_CPLEX "Use CPLEX" OFF)
option(MATHOPTSOLVERSCMAKE_USE_GUROBI "Use GUROBI" OFF)
option(MATHOPTSOLVERSCMAKE_USE_KNITRO "Use Knitro" OFF)
option(MATHOPTSOLVERSCMAKE_USE_DLIB "Use dlib" OFF)
option(MATHOPTSOLVERSCMAKE_USE_CONICBUNDLE "Use ConicBundle" OFF)

# Avoid FetchContent warning.
cmake_policy(SET CMP0135 NEW)

# Require C++11.
set(CMAKE_CXX_STANDARD 11)

# Enable output of compile commands during generation.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Add sub-directories.
add_subdirectory(extern)
add_subdirectory(src)

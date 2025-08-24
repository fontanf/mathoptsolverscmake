cmake_minimum_required(VERSION 3.15.0)

project(MathOptSolversCMake LANGUAGES CXX)

# Solver options.
option(MATHOPTSOLVERSCMAKE_USE_CLP "Use Clp" ON)
option(MATHOPTSOLVERSCMAKE_USE_CBC "Use Cbc" ON)
option(MATHOPTSOLVERSCMAKE_USE_HIGHS "Use HiGHS" ON)
option(MATHOPTSOLVERSCMAKE_USE_XPRESS "Use XPRESS" ON)
option(MATHOPTSOLVERSCMAKE_USE_CPLEX "Use CPLEX" OFF)
option(MATHOPTSOLVERSCMAKE_USE_GUROBI "Use GUROBI" OFF)

# Avoid FetchContent warning.
cmake_policy(SET CMP0135 NEW)

# Require C++11.
set(CMAKE_CXX_STANDARD 11)

# Enable output of compile commands during generation.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Add sub-directories.
add_subdirectory(extern)
add_subdirectory(src)

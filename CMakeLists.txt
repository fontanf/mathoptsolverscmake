cmake_minimum_required(VERSION 3.15.0)

project(MathOptSolversCMake LANGUAGES CXX)

option(MATHOPTSOLVERSCMAKE_USE_CLP "Use CLP" OFF)
option(MATHOPTSOLVERSCMAKE_USE_CBC "Use CBC" OFF)
option(MATHOPTSOLVERSCMAKE_USE_XPRESS "Use XPRESS" OFF)
option(MATHOPTSOLVERSCMAKE_USE_CPLEX "Use CPLEX" OFF)
option(MATHOPTSOLVERSCMAKE_USE_GUROBI "Use GUROBI" OFF)

# Require C++11.
set(CMAKE_CXX_STANDARD 11)

# Enable output of compile commands during generation.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Add sub-directories.
add_subdirectory(extern)

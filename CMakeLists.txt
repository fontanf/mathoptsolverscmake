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
option(MATHOPTSOLVERSCMAKE_USE_IPOPT "Use Ipopt" OFF)
option(MATHOPTSOLVERSCMAKE_USE_SCIP "Use SCIP" OFF)
option(MATHOPTSOLVERSCMAKE_SCIP_USE_IPOPT "Build SCIP with Ipopt as NLP solver" ON)
option(MATHOPTSOLVERSCMAKE_USE_DLIB "Use dlib" OFF)
option(MATHOPTSOLVERSCMAKE_USE_CONICBUNDLE "Use ConicBundle" OFF)
option(MATHOPTSOLVERSCMAKE_USE_SMSPP "Use SMS++ (BundleSolver)" OFF)

# Test option.
option(MATHOPTSOLVERSCMAKE_BUILD_TESTS "Build tests" ON)

# App option.
option(MATHOPTSOLVERSCMAKE_BUILD_APP "Build the mathoptsolverscmake_solve CLI app" ON)

# Python bindings option (set by pyproject.toml when building the Python package).
option(MATHOPTSOLVERSCMAKE_BUILD_PYTHON "Build the Python bindings" OFF)

# C API option: the "mathoptmodel" shared library, used by the Java bindings.
option(MATHOPTSOLVERSCMAKE_BUILD_C_API "Build the C API shared library" OFF)

# Avoid FetchContent warning.
cmake_policy(SET CMP0135 NEW)

# Require C++11.
set(CMAKE_CXX_STANDARD 11)

# Enable output of compile commands during generation.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# The static library is linked into the Python extension module and into the
# C API shared library.
if(MATHOPTSOLVERSCMAKE_BUILD_PYTHON OR MATHOPTSOLVERSCMAKE_BUILD_C_API)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Add sub-directories.
add_subdirectory(extern)
add_subdirectory(src)
if(MATHOPTSOLVERSCMAKE_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()
if(MATHOPTSOLVERSCMAKE_BUILD_APP)
    add_subdirectory(app)
endif()
if(MATHOPTSOLVERSCMAKE_BUILD_PYTHON)
    add_subdirectory(python)
endif()

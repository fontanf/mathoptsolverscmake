find_package(Python 3.8
    COMPONENTS Interpreter Development.Module
    OPTIONAL_COMPONENTS Development.SABIModule
    REQUIRED)

# Fetch nanobind.
include(FetchContent)
FetchContent_Declare(
    nanobind
    GIT_REPOSITORY https://github.com/wjakob/nanobind.git
    GIT_TAG v3.1.0
    EXCLUDE_FROM_ALL)
FetchContent_MakeAvailable(nanobind)

# nanobind requires C++17; only the binding sources are compiled with it, the
# library itself stays C++11.
nanobind_add_module(_core
    NB_STATIC
    STABLE_ABI
    src/mathoptmodel_ext.cpp)
target_compile_features(_core PRIVATE cxx_std_17)
target_link_libraries(_core PRIVATE
    MathOptSolversCMake::mathopt)

install(TARGETS _core LIBRARY DESTINATION mathoptmodel)

# HiGHS's license, next to the bundled HiGHS.
if(MATHOPTSOLVERSCMAKE_USE_HIGHS)
    FetchContent_GetProperties(highs)
    install(FILES "${highs_SOURCE_DIR}/LICENSE.txt"
        DESTINATION mathoptmodel
        RENAME LICENSE-HiGHS.txt)
endif()

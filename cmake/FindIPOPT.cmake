# Find module used by SCIP's find_package(IPOPT) (see
# mathoptsolverscmake_fetch_scip in extern/CMakeLists.txt) so that SCIP links
# the Ipopt fetched by MathOptSolversCMake instead of searching the system.
#
# Defines the variables SCIP's CMakeLists.txt uses: IPOPT_FOUND,
# IPOPT_INCLUDE_DIRS, IPOPT_LIBRARIES and IPOPT_VERSION.

if(TARGET MathOptSolversCMake_ipopt)
    get_target_property(IPOPT_INCLUDE_DIRS MathOptSolversCMake_ipopt INTERFACE_INCLUDE_DIRECTORIES)
    set(IPOPT_LIBRARIES MathOptSolversCMake_ipopt)
    if(EXISTS "${IPOPT_INCLUDE_DIRS}/IpoptConfig.h")
        file(STRINGS "${IPOPT_INCLUDE_DIRS}/IpoptConfig.h" _ipopt_version_line
            REGEX "^#define IPOPT_VERSION ")
        string(REGEX REPLACE "^#define IPOPT_VERSION \"([^\"]*)\".*" "\\1"
            IPOPT_VERSION "${_ipopt_version_line}")
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(IPOPT
    REQUIRED_VARS IPOPT_LIBRARIES IPOPT_INCLUDE_DIRS
    VERSION_VAR IPOPT_VERSION)

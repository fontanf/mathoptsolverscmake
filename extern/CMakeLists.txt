# Enable FetchContent.
include(FetchContent)

# COIN-OR CLP/CBC

if(MATHOPTSOLVERSCMAKE_USE_CLP OR MATHOPTSOLVERSCMAKE_USE_CBC OR MATHOPTSOLVERSCMAKE_USE_SMSPP)

    # CoinUtils

    add_library(MathOptSolversCMake_coinutils STATIC IMPORTED GLOBAL)

    if(UNIX AND NOT APPLE)
        FetchContent_Declare(coinutils
            URL https://github.com/coin-or/CoinUtils/releases/download/releases%2F2.11.9/CoinUtils-releases.2.11.9-x86_64-ubuntu20-gcc940-static.tar.gz
            URL_HASH SHA256=14d07de1b7961f68e037da6f0c57844fd67d4cc1a4b125642f42cd134b228094
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            OVERRIDE_FIND_PACKAGE)
        FetchContent_MakeAvailable(coinutils)

        target_include_directories(MathOptSolversCMake_coinutils INTERFACE
            "${coinutils_SOURCE_DIR}/include")
        set_target_properties(MathOptSolversCMake_coinutils PROPERTIES
            IMPORTED_LOCATION "${coinutils_SOURCE_DIR}/lib/libCoinUtils.a")
        find_package(LAPACK)
        find_package(BLAS)
        find_package(BZip2)
        find_package(ZLIB)
        target_link_libraries(MathOptSolversCMake_coinutils INTERFACE
            LAPACK::LAPACK
            BLAS::BLAS
            BZip2::BZip2
            ZLIB::ZLIB)
    endif()

    if(APPLE)
        FetchContent_Declare(coinutils
            URL https://github.com/fontanf/coinor-builds/releases/download/v0.1/CoinUtils-release.2.11.12-universal-macos.tar.gz
            URL_HASH SHA256=8cd0b09e2eed185b5b37f81ef66554e0d41f29e9b5ea87b128538ce158fa6bfa
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            OVERRIDE_FIND_PACKAGE)
        FetchContent_MakeAvailable(coinutils)

        target_include_directories(MathOptSolversCMake_coinutils INTERFACE
            "${coinutils_SOURCE_DIR}/include")
        set_target_properties(MathOptSolversCMake_coinutils PROPERTIES
            IMPORTED_LOCATION "${coinutils_SOURCE_DIR}/lib/libCoinUtils.a")
    endif()

    if(WIN32)
        FetchContent_Declare(coinutils
            URL https://github.com/coin-or/CoinUtils/releases/download/releases%2F2.11.9/CoinUtils-releases.2.11.9-w64-msvc16-md.zip
            URL_HASH SHA256=2bc64f0afd80113571697e949b2663e9047272decf90d5f62e452c2628d33ca6
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            OVERRIDE_FIND_PACKAGE)
        FetchContent_MakeAvailable(coinUtils)

        target_include_directories(MathOptSolversCMake_coinutils INTERFACE
            "${coinutils_SOURCE_DIR}/include")
        set_target_properties(MathOptSolversCMake_coinutils PROPERTIES
            IMPORTED_LOCATION "${coinutils_SOURCE_DIR}/lib/libCoinUtils.lib")
    endif()

    add_dependencies(MathOptSolversCMake_coinutils coinutils)
    add_library(MathOptSolversCMake::coinutils ALIAS MathOptSolversCMake_coinutils)

    # Osi

    add_library(MathOptSolversCMake_osi STATIC IMPORTED GLOBAL)

    if(UNIX AND NOT APPLE)
        FetchContent_Declare(osi
            URL https://github.com/coin-or/Osi/releases/download/releases%2F0.108.8/Osi-releases.0.108.8-x86_64-ubuntu20-gcc940-static.tar.gz
            URL_HASH SHA256=bd5a5bf1e6b6a28d13d41ab1554becd9f3992afe775785e51a88c9405cf2853e
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            OVERRIDE_FIND_PACKAGE)
        FetchContent_MakeAvailable(osi)

        add_library(MathOptSolversCMake_osicommontests STATIC IMPORTED GLOBAL)
        target_include_directories(MathOptSolversCMake_osicommontests INTERFACE
            "${osi_SOURCE_DIR}/include")
        set_target_properties(MathOptSolversCMake_osicommontests PROPERTIES
            IMPORTED_LOCATION "${osi_SOURCE_DIR}/lib/libOsiCommonTests.a")

        target_include_directories(MathOptSolversCMake_osi INTERFACE
            "${osi_SOURCE_DIR}/include")
        target_link_directories(MathOptSolversCMake_osi INTERFACE
            "${osi_SOURCE_DIR}/lib")
        set_target_properties(MathOptSolversCMake_osi PROPERTIES
            IMPORTED_LOCATION "${osi_SOURCE_DIR}/lib/libOsi.a")
        target_link_libraries(MathOptSolversCMake_osi INTERFACE
            MathOptSolversCMake_osicommontests)
    endif()

    if(APPLE)
        FetchContent_Declare(osi
            URL https://github.com/fontanf/coinor-builds/releases/download/v0.1/Osi-release.0.108.11-universal-macos.tar.gz
            URL_HASH SHA256=fc6b7b93106e0d4d0d8eea1b2fb9f9db3c345580444409655133fc984391682b
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            OVERRIDE_FIND_PACKAGE)
        FetchContent_MakeAvailable(osi)

        add_library(MathOptSolversCMake_osicommontests STATIC IMPORTED GLOBAL)
        target_include_directories(MathOptSolversCMake_osicommontests INTERFACE
            "${osi_SOURCE_DIR}/include")
        set_target_properties(MathOptSolversCMake_osicommontests PROPERTIES
            IMPORTED_LOCATION "${osi_SOURCE_DIR}/lib/libOsiCommonTests.a")

        target_include_directories(MathOptSolversCMake_osi INTERFACE
            "${osi_SOURCE_DIR}/include")
        target_link_directories(MathOptSolversCMake_osi INTERFACE
            "${osi_SOURCE_DIR}/lib")
        set_target_properties(MathOptSolversCMake_osi PROPERTIES
            IMPORTED_LOCATION "${osi_SOURCE_DIR}/lib/libOsi.a")
        target_link_libraries(MathOptSolversCMake_osi INTERFACE
            MathOptSolversCMake_osicommontests)
    endif()

    if(WIN32)
        FetchContent_Declare(osi
            URL https://github.com/coin-or/Osi/releases/download/releases%2F0.108.8/Osi-releases.0.108.8-w64-msvc16-md.zip
            URL_HASH SHA256=a61fc462cb598139d205cd2323522581a01900575d0d6bccf660b6c7e1b0b71c
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            OVERRIDE_FIND_PACKAGE)
        FetchContent_MakeAvailable(osi)

        add_library(MathOptSolversCMake_osicommontests STATIC IMPORTED GLOBAL)
        target_include_directories(MathOptSolversCMake_osicommontests INTERFACE
            "${osi_SOURCE_DIR}/include")
        set_target_properties(MathOptSolversCMake_osicommontests PROPERTIES
            IMPORTED_LOCATION "${osi_SOURCE_DIR}/lib/libOsiCommonTests.lib")

        target_include_directories(MathOptSolversCMake_osi INTERFACE
            "${osi_SOURCE_DIR}/include")
        target_link_directories(MathOptSolversCMake_osi INTERFACE
            "${osi_SOURCE_DIR}/lib")
        set_target_properties(MathOptSolversCMake_osi PROPERTIES
            IMPORTED_LOCATION "${osi_SOURCE_DIR}/lib/libOsi.lib")
        target_link_libraries(MathOptSolversCMake_osi INTERFACE
            MathOptSolversCMake_osicommontests)
    endif()

    add_dependencies(MathOptSolversCMake_osi osi)
    target_link_libraries(MathOptSolversCMake_osi INTERFACE
        MathOptSolversCMake_coinutils)
    add_library(MathOptSolversCMake::osi ALIAS MathOptSolversCMake_osi)

    # Clp

    add_library(MathOptSolversCMake_clp STATIC IMPORTED GLOBAL)
    if(UNIX AND NOT APPLE)
        FetchContent_Declare(clp
            URL https://github.com/coin-or/Clp/releases/download/releases%2F1.17.8/Clp-releases.1.17.8-x86_64-ubuntu20-gcc940-static.tar.gz
            URL_HASH SHA256=d569b04d19c25876e55d2557a1d9739df8eb50ec8ca11a98ad387fd8b90212c9
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            OVERRIDE_FIND_PACKAGE)
        FetchContent_MakeAvailable(clp)

        add_library(MathOptSolversCMake_osiclp STATIC IMPORTED GLOBAL)
        target_include_directories(MathOptSolversCMake_osiclp INTERFACE
            "${clp_SOURCE_DIR}/include")
        set_target_properties(MathOptSolversCMake_osiclp PROPERTIES
            IMPORTED_LOCATION "${clp_SOURCE_DIR}/lib/libOsiClp.a")

        target_include_directories(MathOptSolversCMake_clp INTERFACE
            "${clp_SOURCE_DIR}/include")
        set_target_properties(MathOptSolversCMake_clp PROPERTIES
            IMPORTED_LOCATION "${clp_SOURCE_DIR}/lib/libClp.a")
        target_link_libraries(MathOptSolversCMake_clp INTERFACE
            MathOptSolversCMake_osiclp)
    endif()

    if(APPLE)
        FetchContent_Declare(clp
            URL https://github.com/fontanf/coinor-builds/releases/download/v0.1/Clp-release.1.17.10-universal-macos.tar.gz
            URL_HASH SHA256=1179ae58b6d0fe827ff8933f72eed2a674273968d12e8125a93f5b03f9e8757b
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            OVERRIDE_FIND_PACKAGE)
        FetchContent_MakeAvailable(clp)

        add_library(MathOptSolversCMake_osiclp STATIC IMPORTED GLOBAL)
        target_include_directories(MathOptSolversCMake_osiclp INTERFACE
            "${clp_SOURCE_DIR}/include")
        set_target_properties(MathOptSolversCMake_osiclp PROPERTIES
            IMPORTED_LOCATION "${clp_SOURCE_DIR}/lib/libOsiClp.a")

        target_include_directories(MathOptSolversCMake_clp INTERFACE
            "${clp_SOURCE_DIR}/include")
        set_target_properties(MathOptSolversCMake_clp PROPERTIES
            IMPORTED_LOCATION "${clp_SOURCE_DIR}/lib/libClp.a")
        target_link_libraries(MathOptSolversCMake_clp INTERFACE
            MathOptSolversCMake_osiclp)
    endif()

    if(WIN32)
        FetchContent_Declare(clp
            URL https://github.com/coin-or/Clp/releases/download/releases%2F1.17.8/Clp-releases.1.17.8-w64-msvc16-md.zip
            URL_HASH SHA256=e37c834aea5c31dfd8620b7d2432cb31fc16ecb0c6ffb398e8f07c9c82dd5028
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            OVERRIDE_FIND_PACKAGE)
        FetchContent_MakeAvailable(clp)

        add_library(MathOptSolversCMake_osiclp STATIC IMPORTED GLOBAL)
        target_include_directories(MathOptSolversCMake_osiclp INTERFACE
            "${clp_SOURCE_DIR}/include")
        set_target_properties(MathOptSolversCMake_osiclp PROPERTIES
            IMPORTED_LOCATION "${clp_SOURCE_DIR}/lib/libOsiClp.lib")

        target_include_directories(MathOptSolversCMake_clp INTERFACE
            "${clp_SOURCE_DIR}/include")
        set_target_properties(MathOptSolversCMake_clp PROPERTIES
            IMPORTED_LOCATION "${clp_SOURCE_DIR}/lib/libClp.lib")
        target_link_libraries(MathOptSolversCMake_clp INTERFACE
            MathOptSolversCMake_osiclp)
    endif()

    add_dependencies(MathOptSolversCMake_clp clp)
    target_link_libraries(MathOptSolversCMake_clp INTERFACE
        MathOptSolversCMake_osi)
    add_library(MathOptSolversCMake::clp ALIAS MathOptSolversCMake_clp)
endif()

# Cbc

if(MATHOPTSOLVERSCMAKE_USE_CBC)

    add_library(MathOptSolversCMake_cbc STATIC IMPORTED GLOBAL)

    if(UNIX AND NOT APPLE)
        FetchContent_Declare(cbc
            URL https://github.com/coin-or/Cbc/releases/download/releases%2F2.10.10/Cbc-releases.2.10.10-x86_64-ubuntu20-gcc940-static.tar.gz
            URL_HASH SHA256=872c78bfcdd1566f134d2f7757b76b2a2479a5b1ade065cdd1d4b303ed6f8006
            CONFIGURE_COMMAND ""
            BUILD_COMMAND "")
        FetchContent_MakeAvailable(cbc)

        add_library(MathOptSolversCMake_osicbc STATIC IMPORTED GLOBAL)
        target_include_directories(MathOptSolversCMake_osicbc INTERFACE
            "${osi_SOURCE_DIR}/include")
        set_target_properties(MathOptSolversCMake_osicbc PROPERTIES
            IMPORTED_LOCATION "${cbc_SOURCE_DIR}/lib/libOsiCbc.a")

        add_library(MathOptSolversCMake_cgl STATIC IMPORTED GLOBAL)
        target_include_directories(MathOptSolversCMake_cgl INTERFACE
            "${osi_SOURCE_DIR}/include")
        set_target_properties(MathOptSolversCMake_cgl PROPERTIES
            IMPORTED_LOCATION "${cbc_SOURCE_DIR}/lib/libCgl.a")

        target_include_directories(MathOptSolversCMake_cbc INTERFACE
            "${cbc_SOURCE_DIR}/include")
        set_target_properties(MathOptSolversCMake_cbc PROPERTIES
            IMPORTED_LOCATION "${cbc_SOURCE_DIR}/lib/libCbc.a")
        target_link_libraries(MathOptSolversCMake_cbc INTERFACE
            nauty
            MathOptSolversCMake_osicbc
            MathOptSolversCMake_cgl)
    endif()

    if(APPLE)
        FetchContent_Declare(cbc
            URL https://github.com/fontanf/coinor-builds/releases/download/v0.1/Cbc-release.2.10.12-universal-macos.tar.gz
            URL_HASH SHA256=1bfbcf30a91561857b15b7b541daddd986007eb7b1bec66b1af3f39a0e139bc8
            CONFIGURE_COMMAND ""
            BUILD_COMMAND "")
        FetchContent_MakeAvailable(cbc)

        add_library(MathOptSolversCMake_osicbc STATIC IMPORTED GLOBAL)
        target_include_directories(MathOptSolversCMake_osicbc INTERFACE
            "${osi_SOURCE_DIR}/include")
        set_target_properties(MathOptSolversCMake_osicbc PROPERTIES
            IMPORTED_LOCATION "${cbc_SOURCE_DIR}/lib/libOsiCbc.a")

        add_library(MathOptSolversCMake_cgl STATIC IMPORTED GLOBAL)
        target_include_directories(MathOptSolversCMake_cgl INTERFACE
            "${osi_SOURCE_DIR}/include")
        set_target_properties(MathOptSolversCMake_cgl PROPERTIES
            IMPORTED_LOCATION "${cbc_SOURCE_DIR}/lib/libCgl.a")

        target_include_directories(MathOptSolversCMake_cbc INTERFACE
            "${cbc_SOURCE_DIR}/include")
        set_target_properties(MathOptSolversCMake_cbc PROPERTIES
            IMPORTED_LOCATION "${cbc_SOURCE_DIR}/lib/libCbc.a")
        target_link_libraries(MathOptSolversCMake_cbc INTERFACE
            MathOptSolversCMake_osicbc
            MathOptSolversCMake_cgl)
    endif()

    if(WIN32)
        FetchContent_Declare(cbc
            URL https://github.com/coin-or/Cbc/releases/download/releases%2F2.10.10/Cbc-releases.2.10.10-w64-msvc16-md.zip
            URL_HASH SHA256=94a951904eb80c2f351785fc9340216e88970a716634bc1ccd8b4fc5024af37c
            CONFIGURE_COMMAND ""
            BUILD_COMMAND "")
        FetchContent_MakeAvailable(cbc)

        add_library(MathOptSolversCMake_osicbc STATIC IMPORTED GLOBAL)
        target_include_directories(MathOptSolversCMake_osicbc INTERFACE
            "${osi_SOURCE_DIR}/include")
        set_target_properties(MathOptSolversCMake_osicbc PROPERTIES
            IMPORTED_LOCATION "${cbc_SOURCE_DIR}/lib/libOsiCbc.lib")

        add_library(MathOptSolversCMake_cgl STATIC IMPORTED GLOBAL)
        target_include_directories(MathOptSolversCMake_cgl INTERFACE
            "${osi_SOURCE_DIR}/include")
        set_target_properties(MathOptSolversCMake_cgl PROPERTIES
            IMPORTED_LOCATION "${cbc_SOURCE_DIR}/lib/libCgl.lib")

        target_include_directories(MathOptSolversCMake_cbc INTERFACE
            "${cbc_SOURCE_DIR}/include")
        set_target_properties(MathOptSolversCMake_cbc PROPERTIES
            IMPORTED_LOCATION "${cbc_SOURCE_DIR}/lib/libCbc.lib")
        target_link_libraries(MathOptSolversCMake_cbc INTERFACE
            MathOptSolversCMake_osicbc
            MathOptSolversCMake_cgl)
    endif()

    add_dependencies(MathOptSolversCMake_cbc cbc)
    target_link_libraries(MathOptSolversCMake_cbc INTERFACE
        MathOptSolversCMake_clp)
    add_library(MathOptSolversCMake::cbc ALIAS MathOptSolversCMake_cbc)
endif()

# Highs
if(MATHOPTSOLVERSCMAKE_USE_HIGHS)
    set(BUILD_SHARED_LIBS OFF)
    FetchContent_Declare(
        highs
        URL https://github.com/ERGO-Code/HiGHS/archive/refs/tags/v1.12.0.zip
        EXCLUDE_FROM_ALL)
    FetchContent_MakeAvailable(highs)
endif()

# FICO Xpress

if(MATHOPTSOLVERSCMAKE_USE_XPRESS)
    add_library(MathOptSolversCMake_xpress INTERFACE)
    target_include_directories(MathOptSolversCMake_xpress INTERFACE
        "$ENV{XPRESSDIR}/include")
    if(UNIX AND NOT APPLE)
        target_link_directories(MathOptSolversCMake_xpress INTERFACE
            "$ENV{XPRESSDIR}/lib/")
        target_link_libraries(MathOptSolversCMake_xpress INTERFACE
            libxprs.so)
    endif()
    if(APPLE)
        # TODO
    endif()
    if(WIN32)
        # TODO
    endif()
    add_library(MathOptSolversCMake::xpress ALIAS MathOptSolversCMake_xpress)
endif()

# IBM ILOG CPLEX Optimizer

if(MATHOPTSOLVERSCMAKE_USE_CPLEX)
    add_library(MathOptSolversCMake_concert INTERFACE)
    target_include_directories(MathOptSolversCMake_concert INTERFACE
        "$ENV{CPLEX_HOME}/concert/include")
    if(UNIX AND NOT APPLE)
        target_compile_definitions(MathOptSolversCMake_concert INTERFACE
            IL_STD=1)
        target_link_directories(MathOptSolversCMake_concert INTERFACE
            "$ENV{CPLEX_HOME}/concert/lib/x86-64_linux/static_pic/")
        target_link_libraries(MathOptSolversCMake_concert INTERFACE
            libconcert.a
            libm.a
            libpthread.a
            libdl.a)
    endif()
    if(APPLE)
        # TODO
    endif()
    if(WIN32)
        # TODO
    endif()
    add_library(MathOptSolversCMake::concert ALIAS MathOptSolversCMake_concert)

    add_library(MathOptSolversCMake_cplex INTERFACE)
    target_include_directories(MathOptSolversCMake_cplex INTERFACE
        "$ENV{CPLEX_HOME}/cplex/include")
    if(UNIX AND NOT APPLE)
        target_link_directories(MathOptSolversCMake_cplex INTERFACE
            "$ENV{CPLEX_HOME}/cplex/lib/x86-64_linux/static_pic/")
        target_link_libraries(MathOptSolversCMake_cplex INTERFACE
            libilocplex.a
            libcplex.a)
    endif()
    if(APPLE)
        # TODO
    endif()
    if(WIN32)
        # TODO
    endif()
    target_link_libraries(MathOptSolversCMake_cplex INTERFACE
        MathOptSolversCMake_concert)
    add_library(MathOptSolversCMake::cplex ALIAS MathOptSolversCMake_cplex)
endif()

# Gurobi

if(MATHOPTSOLVERSCMAKE_USE_GUROBI)
    add_library(MathOptSolversCMake_gurobi INTERFACE)
    target_include_directories(MathOptSolversCMake_gurobi INTERFACE
        "$ENV{GUROBI_HOME}/include")
    target_link_directories(MathOptSolversCMake_gurobi INTERFACE
        "$ENV{GUROBI_HOME}/lib/")
    if(UNIX AND NOT APPLE)
        target_link_libraries(MathOptSolversCMake_gurobi INTERFACE
            libgurobi_c++.a
            libgurobi_g++4.2.a)
    endif()
    if(APPLE)
        # TODO
    endif()
    if(WIN32)
        # TODO
    endif()
    add_library(MathOptSolversCMake::gurobi ALIAS MathOptSolversCMake_gurobi)
endif()

# Knitro
if(MATHOPTSOLVERSCMAKE_USE_KNITRO)
    add_library(MathOptSolversCMake_knitro INTERFACE)
    target_include_directories(MathOptSolversCMake_knitro INTERFACE
        "$ENV{KNITRODIR}/include")
    target_link_directories(MathOptSolversCMake_knitro INTERFACE
        "$ENV{KNITRODIR}/lib")
    if(UNIX AND NOT APPLE)
        target_link_libraries(MathOptSolversCMake_knitro INTERFACE
            knitro.so)
    endif()
    if(APPLE)
        target_link_libraries(MathOptSolversCMake_knitro INTERFACE
            knitro.dylib)
    endif()
    if(WIN32)
        file(GLOB KNITROLIB "$ENV{KNITRODIR}/lib/knitro*.lib")
        target_link_libraries(MathOptSolversCMake_knitro INTERFACE
            ${KNITROLIB})
    endif()
    add_library(MathOptSolversCMake::knitro ALIAS MathOptSolversCMake_knitro)
endif()

# Fetch dlib.
if(MATHOPTSOLVERSCMAKE_USE_DLIB)
    FetchContent_Declare(
        dlib
        GIT_REPOSITORY https://github.com/davisking/dlib.git
        GIT_TAG v20.0)
    FetchContent_MakeAvailable(dlib)
endif()

# Fetch ConicBundle.
if(MATHOPTSOLVERSCMAKE_USE_CONICBUNDLE)
    FetchContent_Declare(
        conicbundle
        GIT_REPOSITORY https://github.com/fontanf/ConicBundle.git
        GIT_TAG ee4f2b177be46e1750edf6093cd2fcdd033eb3c0
        EXCLUDE_FROM_ALL)
    FetchContent_MakeAvailable(conicbundle)
endif()

# Fetch SMS++ (core + BundleSolver, via the smspp-project umbrella, which
# handles the inter-module dependencies via git submodules).
#
# Unlike the other solvers above, SMS++ is not a set of prebuilt static
# libraries: it is fetched and built from source, and its CMakeLists.txt
# (a) requires CMake >= 3.21 (this repo's own cmake_minimum_required is
#     lower, so building with this option needs a newer CMake than the
#     project otherwise does), and
# (b) requires Boost (>= 1.72), Eigen3 and netCDF-cxx4 to already be
#     available on the system for the SMS++ core -- these are located with
#     find_package() and are NOT fetched here. BundleSolver also needs Osi
#     and an OsiClp (Clp) install; those ARE covered below, via the Clp
#     tarball fetched above (note the MATHOPTSOLVERSCMAKE_USE_SMSPP in that
#     if() condition).
if(MATHOPTSOLVERSCMAKE_USE_SMSPP)
    # BundleSolver's vendored FindCoinUtils/FindOsi/FindClp.cmake modules
    # (invoked by its find_package(Osi REQUIRED ...) and
    # find_package(Clp COMPONENTS OsiClp)) each have unconditional hard
    # requirements of their own (e.g. FindCoinUtils.cmake always does
    # find_package(BZip2 REQUIRED), regardless of whether CoinUtils is
    # eventually found via a *_ROOT hint) that aren't necessarily available
    # on the build machine (BZip2 isn't installed by default on Windows CI
    # runners). Redirect find_package(CoinUtils/Osi/Clp) instead (via
    # OVERRIDE_FIND_PACKAGE on the coinutils/osi/clp fetches above) to
    # targets built from the Clp tarball fetched above, which bundles
    # matching CoinUtils, Osi, Clp and OsiClp headers/libraries together
    # (safer than reusing the separately-fetched coinutils/osi targets,
    # which may be a different, potentially ABI-incompatible release).
    # Their headers are installed under "coin/", and are #include'd by
    # BundleSolver as e.g. #include "CoinPackedMatrix.hpp" (no "coin/"
    # prefix), so the include path must point directly at that subdirectory.
    foreach(_smspp_coin_lib CoinUtils Osi Clp OsiClp)
        add_library(MathOptSolversCMake_smspp_${_smspp_coin_lib} STATIC IMPORTED)
        set_target_properties(MathOptSolversCMake_smspp_${_smspp_coin_lib} PROPERTIES
            IMPORTED_LOCATION "${clp_SOURCE_DIR}/lib/lib${_smspp_coin_lib}${CMAKE_STATIC_LIBRARY_SUFFIX}")
        target_include_directories(MathOptSolversCMake_smspp_${_smspp_coin_lib} INTERFACE
            "${clp_SOURCE_DIR}/include/coin")
        add_library(Coin::${_smspp_coin_lib} ALIAS MathOptSolversCMake_smspp_${_smspp_coin_lib})
    endforeach()
    target_link_libraries(MathOptSolversCMake_smspp_Osi INTERFACE MathOptSolversCMake_smspp_CoinUtils)
    target_link_libraries(MathOptSolversCMake_smspp_Clp INTERFACE MathOptSolversCMake_smspp_CoinUtils)
    if(UNIX AND NOT APPLE)
        # CoinUtils uses BLAS/LAPACK routines (e.g. dgetrf_) that a static
        # library doesn't pull in on its own; LAPACK::LAPACK/BLAS::BLAS are
        # already found above, by the CoinUtils fetch in the COIN-OR section.
        target_link_libraries(MathOptSolversCMake_smspp_CoinUtils INTERFACE
            LAPACK::LAPACK
            BLAS::BLAS)
    endif()
    target_link_libraries(MathOptSolversCMake_smspp_OsiClp INTERFACE
        MathOptSolversCMake_smspp_Clp MathOptSolversCMake_smspp_Osi)

    # SMS++ core does find_package(Eigen3 NO_MODULE) (falling back to
    # find_package(Eigen3 REQUIRED) module mode if that fails), so fetch it
    # here instead of relying on a system install -- Windows CI runners
    # don't have one by default. Eigen is header-only, so this is cheap, and
    # its own CMakeLists.txt already defines an "Eigen3::Eigen" target for
    # exactly this build-tree (not installed) usage, so no compatibility
    # shim is needed here, unlike Boost/HiGHS/CoinUtils above.
    set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(
        eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
        EXCLUDE_FROM_ALL
        OVERRIDE_FIND_PACKAGE)
    FetchContent_MakeAvailable(eigen3)

    # BundleSolver links MILPSolver but we never call into it; disable its
    # commercial-solver backends so a locally installed but older/newer
    # CPLEX or Gurobi than what MILPSolver's code expects (mismatched
    # callback API symbols, say) can't break this unrelated build.
    set(MILPSolver_USE_CPLEX OFF CACHE BOOL "" FORCE)
    set(MILPSolver_USE_GUROBI OFF CACHE BOOL "" FORCE)

    # MILPSolver's vendored FindHiGHS.cmake unconditionally does
    # find_package(ZLIB REQUIRED), even though MILPSolver's own
    # find_package(HiGHS) call is optional -- so on a machine without ZLIB
    # (e.g. a default Windows CI runner), configuration fails hard even
    # though nothing here calls into MILPSolver's HiGHS backend (BundleSolver
    # only instantiates the plain base MILPSolver class, never
    # HiGHSMILPSolver -- its own numeric master-problem solver is CPLEX or
    # Gurobi, chosen via WHICH_OSI_QP, not MILPSolver at all).
    set(MILPSolver_USE_HiGHS OFF CACHE BOOL "" FORCE)

    # SMS++ core's own CMakeLists.txt unconditionally does
    # find_library(GETOPT_LIBRARY getopt) and links it PUBLIC on MSVC, even
    # though getopt() isn't called anywhere in SMS++'s core src/include
    # (only by its tools/tests, already disabled below) -- so on a Windows
    # CI runner without a getopt library installed, GETOPT_LIBRARY resolves
    # to "GETOPT_LIBRARY-NOTFOUND", and every consumer of SMS++ fails to
    # link against the literal file "GETOPT_LIBRARY-NOTFOUND.lib". Pre-set
    # the cache variable so that find_library() call is a no-op instead,
    # linking an empty string (itself a no-op for target_link_libraries()).
    if(MSVC)
        set(GETOPT_LIBRARY "" CACHE FILEPATH
            "Unused by SMS++'s core library; left empty to avoid a broken NOTFOUND link." FORCE)
    endif()

    if(WIN32)
        # SMS++ core's own CMakeLists.txt only tries the well-behaved,
        # Config-mode netCDFCxx lookup "if (UNIX)", so on Windows it instead
        # goes straight to its own vendored FindnetCDFCxx.cmake (Module
        # mode), which in turn uses its own vendored FindnetCDF.cmake for
        # the plain C library. Both modules' WIN32 search paths key off
        # VCPKG_INSTALLED_DIR/VCPKG_TARGET_TRIPLET, which are only set by
        # the full vcpkg toolchain file -- not used here, since it globally
        # wraps add_library()/find_package() and breaks unrelated
        # dependencies. They do also unconditionally honor a
        # netCDF(Cxx)_ROOT hint, though, and CMake's own bundled
        # FindZLIB.cmake honors a ZLIB_ROOT hint the same way, so point all
        # three at the same vcpkg-installed prefix the CI workflow already
        # passes via CMAKE_PREFIX_PATH.
        #
        # These need to be CACHE variables, not plain ones: this src/
        # directory's own later find_package(ZLIB REQUIRED) (see
        # src/CMakeLists.txt) is a *sibling* of this extern/ directory, not
        # a descendant of it, so a plain set() here wouldn't be visible
        # there -- unlike smspp_project-src, which is fetched as a
        # descendant of this same extern/ scope and so does see plain
        # variables set here.
        set(netCDF_ROOT "${CMAKE_PREFIX_PATH}" CACHE PATH "" FORCE)
        set(netCDFCxx_ROOT "${CMAKE_PREFIX_PATH}" CACHE PATH "" FORCE)

        # ZLIB_ROOT alone isn't enough: the vcpkg zlib port, built under the
        # x64-windows-static triplet, names its static library "zs.lib" --
        # not any of the names CMake's own bundled FindZLIB.cmake searches
        # for ("z", "zlib", "zlibstatic", etc.), so its NAMES-based
        # find_library() call can't find it regardless of the search path.
        # Bypass that search entirely via FindZLIB.cmake's own documented
        # manual-override mechanism ("Allow ZLIB_LIBRARY to be set
        # manually"): pre-setting ZLIB_LIBRARY skips the search, and
        # ZLIB_LIBRARY_RELEASE is what actually ends up as the ZLIB::ZLIB
        # target's IMPORTED_LOCATION_RELEASE.
        set(ZLIB_ROOT "${CMAKE_PREFIX_PATH}" CACHE PATH "" FORCE)
        set(ZLIB_LIBRARY "${CMAKE_PREFIX_PATH}/lib/zs.lib" CACHE FILEPATH "" FORCE)
        set(ZLIB_LIBRARY_RELEASE "${CMAKE_PREFIX_PATH}/lib/zs.lib" CACHE FILEPATH "" FORCE)
    endif()

    # Only build the two modules we need (SMS++ core + BundleSolver, which
    # in turn needs MILPSolver); every other Block/Solver in the umbrella
    # is switched off to avoid pulling in unrelated dependencies.
    set(BUILD_BinaryKnapsackBlock OFF CACHE BOOL "" FORCE)
    set(BUILD_UCBlock OFF CACHE BOOL "" FORCE)
    set(BUILD_MCFBlock OFF CACHE BOOL "" FORCE)
    set(BUILD_MMCFBlock OFF CACHE BOOL "" FORCE)
    set(BUILD_LukFiBlock OFF CACHE BOOL "" FORCE)
    set(BUILD_StochasticBlock OFF CACHE BOOL "" FORCE)
    set(BUILD_SDDPBlock OFF CACHE BOOL "" FORCE)
    set(BUILD_CapacitatedFacilityLocationBlock OFF CACHE BOOL "" FORCE)
    set(BUILD_InvestmentBlock OFF CACHE BOOL "" FORCE)
    set(BUILD_TwoStageStochasticBlock OFF CACHE BOOL "" FORCE)
    set(BUILD_MultiStageStochasticBlock OFF CACHE BOOL "" FORCE)
    set(BUILD_SVMBlock OFF CACHE BOOL "" FORCE)
    set(BUILD_MILPSolver ON CACHE BOOL "" FORCE)
    set(BUILD_BranchAndXSolver OFF CACHE BOOL "" FORCE)
    set(BUILD_FrankWolfeSolver OFF CACHE BOOL "" FORCE)
    set(BUILD_BundleSolver ON CACHE BOOL "" FORCE)
    set(BUILD_LagrangianDualSolver OFF CACHE BOOL "" FORCE)
    set(BUILD_MCFClassSolver OFF CACHE BOOL "" FORCE)
    set(BUILD_MCFLemonSolver OFF CACHE BOOL "" FORCE)
    set(BUILD_ScenarioReductionSolver OFF CACHE BOOL "" FORCE)
    set(BUILD_tools OFF CACHE BOOL "" FORCE)
    set(BUILD_tests OFF CACHE BOOL "" FORCE)

    # smspp-project's own CMakeSettings.txt (auto-included, see its
    # CMakeLists.txt) sets BUILD_SHARED_LIBS ON as a *plain* (non-cache)
    # variable, which shadows any value set beforehand (a CACHE ... FORCE
    # included). Appending an OFF override -- plain sets are last-write-wins
    # within that file -- is the only way to get the static build every
    # other solver here uses, needed because the CoinUtils/Osi/Clp tarballs
    # fetched above are non-PIC static archives that cannot link into a .so.
    FetchContent_Declare(
        smspp_project
        GIT_REPOSITORY https://gitlab.com/smspp/smspp-project.git
        GIT_TAG 365d6d3ace24d32e2b9037782f2ef560ffd9fe85
        GIT_SUBMODULES_RECURSE TRUE
        PATCH_COMMAND sh -c "echo 'set(BUILD_SHARED_LIBS OFF)' >> CMakeSettings.txt"
        EXCLUDE_FROM_ALL)
    FetchContent_MakeAvailable(smspp_project)
endif()

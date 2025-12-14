# Enable FetchContent.
include(FetchContent)

# COIN-OR CLP/CBC

if(MATHOPTSOLVERSCMAKE_USE_CLP OR MATHOPTSOLVERSCMAKE_USE_CBC)

    # CoinUtils

    add_library(MathOptSolversCMake_coinutils STATIC IMPORTED GLOBAL)

    if(UNIX AND NOT APPLE)
        FetchContent_Declare(coinutils
            URL https://github.com/coin-or/CoinUtils/releases/download/releases%2F2.11.9/CoinUtils-releases.2.11.9-x86_64-ubuntu20-gcc940-static.tar.gz
            URL_HASH SHA256=14d07de1b7961f68e037da6f0c57844fd67d4cc1a4b125642f42cd134b228094
            CONFIGURE_COMMAND ""
            BUILD_COMMAND "")
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
            BUILD_COMMAND "")
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
            BUILD_COMMAND "")
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
            BUILD_COMMAND "")
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
            BUILD_COMMAND "")
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
            BUILD_COMMAND "")
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
            BUILD_COMMAND "")
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
            BUILD_COMMAND "")
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
            BUILD_COMMAND "")
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
        URL https://github.com/ERGO-Code/HiGHS/archive/refs/tags/v1.11.0.zip
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
        GIT_TAG b256eda2a7ac6fe5fcc92d479f95ca606f361ad1
        EXCLUDE_FROM_ALL)
    FetchContent_MakeAvailable(conicbundle)
endif()

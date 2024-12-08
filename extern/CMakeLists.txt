# Enable ExternalProject.
include(ExternalProject)

# COIN-OR CLP/CBC

if(MATHOPTSOLVERSCMAKE_USE_CLP OR MATHOPTSOLVERSCMAKE_USE_CBC OR MATHOPTSOLVERSCMAKE_USE_ORTOOLS)
    add_library(MathOptSolversCMake_coinutils INTERFACE)
    if(UNIX AND NOT APPLE)
        ExternalProject_Add(CoinUtils
            URL https://github.com/coin-or/CoinUtils/releases/download/releases%2F2.11.9/CoinUtils-releases.2.11.9-x86_64-ubuntu20-gcc940-static.tar.gz
            URL_HASH SHA256=14d07de1b7961f68e037da6f0c57844fd67d4cc1a4b125642f42cd134b228094
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            INSTALL_COMMAND ""
            TEST_COMMAND "")
        target_include_directories(MathOptSolversCMake_coinutils INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/CoinUtils-prefix/src/CoinUtils/include")
        target_link_directories(MathOptSolversCMake_coinutils INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/CoinUtils-prefix/src/CoinUtils/lib")
        target_link_libraries(MathOptSolversCMake_coinutils INTERFACE
            libCoinUtils.a
            liblapack.a
            libblas.a
            libbz2.a
            libz.a
            libgfortran.a
            libquadmath.a)
    endif()
    if(APPLE)
        ExternalProject_Add(CoinUtils
                URL https://github.com/fontanf/coinor-builds/releases/download/v0.1/CoinUtils-release.2.11.12-universal-macos.tar.gz
                URL_HASH SHA256=8cd0b09e2eed185b5b37f81ef66554e0d41f29e9b5ea87b128538ce158fa6bfa
                CONFIGURE_COMMAND ""
                BUILD_COMMAND ""
                INSTALL_COMMAND ""
                TEST_COMMAND "")
        target_include_directories(MathOptSolversCMake_coinutils INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/CoinUtils-prefix/src/CoinUtils/include")
        target_link_directories(MathOptSolversCMake_coinutils INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/CoinUtils-prefix/src/CoinUtils/lib")
        target_link_libraries(MathOptSolversCMake_coinutils INTERFACE
                libCoinUtils.a)
    endif()
    if(WIN32)
        ExternalProject_Add(CoinUtils
            URL https://github.com/coin-or/CoinUtils/releases/download/releases%2F2.11.9/CoinUtils-releases.2.11.9-w64-msvc16-md.zip
            URL_HASH SHA256=2bc64f0afd80113571697e949b2663e9047272decf90d5f62e452c2628d33ca6
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            INSTALL_COMMAND ""
            TEST_COMMAND "")
        target_include_directories(MathOptSolversCMake_coinutils INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/CoinUtils-prefix/src/CoinUtils/include")
        target_link_directories(MathOptSolversCMake_coinutils INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/CoinUtils-prefix/src/CoinUtils/lib")
        target_link_libraries(MathOptSolversCMake_coinutils INTERFACE
            libCoinUtils.lib)
    endif()
    add_dependencies(MathOptSolversCMake_coinutils CoinUtils)
    add_library(MathOptSolversCMake::coinutils ALIAS MathOptSolversCMake_coinutils)

    add_library(MathOptSolversCMake_osi INTERFACE)
    if(UNIX AND NOT APPLE)
        ExternalProject_Add(Osi
            URL https://github.com/coin-or/Osi/releases/download/releases%2F0.108.8/Osi-releases.0.108.8-x86_64-ubuntu20-gcc940-static.tar.gz
            URL_HASH SHA256=bd5a5bf1e6b6a28d13d41ab1554becd9f3992afe775785e51a88c9405cf2853e
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            INSTALL_COMMAND ""
            TEST_COMMAND "")
        target_include_directories(MathOptSolversCMake_osi INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/Osi-prefix/src/Osi/include")
        target_link_directories(MathOptSolversCMake_osi INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/Osi-prefix/src/Osi/lib")
        target_link_libraries(MathOptSolversCMake_osi INTERFACE
            libOsi.a
            libOsiCommonTests.a)
    endif()
    if(APPLE)
        ExternalProject_Add(Osi
                URL https://github.com/fontanf/coinor-builds/releases/download/v0.1/Osi-release.0.108.11-universal-macos.tar.gz
                URL_HASH SHA256=fc6b7b93106e0d4d0d8eea1b2fb9f9db3c345580444409655133fc984391682b
                CONFIGURE_COMMAND ""
                BUILD_COMMAND ""
                INSTALL_COMMAND ""
                TEST_COMMAND "")
        target_include_directories(MathOptSolversCMake_osi INTERFACE
                "${CMAKE_CURRENT_BINARY_DIR}/Osi-prefix/src/Osi/include")
        target_link_directories(MathOptSolversCMake_osi INTERFACE
                "${CMAKE_CURRENT_BINARY_DIR}/Osi-prefix/src/Osi/lib")
        target_link_libraries(MathOptSolversCMake_osi INTERFACE
                libOsi.a
                libOsiCommonTests.a)
    endif()
    if(WIN32)
        ExternalProject_Add(Osi
            URL https://github.com/coin-or/Osi/releases/download/releases%2F0.108.8/Osi-releases.0.108.8-w64-msvc16-md.zip
            URL_HASH SHA256=a61fc462cb598139d205cd2323522581a01900575d0d6bccf660b6c7e1b0b71c
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            INSTALL_COMMAND ""
            TEST_COMMAND "")
        target_include_directories(MathOptSolversCMake_osi INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/Osi-prefix/src/Osi/include")
        target_link_directories(MathOptSolversCMake_osi INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/Osi-prefix/src/Osi/lib")
        target_link_libraries(MathOptSolversCMake_osi INTERFACE
            libOsi.lib
            libOsiCommonTests.lib)
    endif()
    add_dependencies(MathOptSolversCMake_osi Osi)
    target_link_libraries(MathOptSolversCMake_osi INTERFACE
        MathOptSolversCMake_coinutils)
    add_library(MathOptSolversCMake::osi ALIAS MathOptSolversCMake_osi)

    add_library(MathOptSolversCMake_clp INTERFACE)
    if(UNIX AND NOT APPLE)
        ExternalProject_Add(Clp
            URL https://github.com/coin-or/Clp/releases/download/releases%2F1.17.8/Clp-releases.1.17.8-x86_64-ubuntu20-gcc940-static.tar.gz
            URL_HASH SHA256=d569b04d19c25876e55d2557a1d9739df8eb50ec8ca11a98ad387fd8b90212c9
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            INSTALL_COMMAND ""
            TEST_COMMAND "")
        target_include_directories(MathOptSolversCMake_clp INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/Clp-prefix/src/Clp/include")
        target_link_directories(MathOptSolversCMake_clp INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/Clp-prefix/src/Clp/lib")
        target_link_libraries(MathOptSolversCMake_clp INTERFACE
            libClp.a
            libOsiClp.a)
    endif()
    if(APPLE)
        ExternalProject_Add(Clp
            URL https://github.com/fontanf/coinor-builds/releases/download/v0.1/Clp-release.1.17.10-universal-macos.tar.gz
            URL_HASH SHA256=1179ae58b6d0fe827ff8933f72eed2a674273968d12e8125a93f5b03f9e8757b
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            INSTALL_COMMAND ""
            TEST_COMMAND "")
        target_include_directories(MathOptSolversCMake_clp INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/Clp-prefix/src/Clp/include")
        target_link_directories(MathOptSolversCMake_clp INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/Clp-prefix/src/Clp/lib")
        target_link_libraries(MathOptSolversCMake_clp INTERFACE
            libClp.a
            libOsiClp.a)
    endif()
    if(WIN32)
        ExternalProject_Add(Clp
            URL https://github.com/coin-or/Clp/releases/download/releases%2F1.17.8/Clp-releases.1.17.8-w64-msvc16-md.zip
            URL_HASH SHA256=e37c834aea5c31dfd8620b7d2432cb31fc16ecb0c6ffb398e8f07c9c82dd5028
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            INSTALL_COMMAND ""
            TEST_COMMAND "")
        target_include_directories(MathOptSolversCMake_clp INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/Clp-prefix/src/Clp/include")
        target_link_directories(MathOptSolversCMake_clp INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/Clp-prefix/src/Clp/lib")
        target_link_libraries(MathOptSolversCMake_clp INTERFACE
            libClp.lib
            libOsiClp.lib)
    endif()
    add_dependencies(MathOptSolversCMake_clp Clp)
    target_link_libraries(MathOptSolversCMake_clp INTERFACE
        MathOptSolversCMake_osi)
    add_library(MathOptSolversCMake::clp ALIAS MathOptSolversCMake_clp)
endif()

if(MATHOPTSOLVERSCMAKE_USE_CBC OR MATHOPTSOLVERSCMAKE_USE_ORTOOLS)
    add_library(MathOptSolversCMake_cbc INTERFACE)
    if(UNIX AND NOT APPLE)
        ExternalProject_Add(Cbc
            URL https://github.com/coin-or/Cbc/releases/download/releases%2F2.10.10/Cbc-releases.2.10.10-x86_64-ubuntu20-gcc940-static.tar.gz
            URL_HASH SHA256=872c78bfcdd1566f134d2f7757b76b2a2479a5b1ade065cdd1d4b303ed6f8006
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            INSTALL_COMMAND ""
            TEST_COMMAND "")
        target_include_directories(MathOptSolversCMake_cbc INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/Cbc-prefix/src/Cbc/include")
        target_link_directories(MathOptSolversCMake_cbc INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/Cbc-prefix/src/Cbc/lib")
        target_link_libraries(MathOptSolversCMake_cbc INTERFACE
            libCbc.a
            libOsiCbc.a
            libCgl.a
            libnauty.a)
    endif()
    if(APPLE)
        ExternalProject_Add(Cbc
            URL https://github.com/fontanf/coinor-builds/releases/download/v0.1/Cbc-release.2.10.12-universal-macos.tar.gz
            URL_HASH SHA256=501199f4d14f33555dedd6d233d7c51a097088edcef15b74d8304469e32ecade
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            INSTALL_COMMAND ""
            TEST_COMMAND "")
        target_include_directories(MathOptSolversCMake_cbc INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/Cbc-prefix/src/Cbc/include")
        target_link_directories(MathOptSolversCMake_cbc INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/Cbc-prefix/src/Cbc/lib")
        target_link_libraries(MathOptSolversCMake_cbc INTERFACE
            libCbc.a
            libOsiCbc.a
            libCgl.a)
    endif()
    if(WIN32)
        ExternalProject_Add(Cbc
            URL https://github.com/coin-or/Cbc/releases/download/releases%2F2.10.10/Cbc-releases.2.10.10-w64-msvc16-md.zip
            URL_HASH SHA256=94a951904eb80c2f351785fc9340216e88970a716634bc1ccd8b4fc5024af37c
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            INSTALL_COMMAND ""
            TEST_COMMAND "")
        target_include_directories(MathOptSolversCMake_cbc INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/Cbc-prefix/src/Cbc/include")
        target_link_directories(MathOptSolversCMake_cbc INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/Cbc-prefix/src/Cbc/lib")
        target_link_libraries(MathOptSolversCMake_cbc INTERFACE
            libCbc.lib
            libOsiCbc.lib
            libCgl.lib)
    endif()
    add_dependencies(MathOptSolversCMake_cbc Cbc)
    target_link_libraries(MathOptSolversCMake_cbc INTERFACE
        MathOptSolversCMake_clp)
    add_library(MathOptSolversCMake::cbc ALIAS MathOptSolversCMake_cbc)
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

# OR-Tools
if(MATHOPTSOLVERSCMAKE_USE_ORTOOLS)
    add_library(MathOptSolversCMake_ortools INTERFACE)
    if(UNIX AND NOT APPLE)
        ExternalProject_Add(ORTools
            URL https://github.com/google/or-tools/releases/download/v9.11/or-tools_amd64_ubuntu-24.04_cpp_v9.11.4210.tar.gz
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            INSTALL_COMMAND ""
            TEST_COMMAND "")
        target_include_directories(MathOptSolversCMake_ortools INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/ORTools-prefix/src/ORTools/include")
        target_link_directories(MathOptSolversCMake_ortools INTERFACE
            "${CMAKE_CURRENT_BINARY_DIR}/ORTools-prefix/src/ORTools/lib")
        target_link_libraries(MathOptSolversCMake_ortools INTERFACE
            libabsl_bad_any_cast_impl.a
            libabsl_bad_optional_access.a
            libabsl_bad_variant_access.a
            libabsl_base.a
            libabsl_city.a
            libabsl_civil_time.a
            libabsl_cord.a
            libabsl_cord_internal.a
            libabsl_cordz_functions.a
            libabsl_cordz_handle.a
            libabsl_cordz_info.a
            libabsl_cordz_sample_token.a
            libabsl_crc32c.a
            libabsl_crc_cord_state.a
            libabsl_crc_cpu_detect.a
            libabsl_crc_internal.a
            libabsl_debugging_internal.a
            libabsl_decode_rust_punycode.a
            libabsl_demangle_internal.a
            libabsl_demangle_rust.a
            libabsl_die_if_null.a
            libabsl_examine_stack.a
            libabsl_exponential_biased.a
            libabsl_failure_signal_handler.a
            libabsl_flags_commandlineflag.a
            libabsl_flags_commandlineflag_internal.a
            libabsl_flags_config.a
            libabsl_flags_internal.a
            libabsl_flags_marshalling.a
            libabsl_flags_parse.a
            libabsl_flags_private_handle_accessor.a
            libabsl_flags_program_name.a
            libabsl_flags_reflection.a
            libabsl_flags_usage.a
            libabsl_flags_usage_internal.a
            libabsl_graphcycles_internal.a
            libabsl_hash.a
            libabsl_hashtablez_sampler.a
            libabsl_int128.a
            libabsl_kernel_timeout_internal.a
            libabsl_leak_check.a
            libabsl_log_entry.a
            libabsl_log_flags.a
            libabsl_log_globals.a
            libabsl_log_initialize.a
            libabsl_log_internal_check_op.a
            libabsl_log_internal_conditions.a
            libabsl_log_internal_fnmatch.a
            libabsl_log_internal_format.a
            libabsl_log_internal_globals.a
            libabsl_log_internal_log_sink_set.a
            libabsl_log_internal_message.a
            libabsl_log_internal_nullguard.a
            libabsl_log_internal_proto.a
            libabsl_log_severity.a
            libabsl_log_sink.a
            libabsl_low_level_hash.a
            libabsl_malloc_internal.a
            libabsl_periodic_sampler.a
            libabsl_poison.a
            libabsl_random_distributions.a
            libabsl_random_internal_distribution_test_util.a
            libabsl_random_internal_platform.a
            libabsl_random_internal_pool_urbg.a
            libabsl_random_internal_randen.a
            libabsl_random_internal_randen_hwaes.a
            libabsl_random_internal_randen_hwaes_impl.a
            libabsl_random_internal_randen_slow.a
            libabsl_random_internal_seed_material.a
            libabsl_random_seed_gen_exception.a
            libabsl_random_seed_sequences.a
            libabsl_raw_hash_set.a
            libabsl_raw_logging_internal.a
            libabsl_scoped_set_env.a
            libabsl_spinlock_wait.a
            libabsl_stacktrace.a
            libabsl_status.a
            libabsl_statusor.a
            libabsl_strerror.a
            libabsl_str_format_internal.a
            libabsl_strings.a
            libabsl_strings_internal.a
            libabsl_string_view.a
            libabsl_symbolize.a
            libabsl_synchronization.a
            libabsl_throw_delegate.a
            libabsl_time.a
            libabsl_time_zone.a
            libabsl_utf8_for_code_point.a
            libabsl_vlog_config_internal.a
            #libCbc.a
            #libCbcSolver.a
            #libCgl.a
            #libClp.a
            #libClpSolver.a
            #libCoinUtils.a
            #libOsi.a
            #libOsiCbc.a
            #libOsiClp.a
            libprotobuf.a
            libprotobuf-lite.a
            libprotoc.a
            libre2.a
            libscip.a
            libutf8_range.a
            libutf8_validity.a
            libz.a)
    endif()
    add_dependencies(MathOptSolversCMake_ortools ORTools)
    target_link_libraries(MathOptSolversCMake_ortools INTERFACE
        MathOptSolversCMake_cbc)
    add_library(MathOptSolversCMake::ortools ALIAS MathOptSolversCMake_ortools)
endif()

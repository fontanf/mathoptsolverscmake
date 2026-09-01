add_library(MathOptSolversCMake_mathopt)
target_sources(MathOptSolversCMake_mathopt PRIVATE
    mathopt.cpp
    mathopt_nonlinear.cpp)
if(MATHOPTSOLVERSCMAKE_USE_CBC)
    target_sources(MathOptSolversCMake_mathopt PRIVATE mathopt_cbc.cpp)
endif()
if(MATHOPTSOLVERSCMAKE_USE_HIGHS)
    target_sources(MathOptSolversCMake_mathopt PRIVATE mathopt_highs.cpp)
endif()
if(MATHOPTSOLVERSCMAKE_USE_XPRESS)
    target_sources(MathOptSolversCMake_mathopt PRIVATE mathopt_xpress.cpp)
endif()
if(MATHOPTSOLVERSCMAKE_USE_KNITRO)
    target_sources(MathOptSolversCMake_mathopt PRIVATE mathopt_knitro.cpp)
endif()
if(MATHOPTSOLVERSCMAKE_USE_DLIB)
    target_sources(MathOptSolversCMake_mathopt PRIVATE mathopt_dlib.cpp)
endif()
if(MATHOPTSOLVERSCMAKE_USE_CONICBUNDLE)
    target_sources(MathOptSolversCMake_mathopt PRIVATE mathopt_conicbundle.cpp)
endif()
if(MATHOPTSOLVERSCMAKE_USE_SMSPP)
    target_sources(MathOptSolversCMake_mathopt PRIVATE mathopt_smspp.cpp)
endif()
target_include_directories(MathOptSolversCMake_mathopt PUBLIC
    ${PROJECT_SOURCE_DIR}/include)
if(MATHOPTSOLVERSCMAKE_USE_CBC)
    target_link_libraries(MathOptSolversCMake_mathopt PUBLIC
        MathOptSolversCMake_cbc)
    target_compile_definitions(MathOptSolversCMake_mathopt PUBLIC
        CBC_FOUND=1)
endif()
if(MATHOPTSOLVERSCMAKE_USE_HIGHS)
    target_link_libraries(MathOptSolversCMake_mathopt PUBLIC
        highs)
    target_compile_definitions(MathOptSolversCMake_mathopt PUBLIC
        HIGHS_FOUND=1)
endif()
if(MATHOPTSOLVERSCMAKE_USE_XPRESS)
    target_link_libraries(MathOptSolversCMake_mathopt PUBLIC
        MathOptSolversCMake_xpress)
    target_compile_definitions(MathOptSolversCMake_mathopt PUBLIC
        XPRESS_FOUND=1)
endif()
if(MATHOPTSOLVERSCMAKE_USE_DLIB)
    target_link_libraries(MathOptSolversCMake_mathopt PUBLIC
        dlib::dlib)
    target_compile_definitions(MathOptSolversCMake_mathopt PUBLIC
        DLIB_FOUND=1)
endif()
if(MATHOPTSOLVERSCMAKE_USE_KNITRO)
    target_link_libraries(MathOptSolversCMake_mathopt PUBLIC
        MathOptSolversCMake_knitro)
    target_compile_definitions(MathOptSolversCMake_mathopt PUBLIC
        KNITRO_FOUND=1)
endif()
if(MATHOPTSOLVERSCMAKE_USE_CONICBUNDLE)
    target_link_libraries(MathOptSolversCMake_mathopt PUBLIC
        cb)
    target_compile_definitions(MathOptSolversCMake_mathopt PUBLIC
        CONICBUNDLE_FOUND=1)
endif()
if(MATHOPTSOLVERSCMAKE_USE_SMSPP)
    # BundleSolver's own Coin::CoinUtils/Coin::Clp imported targets (built
    # from the CoinUtils/Clp tarball fetched in extern/CMakeLists.txt) don't
    # themselves declare their ZLIB/BZip2 link dependency, so it has to be
    # added explicitly here too, or the final link fails on the static
    # archive's unresolved gzread. The ZLIB::ZLIB/BZip2::BZip2 targets
    # extern/CMakeLists.txt finds aren't visible here (directory-scoped,
    # and extern/ and src/ are sibling subdirectories), hence re-finding.
    find_package(ZLIB REQUIRED)
    find_package(BZip2 REQUIRED)
    target_link_libraries(MathOptSolversCMake_mathopt PUBLIC
        SMS++::SMS++
        SMS++::BundleSolver
        ZLIB::ZLIB
        BZip2::BZip2)
    target_compile_definitions(MathOptSolversCMake_mathopt PUBLIC
        SMSPP_FOUND=1)
endif()
add_library(MathOptSolversCMake::mathopt ALIAS MathOptSolversCMake_mathopt)

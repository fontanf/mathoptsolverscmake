add_library(MathOptSolversCMake_mathopt)
target_sources(MathOptSolversCMake_mathopt PRIVATE
    mathopt.cpp)
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
        KnitroCpp::knitrocpp)
    target_compile_definitions(MathOptSolversCMake_mathopt PUBLIC
        KNITRO_FOUND=1)
endif()
if(MATHOPTSOLVERSCMAKE_USE_CONICBUNDLE)
    target_link_libraries(MathOptSolversCMake_mathopt PUBLIC
        cb)
    target_compile_definitions(MathOptSolversCMake_mathopt PUBLIC
        CONICBUNDLE_FOUND=1)
endif()
add_library(MathOptSolversCMake::mathopt ALIAS MathOptSolversCMake_mathopt)

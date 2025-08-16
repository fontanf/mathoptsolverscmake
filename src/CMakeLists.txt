add_library(MathOptSolversCMake_milp)
target_sources(MathOptSolversCMake_milp PRIVATE
    milp.cpp)
target_include_directories(MathOptSolversCMake_milp PUBLIC
    ${PROJECT_SOURCE_DIR}/include)
if(MATHOPTSOLVERSCMAKE_USE_CBC)
    target_link_libraries(MathOptSolversCMake_milp PUBLIC
        MathOptSolversCMake_cbc)
    target_compile_definitions(MathOptSolversCMake_milp PUBLIC
        CBC_FOUND=1)
endif()
if(MATHOPTSOLVERSCMAKE_USE_HIGHS)
    target_link_libraries(MathOptSolversCMake_milp PUBLIC
        highs)
    target_compile_definitions(MathOptSolversCMake_milp PUBLIC
        HIGHS_FOUND=1)
endif()
if(MATHOPTSOLVERSCMAKE_USE_XPRESS)
    target_link_libraries(MathOptSolversCMake_milp PUBLIC
        MathOptSolversCMake_xpress)
    target_compile_definitions(MathOptSolversCMake_milp PUBLIC
        XPRESS_FOUND=1)
endif()
add_library(MathOptSolversCMake::milp ALIAS MathOptSolversCMake_milp)

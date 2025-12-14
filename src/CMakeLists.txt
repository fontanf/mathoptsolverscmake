add_library(MathOptSolversCMake_common)
target_sources(MathOptSolversCMake_common PRIVATE
    common.cpp)
target_include_directories(MathOptSolversCMake_common PUBLIC
    ${PROJECT_SOURCE_DIR}/include)
add_library(MathOptSolversCMake::common ALIAS MathOptSolversCMake_common)

add_library(MathOptSolversCMake_milp)
target_sources(MathOptSolversCMake_milp PRIVATE
    milp.cpp)
target_include_directories(MathOptSolversCMake_milp PUBLIC
    ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(MathOptSolversCMake_milp PUBLIC
    MathOptSolversCMake_common)
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

add_library(MathOptSolversCMake_box_constrained_nlp)
target_sources(MathOptSolversCMake_box_constrained_nlp PRIVATE
    box_constrained_nlp.cpp)
target_include_directories(MathOptSolversCMake_box_constrained_nlp PUBLIC
    ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(MathOptSolversCMake_box_constrained_nlp PUBLIC
    MathOptSolversCMake_common)
if(MATHOPTSOLVERSCMAKE_USE_DLIB)
    target_link_libraries(MathOptSolversCMake_box_constrained_nlp PUBLIC
        dlib::dlib)
    target_compile_definitions(MathOptSolversCMake_box_constrained_nlp PUBLIC
        DLIB_FOUND=1)
endif()
if(MATHOPTSOLVERSCMAKE_USE_CONICBUNDLE)
    target_link_libraries(MathOptSolversCMake_box_constrained_nlp PUBLIC
        cb)
    target_compile_definitions(MathOptSolversCMake_box_constrained_nlp PUBLIC
        CONICBUNDLE_FOUND=1)
endif()
add_library(MathOptSolversCMake::box_constrained_nlp ALIAS MathOptSolversCMake_box_constrained_nlp)

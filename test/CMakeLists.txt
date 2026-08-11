# GoogleTest requires a newer standard than the library's public C++11; this
# override is scoped to the test/ subdirectory only.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.15.2)
FetchContent_MakeAvailable(googletest)

add_library(MathOptSolversCMake_test_models STATIC
    test_models.cpp)
target_include_directories(MathOptSolversCMake_test_models PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(MathOptSolversCMake_test_models PUBLIC
    MathOptSolversCMake::mathopt)

if(MATHOPTSOLVERSCMAKE_USE_KNITRO)
    add_executable(mathopt_knitro_test
        mathopt_knitro_test.cpp)
    target_link_libraries(mathopt_knitro_test PRIVATE
        MathOptSolversCMake_test_models
        MathOptSolversCMake::mathopt
        GTest::gtest_main)

    include(GoogleTest)
    gtest_discover_tests(mathopt_knitro_test)
endif()

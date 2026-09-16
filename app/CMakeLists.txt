# boost::program_options requires a newer standard than the library's public
# C++11; this override is scoped to the app/ subdirectory only.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(mathoptsolverscmake_solve
    solve.cpp)
target_link_libraries(mathoptsolverscmake_solve PRIVATE
    MathOptSolversCMake::mathopt
    Boost::program_options)

install(TARGETS mathoptsolverscmake_solve)

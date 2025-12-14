#pragma once

#include <istream>

namespace mathoptsolverscmake
{

enum class SolverName
{
    Cbc,
    Highs,
    Xpress,
    Dlib,
    ConicBundle,
};

std::istream& operator>>(
        std::istream& in,
        SolverName& solver_name);

std::ostream& operator<<(
        std::ostream& os,
        SolverName solver_name);

enum class ObjectiveDirection
{
    Minimize,
    Maximize,
};

std::istream& operator>>(
        std::istream& in,
        ObjectiveDirection& objective_direction);

std::ostream& operator<<(
        std::ostream& os,
        ObjectiveDirection objective_direction);

}

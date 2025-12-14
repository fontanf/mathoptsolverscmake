#include "mathoptsolverscmake/common.hpp"

#include <iostream>

using namespace mathoptsolverscmake;

std::istream& mathoptsolverscmake::operator>>(
        std::istream& in,
        SolverName& solver_name)
{
    std::string token;
    in >> token;
    if (token == "cbc"
            || token == "Cbc"
            || token == "CBC") {
        solver_name = SolverName::Cbc;
    } else if (token == "highs"
            || token == "Highs"
            || token == "HiGHS"
            || token == "HIGHS") {
        solver_name = SolverName::Highs;
    } else if (token == "xpress"
            || token == "Xpress"
            || token == "XPRESS") {
        solver_name = SolverName::Xpress;
    } else if (token == "dlib") {
        solver_name = SolverName::Dlib;
    } else  {
        in.setstate(std::ios_base::failbit);
    }
    return in;
}

std::ostream& mathoptsolverscmake::operator<<(
        std::ostream& os,
        SolverName solver_name)
{
    switch (solver_name) {
    case SolverName::Cbc: {
        os << "Cbc";
        break;
    } case SolverName::Highs: {
        os << "HiGHS";
        break;
    } case SolverName::Xpress: {
        os << "XPRESS";
        break;
    } case SolverName::Dlib: {
        os << "dlib";
        break;
    }
    }
    return os;
}

std::istream& mathoptsolverscmake::operator>>(
        std::istream& in,
        ObjectiveDirection& objective_direction)
{
    std::string token;
    std::getline(in, token);
    if (token == "min"
            || token == "Min"
            || token == "minimize"
            || token == "Minimize") {
        objective_direction = ObjectiveDirection::Minimize;
    } else if (token == "max"
            || token == "Max"
            || token == "maximize"
            || token == "Maximize") {
        objective_direction = ObjectiveDirection::Maximize;
    } else  {
        //throw std::invalid_argument(
        //        FUNC_SIGNATURE + ": "
        //        "invalid input; "
        //        "in: " + token + ".");
        in.setstate(std::ios_base::failbit);
    }
    return in;
}

std::ostream& mathoptsolverscmake::operator<<(
        std::ostream& os,
        ObjectiveDirection objective_direction)
{
    switch (objective_direction) {
    case ObjectiveDirection::Minimize: {
        os << "Minimize";
        break;
    } case ObjectiveDirection::Maximize: {
        os << "Maximize";
        break;
    }
    }
    return os;
}

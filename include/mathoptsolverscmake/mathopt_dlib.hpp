#pragma once

#include "mathoptsolverscmake/mathopt.hpp"

#include <vector>

namespace mathoptsolverscmake
{

struct DlibOutput
{
    double objective_value = 0;

    std::vector<double> solution;
};

DlibOutput solve_dlib(
        const MathOptModel& model);

}

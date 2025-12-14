#pragma once

#include "mathoptsolverscmake/common.hpp"

#include <vector>
#include <limits>
#include <functional>

namespace mathoptsolverscmake
{

struct BoxConstrainedNlpFunctionOutput
{
    double objective_value = 0;

    std::vector<double> gradient;
};

using BoxConstrainedNlpFunction = std::function<BoxConstrainedNlpFunctionOutput(const std::vector<double>&)>;

struct BoxConstrainedNlpModel
{
    BoxConstrainedNlpModel(
            int number_of_variables = 0):
        variables_lower_bounds(number_of_variables, -std::numeric_limits<double>::infinity()),
        variables_upper_bounds(number_of_variables, +std::numeric_limits<double>::infinity())
    { }

    /** Get the number of variables. */
    int number_of_variables() const { return variables_lower_bounds.size(); }


    ObjectiveDirection objective_direction;

    std::vector<double> variables_initial_values;
    std::vector<double> variables_lower_bounds;
    std::vector<double> variables_upper_bounds;

    BoxConstrainedNlpFunction objective_function = [](const std::vector<double>&) { return BoxConstrainedNlpFunctionOutput(); };
};

#ifdef DLIB_FOUND

struct BoxConstrainedNlpDlibOutput
{
    double objective_value = 0;

    std::vector<double> solution;
};

BoxConstrainedNlpDlibOutput solve_dlib(
        BoxConstrainedNlpModel& dlib_model);

#endif

}

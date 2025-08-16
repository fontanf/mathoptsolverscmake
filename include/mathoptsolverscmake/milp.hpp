#pragma once

#include <vector>

#ifdef HIGHS_FOUND
#include "Highs.h"
#endif

#ifdef CBC_FOUND
#include "coin/OsiCbcSolverInterface.hpp"
#endif

#ifdef XPRESS_FOUND
#include "xprs.h"
#endif

namespace mathoptsolverscmake
{

enum class ObjectiveDirection
{
    Minimize,
    Maximize,
};

enum class VariableType
{
    Continuous,
    Binary,
    Integer,
};

struct MilpModel
{
    MilpModel(
            int number_of_variables = 0,
            int number_of_constraints = 0,
            int number_of_elements = 0):
        variables_lower_bounds(number_of_variables, -std::numeric_limits<double>::infinity()),
        variables_upper_bounds(number_of_variables, +std::numeric_limits<double>::infinity()),
        variables_types(number_of_variables),
        objective_coefficients(number_of_variables),
        constraints_lower_bounds(number_of_constraints, -std::numeric_limits<double>::infinity()),
        constraints_upper_bounds(number_of_constraints, +std::numeric_limits<double>::infinity()),
        constraints_starts(number_of_constraints),
        elements_variables(number_of_elements),
        elements_coefficients(number_of_elements)
    { }

    int number_of_variables() const { return variables_lower_bounds.size(); }
    int number_of_constraints() const { return constraints_lower_bounds.size(); }
    int number_of_elements() const { return elements_variables.size(); }

    ObjectiveDirection objective_direction;

    std::vector<double> variables_lower_bounds;
    std::vector<double> variables_upper_bounds;
    std::vector<VariableType> variables_types;

    std::vector<double> objective_coefficients;

    std::vector<double> constraints_lower_bounds;
    std::vector<double> constraints_upper_bounds;
    std::vector<int> constraints_starts;
    std::vector<int> elements_variables;
    std::vector<double> elements_coefficients;
};

#ifdef CBC_FOUND

namespace cbc
{

void load(
        OsiCbcSolverInterface& cbc_model,
        const MilpModel& model);

}

#endif

#ifdef HIGHS_FOUND

namespace highs
{

void load(
        Highs& highs_model,
        const MilpModel& model);

}

#endif

#ifdef XPRESS_FOUND

namespace xpress
{

void load(
        XPRSprob& xpress_model,
        const MilpModel& model);

}

#endif

}

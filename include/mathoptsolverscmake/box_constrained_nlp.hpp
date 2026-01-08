#pragma once

#include "mathoptsolverscmake/common.hpp"

#ifdef KNITRO_FOUND
#include "knitrocpp/knitro.hpp"
#endif

#ifdef CONICBUNDLE_FOUND
#include "CBSolver.hxx"
#endif

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

#ifdef KNITRO_FOUND

void solve(
        const BoxConstrainedNlpModel& model,
        knitrocpp::Context& knitro_context);

double get_solution_value(
        const knitrocpp::Context& knitro_context);

std::vector<double> get_solution(
        const knitrocpp::Context& knitro_context);

#endif

#ifdef DLIB_FOUND

struct BoxConstrainedNlpDlibOutput
{
    double objective_value = 0;

    std::vector<double> solution;
};

BoxConstrainedNlpDlibOutput solve_dlib(
        const BoxConstrainedNlpModel& model);

#endif

#ifdef CONICBUNDLE_FOUND

void solve(
        const BoxConstrainedNlpModel& model,
        ConicBundle::CBSolver& solver);

double get_solution_value(
        const BoxConstrainedNlpModel& model,
        const ConicBundle::CBSolver& solver);

std::vector<double> get_solution(
        const BoxConstrainedNlpModel& model,
        const ConicBundle::CBSolver& solver);

#endif

}

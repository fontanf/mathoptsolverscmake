#include "mathoptsolverscmake/box_constrained_nlp.hpp"

#ifdef DLIB_FOUND
#include "dlib/optimization.h"
#endif

using namespace mathoptsolverscmake;

#ifdef DLIB_FOUND

typedef dlib::matrix<double,0,1> column_vector;

class BoxConstrainedNlpDlibFunction
{

public:

    BoxConstrainedNlpDlibFunction(
            const BoxConstrainedNlpModel& model):
        model_(model),
        grad_(model.number_of_variables()) { }

    double f(const column_vector& x)
    {
        std::vector<double> variables(model_.number_of_variables());
        for (int variable_id = 0;
                variable_id < model_.number_of_variables();
                ++variable_id) {
            variables[variable_id] = x(variable_id);
        }
        BoxConstrainedNlpFunctionOutput output = model_.objective_function(variables);

        for (int variable_id = 0;
                variable_id < model_.number_of_variables();
                ++variable_id) {
            grad_(variable_id) = output.gradient[variable_id];
        }
        return output.objective_value;
    }

    const column_vector der(const column_vector& x) const { (void)x; return grad_; }

private:

    const BoxConstrainedNlpModel& model_;

    column_vector grad_;

};

BoxConstrainedNlpDlibOutput mathoptsolverscmake::solve_dlib(
        BoxConstrainedNlpModel& model)
{
    BoxConstrainedNlpDlibOutput output;

    column_vector dlib_variables(model.number_of_variables());
    if (!model.variables_initial_values.empty()) {
        for (int variable_id = 0;
                variable_id < model.number_of_variables();
                ++variable_id) {
            dlib_variables(variable_id) = model.variables_initial_values[variable_id];
        }
    }

    column_vector dlib_variables_lower_bounds(model.number_of_variables());
    column_vector dlib_variables_upper_bounds(model.number_of_variables());
    for (int variable_id = 0;
            variable_id < model.number_of_variables();
            ++variable_id) {
        dlib_variables_lower_bounds(variable_id) = model.variables_lower_bounds[variable_id];
        dlib_variables_upper_bounds(variable_id) = model.variables_upper_bounds[variable_id];
    }

    // Solve
    BoxConstrainedNlpDlibFunction dlib_function(model);
    auto f = [&dlib_function](const column_vector& x) { return dlib_function.f(x); };
    auto def = [&dlib_function](const column_vector& x) { return dlib_function.der(x); };
    //auto stop_strategy = dlib::objective_delta_stop_strategy(0.0001);
    auto stop_strategy = dlib::objective_delta_stop_strategy(0.0001).be_verbose();
    output.objective_value = (model.objective_direction == ObjectiveDirection::Minimize)?
        dlib::find_min_box_constrained(
            dlib::lbfgs_search_strategy(256),
            stop_strategy,
            f,
            def,
            dlib_variables,
            dlib_variables_lower_bounds,
            dlib_variables_upper_bounds):
        dlib::find_max_box_constrained(
            dlib::lbfgs_search_strategy(256),
            stop_strategy,
            f,
            def,
            dlib_variables,
            dlib_variables_lower_bounds,
            dlib_variables_upper_bounds);

    // Compute output parameters
    output.solution = std::vector<double>(model.number_of_variables());
    for (int variable_id = 0;
            variable_id < model.number_of_variables();
            ++variable_id) {
        output.solution[variable_id] = dlib_variables(variable_id);
    }
    return output;
}

#endif

#include "mathoptsolverscmake/mathopt_statistics.hpp"

#include "mathoptsolverscmake/mathopt_json.hpp"
#include "mathoptsolverscmake/mathopt_lp.hpp"
#include "mathoptsolverscmake/mathopt_mps.hpp"
#include "mathoptsolverscmake/mathopt_nl.hpp"

#include <cctype>
#include <stdexcept>

using namespace mathoptsolverscmake;

MathOptModelStatistics mathoptsolverscmake::compute_statistics(const MathOptModel& model)
{
    MathOptModelStatistics statistics;

    statistics.number_of_variables = model.number_of_variables();
    for (int variable_id = 0; variable_id < model.number_of_variables(); ++variable_id) {
        switch (model.variables_types[variable_id]) {
        case VariableType::Continuous: {
            ++statistics.number_of_continuous_variables;
            break;
        } case VariableType::Binary: {
            ++statistics.number_of_binary_variables;
            break;
        } case VariableType::Integer: {
            if (model.variables_lower_bounds[variable_id] == 0.0
                    && model.variables_upper_bounds[variable_id] == 1.0) {
                ++statistics.number_of_binary_variables;
            } else {
                ++statistics.number_of_integer_variables;
            }
            break;
        } case VariableType::SemiContinuous: {
            ++statistics.number_of_semi_continuous_variables;
            break;
        } case VariableType::SemiInteger: {
            ++statistics.number_of_semi_integer_variables;
            break;
        }
        }
    }

    statistics.number_of_constraints = model.number_of_constraints();
    statistics.number_of_elements = model.number_of_elements();
    bool has_quadratic = !model.quadratic_elements_constraints_starts.empty();
    bool has_nonlinear = !model.nonlinear_elements_constraints_starts.empty();
    bool has_functions = !model.constraints_functions.empty();
    for (int constraint_id = 0; constraint_id < model.number_of_constraints(); ++constraint_id) {
        bool nonlinear =
            (has_nonlinear
             && model.nonlinear_elements_constraints_starts[constraint_id]
                < model.nonlinear_constraint_end(constraint_id))
            || (has_functions && bool(model.constraints_functions[constraint_id]));
        bool quadratic =
            has_quadratic
            && model.quadratic_elements_constraints_starts[constraint_id]
                < model.quadratic_constraint_end(constraint_id);
        if (nonlinear)
            ++statistics.number_of_nonlinear_constraints;
        else if (quadratic)
            ++statistics.number_of_quadratic_constraints;
    }

    bool objective_linear = false;
    for (double coefficient: model.objective_coefficients)
        if (coefficient != 0.0)
            objective_linear = true;
    bool objective_quadratic = !model.objective_quadratic_elements_coefficients.empty();
    bool objective_nonlinear = !model.objective_nonlinear_elements_operators.empty()
        || bool(model.objective_function);
    if (objective_linear || objective_quadratic || objective_nonlinear)
        statistics.number_of_objectives = 1;
    if (objective_nonlinear)
        statistics.number_of_nonlinear_objectives = 1;
    else if (objective_quadratic)
        statistics.number_of_quadratic_objectives = 1;

    statistics.number_of_defined_variables = model.number_of_defined_variables();
    return statistics;
}

MathOptModelStatistics mathoptsolverscmake::read_statistics(const std::string& file_path)
{
    size_t dot = file_path.find_last_of('.');
    std::string extension = (dot == std::string::npos)? "": file_path.substr(dot + 1);
    for (char& c: extension)
        c = (char)std::tolower((unsigned char)c);
    if (extension == "nl")
        return read_nl_statistics(file_path);
    if (extension == "mps")
        return compute_statistics(read_mps(file_path));
    if (extension == "lp")
        return compute_statistics(read_lp(file_path));
    if (extension == "json")
        return compute_statistics(read_json(file_path));
    throw std::invalid_argument(
            FUNC_SIGNATURE + ": "
            "unable to guess the format of \"" + file_path + "\" from its extension "
            "(.mps, .lp, .json or .nl).");
}

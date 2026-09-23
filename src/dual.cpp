#include "mathoptsolverscmake/dual.hpp"

#include <cmath>
#include <stdexcept>

using namespace mathoptsolverscmake;

namespace
{

const double inf = std::numeric_limits<double>::infinity();

/**
 * Dual variables associated with the two sides [lower, upper] of a primal
 * constraint or variable bound.
 *
 * Each side is described by the value of the bound (the dual objective
 * coefficient) and the bounds of the dual variable.
 */
struct Side
{
    double value;
    double lower_bound;
    double upper_bound;
    std::string suffix;
};

std::vector<Side> sides(
        double lower,
        double upper,
        double sign)
{
    std::vector<Side> output;
    if (lower == upper) {
        if (std::isfinite(lower))
            output.push_back({lower, -inf, inf, "_eq"});
        return output;
    }
    if (std::isfinite(lower)) {
        if (sign > 0) {
            output.push_back({lower, 0, inf, "_lb"});
        } else {
            output.push_back({lower, -inf, 0, "_lb"});
        }
    }
    if (std::isfinite(upper)) {
        if (sign > 0) {
            output.push_back({upper, -inf, 0, "_ub"});
        } else {
            output.push_back({upper, 0, inf, "_ub"});
        }
    }
    return output;
}

}

MathOptModel mathoptsolverscmake::dual(const MathOptModel& model)
{
    if (!model.is_lp()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "only linear programs (continuous variables, linear "
                "objective and constraints) are supported.");
    }

    // +1 for a minimization primal, -1 for a maximization primal.
    double sign = (model.objective_direction == ObjectiveDirection::Minimize)? 1: -1;

    MathOptModel dual_model;
    dual_model.objective_direction = (sign > 0)?
        ObjectiveDirection::Maximize:
        ObjectiveDirection::Minimize;
    dual_model.feasibility_tolerance = model.feasibility_tolerance;
    dual_model.integrality_tolerance = model.integrality_tolerance;

    auto add_variable = [&dual_model](
            const Side& side,
            const std::string& name)
    {
        dual_model.variables_lower_bounds.push_back(side.lower_bound);
        dual_model.variables_upper_bounds.push_back(side.upper_bound);
        dual_model.variables_types.push_back(VariableType::Continuous);
        dual_model.variables_names.push_back(name + side.suffix);
        dual_model.objective_coefficients.push_back(side.value);
    };

    // Dual variables of the primal constraints.
    std::vector<std::vector<int>> constraints_dual_variables(model.number_of_constraints());
    for (int constraint_id = 0;
            constraint_id < model.number_of_constraints();
            ++constraint_id) {
        for (const Side& side: sides(
                    model.constraints_lower_bounds[constraint_id],
                    model.constraints_upper_bounds[constraint_id],
                    sign)) {
            constraints_dual_variables[constraint_id].push_back(dual_model.number_of_variables());
            add_variable(side, model.constraint_name(constraint_id));
        }
    }

    // Dual variables of the primal variable bounds. Sides with a zero bound
    // are not added as variables, but eliminated into the sides of the
    // corresponding dual constraint.
    std::vector<std::vector<int>> variables_dual_variables(model.number_of_variables());
    dual_model.constraints_lower_bounds = model.objective_coefficients;
    dual_model.constraints_upper_bounds = model.objective_coefficients;
    for (int variable_id = 0;
            variable_id < model.number_of_variables();
            ++variable_id) {
        for (const Side& side: sides(
                    model.variables_lower_bounds[variable_id],
                    model.variables_upper_bounds[variable_id],
                    sign)) {
            if (side.value == 0) {
                // expr + z = c, z in [zl, zu] <=> expr in [c - zu, c - zl].
                dual_model.constraints_lower_bounds[variable_id] -= side.upper_bound;
                dual_model.constraints_upper_bounds[variable_id] -= side.lower_bound;
                continue;
            }
            variables_dual_variables[variable_id].push_back(dual_model.number_of_variables());
            add_variable(side, model.variable_name(variable_id));
        }
    }

    // Dual constraints: transpose of the primal constraint matrix, plus the
    // non-eliminated variable bound dual variables.
    std::vector<std::vector<int>> rows_variables(model.number_of_variables());
    std::vector<std::vector<double>> rows_coefficients(model.number_of_variables());
    for (int constraint_id = 0;
            constraint_id < model.number_of_constraints();
            ++constraint_id) {
        for (int element_id = model.constraints_starts[constraint_id];
                element_id < model.constraint_end(constraint_id);
                ++element_id) {
            int variable_id = model.elements_variables[element_id];
            double coefficient = model.elements_coefficients[element_id];
            for (int dual_variable_id: constraints_dual_variables[constraint_id]) {
                rows_variables[variable_id].push_back(dual_variable_id);
                rows_coefficients[variable_id].push_back(coefficient);
            }
        }
    }
    for (int variable_id = 0;
            variable_id < model.number_of_variables();
            ++variable_id) {
        dual_model.constraints_starts.push_back(dual_model.number_of_elements());
        dual_model.constraints_names.push_back(model.variable_name(variable_id));
        dual_model.elements_variables.insert(
                dual_model.elements_variables.end(),
                rows_variables[variable_id].begin(),
                rows_variables[variable_id].end());
        dual_model.elements_coefficients.insert(
                dual_model.elements_coefficients.end(),
                rows_coefficients[variable_id].begin(),
                rows_coefficients[variable_id].end());
        for (int dual_variable_id: variables_dual_variables[variable_id]) {
            dual_model.elements_variables.push_back(dual_variable_id);
            dual_model.elements_coefficients.push_back(1);
        }
    }
    dual_model.quadratic_elements_constraints_starts = std::vector<int>(model.number_of_variables(), 0);
    dual_model.constraints_functions = std::vector<BlackBoxFunction>(model.number_of_variables());

    return dual_model;
}

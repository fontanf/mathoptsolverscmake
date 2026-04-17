#include "mathoptsolverscmake/mathopt_knitro.hpp"

#include <limits>
#include <vector>
#include <numeric>

using namespace mathoptsolverscmake;

inline double to_knitro(double value)
{
    if (value == std::numeric_limits<double>::infinity())
        return KN_INFINITY;
    if (value == -std::numeric_limits<double>::infinity())
        return -KN_INFINITY;
    return value;
}

void mathoptsolverscmake::solve(
        const MathOptModel& model,
        knitrocpp::Context& knitro_context)
{
    // Variables.
    knitro_context.add_vars(model.number_of_variables());
    if (!model.variables_initial_values.empty())
        knitro_context.set_var_primal_init_values(model.variables_initial_values);
    std::vector<double> knitro_lower_bounds(model.number_of_variables());
    std::vector<double> knitro_upper_bounds(model.number_of_variables());
    for (int variable_id = 0; variable_id < model.number_of_variables(); ++variable_id) {
        knitro_lower_bounds[variable_id] = to_knitro(model.variables_lower_bounds[variable_id]);
        knitro_upper_bounds[variable_id] = to_knitro(model.variables_upper_bounds[variable_id]);
    }
    knitro_context.set_var_lobnds(knitro_lower_bounds);
    knitro_context.set_var_upbnds(knitro_upper_bounds);

    // Variable types (integer/binary).
    if (model.has_non_continuous_variables()) {
        std::vector<int> variable_types(model.number_of_variables(), KN_VARTYPE_CONTINUOUS);
        for (int variable_id = 0; variable_id < model.number_of_variables(); ++variable_id) {
            switch (model.variables_types[variable_id]) {
            case VariableType::Binary:
                variable_types[variable_id] = KN_VARTYPE_BINARY;
                break;
            case VariableType::Integer:
                variable_types[variable_id] = KN_VARTYPE_INTEGER;
                break;
            default:
                break;
            }
        }
        knitro_context.set_var_types(variable_types);
    }

    // Objective direction.
    knitro_context.set_obj_goal(
            model.objective_direction == ObjectiveDirection::Minimize?
            KN_OBJGOAL_MINIMIZE:
            KN_OBJGOAL_MAXIMIZE);

    // Linear objective.
    if (!model.objective_coefficients.empty()) {
        std::vector<int> variables(model.number_of_variables());
        std::iota(variables.begin(), variables.end(), 0);
        knitro_context.add_obj_linear_struct(
                variables,
                model.objective_coefficients);
    }

    // Quadratic objective.
    if (!model.objective_quadratic_elements_variables_1.empty()) {
        knitro_context.add_obj_quadratic_struct(
                model.objective_quadratic_elements_variables_1,
                model.objective_quadratic_elements_variables_2,
                model.objective_quadratic_elements_coefficients);
    }

    // Constraints.
    if (model.number_of_constraints() > 0) {
        knitro_context.add_cons(model.number_of_constraints());

        std::vector<double> knitro_constraints_lower_bounds(model.number_of_constraints());
        std::vector<double> knitro_constraints_upper_bounds(model.number_of_constraints());
        for (int constraint_id = 0; constraint_id < model.number_of_constraints(); ++constraint_id) {
            knitro_constraints_lower_bounds[constraint_id] = to_knitro(model.constraints_lower_bounds[constraint_id]);
            knitro_constraints_upper_bounds[constraint_id] = to_knitro(model.constraints_upper_bounds[constraint_id]);
        }
        knitro_context.set_con_lobnds(knitro_constraints_lower_bounds);
        knitro_context.set_con_upbnds(knitro_constraints_upper_bounds);

        // Linear constraints: convert CSR to COO.
        if (!model.elements_variables.empty()) {
            std::vector<int> constraints;
            constraints.reserve(model.number_of_elements());
            for (int constraint_id = 0; constraint_id < model.number_of_constraints(); ++constraint_id) {
                for (int el = model.constraints_starts[constraint_id];
                        el < model.constraint_end(constraint_id);
                        ++el) {
                    constraints.push_back(constraint_id);
                }
            }
            knitro_context.add_con_linear_struct(
                    constraints,
                    model.elements_variables,
                    model.elements_coefficients);
        }

        // Quadratic constraints: convert CSR to COO.
        if (!model.quadratic_elements_variables_1.empty()) {
            std::vector<int> constraints;
            constraints.reserve(model.quadratic_elements_variables_1.size());
            for (int constraint_id = 0; constraint_id < model.number_of_constraints(); ++constraint_id) {
                for (int el = model.quadratic_elements_constraints_starts[constraint_id];
                        el < model.quadratic_constraint_end(constraint_id);
                        ++el) {
                    constraints.push_back(constraint_id);
                }
            }
            knitro_context.add_con_quadratic_struct(
                    constraints,
                    model.quadratic_elements_variables_1,
                    model.quadratic_elements_variables_2,
                    model.quadratic_elements_coefficients);
        }
    }

    // Black-box objective callback.
    std::vector<double> x(model.number_of_variables(), 0.0);
    BlackBoxFunctionOutput obj_output;
    if (model.objective_function) {
        CB_context* callback_context = knitro_context.add_eval_callback(
                true,
                {},
                [&model, &x, &obj_output](
                    const knitrocpp::Context&,
                    CB_context*,
                    KN_eval_request_ptr const eval_request,
                    KN_eval_result_ptr const eval_result)
                {
                    for (int variable_id = 0; variable_id < model.number_of_variables(); ++variable_id)
                        x[variable_id] = eval_request->x[variable_id];
                    obj_output = model.objective_function(x);
                    *eval_result->obj = obj_output.objective_value;
                    return 0;
                });
        knitro_context.set_cb_grad(
                callback_context,
                [&model, &obj_output](
                    const knitrocpp::Context&,
                    CB_context*,
                    KN_eval_request_ptr const eval_request,
                    KN_eval_result_ptr const eval_result)
                {
                    for (int variable_id = 0; variable_id < model.number_of_variables(); ++variable_id)
                        eval_result->objGrad[variable_id] = obj_output.gradient[variable_id];
                    return 0;
                });
    }

    // Black-box constraint callbacks (one per black-box constraint).
    std::vector<BlackBoxFunctionOutput> con_outputs(model.number_of_constraints());
    for (int constraint_id = 0; constraint_id < model.number_of_constraints(); ++constraint_id) {
        if (model.constraints_functions.empty()
                || !model.constraints_functions[constraint_id]) {
            continue;
        }
        CB_context* callback_context = knitro_context.add_eval_callback(
                false,
                {constraint_id},
                [&model, &x, &con_outputs, constraint_id](
                    const knitrocpp::Context&,
                    CB_context*,
                    KN_eval_request_ptr const eval_request,
                    KN_eval_result_ptr const eval_result)
                {
                    for (int variable_id = 0; variable_id < model.number_of_variables(); ++variable_id)
                        x[variable_id] = eval_request->x[variable_id];
                    con_outputs[constraint_id] = model.constraints_functions[constraint_id](x);
                    eval_result->c[0] = con_outputs[constraint_id].objective_value;
                    return 0;
                });
        knitro_context.set_cb_grad(
                callback_context,
                [&model, &con_outputs, constraint_id](
                    const knitrocpp::Context&,
                    CB_context*,
                    KN_eval_request_ptr const eval_request,
                    KN_eval_result_ptr const eval_result)
                {
                    for (int variable_id = 0; variable_id < model.number_of_variables(); ++variable_id)
                        eval_result->jac[variable_id] = con_outputs[constraint_id].gradient[variable_id];
                    return 0;
                });
    }

    knitro_context.solve();
}

double mathoptsolverscmake::get_solution_value(
        const knitrocpp::Context& knitro_context)
{
    return knitro_context.get_obj_value();
}

std::vector<double> mathoptsolverscmake::get_solution(
        const knitrocpp::Context& knitro_context)
{
    return knitro_context.get_var_primal_values();
}

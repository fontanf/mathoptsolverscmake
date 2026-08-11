#include "mathoptsolverscmake/mathopt_knitro.hpp"

#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using namespace mathoptsolverscmake;

inline double to_knitro(double value)
{
    if (value == std::numeric_limits<double>::infinity())
        return KN_INFINITY;
    if (value == -std::numeric_limits<double>::infinity())
        return -KN_INFINITY;
    return value;
}

static void knitro_check(int rc, const std::string& caller, const char* fn)
{
    if (rc != 0)
        throw std::runtime_error(
                caller + ": "
                + fn + " failed with return code "
                + std::to_string(rc) + ".");
}

/**
 * Convert a nonlinear expression tree (stored as parallel element arrays,
 * see mathopt.hpp) into a Knitro expression graph and return the KNExprId
 * of its root.
 *
 * Mirrors the numerical evaluation performed by evaluate_nonlinear_element
 * in mathopt.cpp: nodes are stored in DFS pre-order (root first), so
 * scanning from the last element back to "start" is a post-order traversal
 * in which every child expression is built before its parent's.
 */
static KNExprId to_knitro_expr(
        KN_context_ptr kc,
        const std::vector<char>& operators,
        const std::vector<double>& values,
        const std::vector<int>& variables,
        const std::vector<int>& parent,
        int start,
        int end)
{
    const int n = end - start;

    // Arity is needed for n-ary + and *; derive it from the parent array.
    std::vector<int> arity(n, 0);
    for (int i = start + 1; i < end; ++i)
        arity[parent[i] - start]++;

    std::vector<KNExprId> stk;
    stk.reserve(n);
    for (int element_id = end - 1; element_id >= start; --element_id) {
        const int i = element_id - start;
        KNExprId expr = 0;
        switch (operators[element_id]) {
        case 'k':
            knitro_check(
                    KN_constant_expr(kc, values[element_id], &expr),
                    FUNC_SIGNATURE, "KN_constant_expr");
            break;
        case 'v':
            knitro_check(
                    KN_var_expr(kc, variables[element_id], &expr),
                    FUNC_SIGNATURE, "KN_var_expr");
            break;
        case 'n': {
            KNExprId a = stk.back(); stk.pop_back();
            knitro_check(KN_neg_expr(kc, a, &expr), FUNC_SIGNATURE, "KN_neg_expr");
            break;
        }
        case 'e': {
            KNExprId a = stk.back(); stk.pop_back();
            knitro_check(KN_exp_expr(kc, a, &expr), FUNC_SIGNATURE, "KN_exp_expr");
            break;
        }
        case 'l': {
            KNExprId a = stk.back(); stk.pop_back();
            knitro_check(KN_log_expr(kc, a, &expr), FUNC_SIGNATURE, "KN_log_expr");
            break;
        }
        case 'q': {
            KNExprId a = stk.back(); stk.pop_back();
            knitro_check(KN_sqrt_expr(kc, a, &expr), FUNC_SIGNATURE, "KN_sqrt_expr");
            break;
        }
        case 's': {
            KNExprId a = stk.back(); stk.pop_back();
            knitro_check(KN_sin_expr(kc, a, &expr), FUNC_SIGNATURE, "KN_sin_expr");
            break;
        }
        case 'c': {
            KNExprId a = stk.back(); stk.pop_back();
            knitro_check(KN_cos_expr(kc, a, &expr), FUNC_SIGNATURE, "KN_cos_expr");
            break;
        }
        case 't': {
            KNExprId a = stk.back(); stk.pop_back();
            knitro_check(KN_tan_expr(kc, a, &expr), FUNC_SIGNATURE, "KN_tan_expr");
            break;
        }
        case '+': {
            expr = stk.back(); stk.pop_back();
            for (int c = 1; c < arity[i]; ++c) {
                KNExprId term = stk.back(); stk.pop_back();
                KNExprId sum = 0;
                knitro_check(KN_add_expr(kc, expr, term, &sum), FUNC_SIGNATURE, "KN_add_expr");
                expr = sum;
            }
            break;
        }
        case '*': {
            expr = stk.back(); stk.pop_back();
            for (int c = 1; c < arity[i]; ++c) {
                KNExprId term = stk.back(); stk.pop_back();
                KNExprId prod = 0;
                knitro_check(KN_mul_expr(kc, expr, term, &prod), FUNC_SIGNATURE, "KN_mul_expr");
                expr = prod;
            }
            break;
        }
        case '-': {
            KNExprId l = stk.back(); stk.pop_back();
            KNExprId r = stk.back(); stk.pop_back();
            knitro_check(KN_sub_expr(kc, l, r, &expr), FUNC_SIGNATURE, "KN_sub_expr");
            break;
        }
        case '/': {
            KNExprId l = stk.back(); stk.pop_back();
            KNExprId r = stk.back(); stk.pop_back();
            knitro_check(KN_div_expr(kc, l, r, &expr), FUNC_SIGNATURE, "KN_div_expr");
            break;
        }
        case 'p': {
            KNExprId l = stk.back(); stk.pop_back();
            KNExprId r = stk.back(); stk.pop_back();
            knitro_check(KN_pow_expr(kc, l, r, &expr), FUNC_SIGNATURE, "KN_pow_expr");
            break;
        }
        default:
            throw std::runtime_error(
                    FUNC_SIGNATURE + ": "
                    "unknown nonlinear operator '"
                    + std::string(1, operators[element_id]) + "'.");
        }
        stk.push_back(expr);
    }
    return stk.back();
}

KnitroContext::KnitroContext()
{
    int rc = KN_new(&kc);
    if (rc != 0)
        throw std::runtime_error(
                FUNC_SIGNATURE + ": "
                "KN_new failed with return code " + std::to_string(rc) + ".");
    if (kc == NULL)
        throw std::runtime_error(
                FUNC_SIGNATURE + ": "
                "failed to find a valid Knitro license.");
}

KnitroContext::~KnitroContext()
{
    if (kc != NULL)
        KN_free(&kc);
}

static int eval_callback(
        KN_context_ptr,
        CB_context_ptr,
        KN_eval_request_ptr const eval_request,
        KN_eval_result_ptr const eval_result,
        void* const user_params)
{
    return static_cast<KnitroContext::EvalCallbackStruct*>(user_params)->eval(eval_request, eval_result);
}

static int grad_callback(
        KN_context_ptr,
        CB_context_ptr,
        KN_eval_request_ptr const eval_request,
        KN_eval_result_ptr const eval_result,
        void* const user_params)
{
    return static_cast<KnitroContext::EvalCallbackStruct*>(user_params)->grad(eval_request, eval_result);
}

void mathoptsolverscmake::reduce_printout(
        KnitroContext& knitro)
{
    knitro_check(
            KN_set_int_param(knitro.kc, KN_PARAM_OUTLEV, KN_OUTLEV_NONE),
            FUNC_SIGNATURE, "KN_set_int_param");
}

void mathoptsolverscmake::set_time_limit(
        KnitroContext& knitro,
        double time_limit)
{
    knitro_check(
            KN_set_double_param(knitro.kc, KN_PARAM_MAXTIMEREAL, time_limit),
            FUNC_SIGNATURE, "KN_set_double_param");
}

void mathoptsolverscmake::load(
        KnitroContext& knitro,
        const MathOptModel& model)
{
    // Variables.
    knitro_check(
            KN_add_vars(knitro.kc, model.number_of_variables(), NULL),
            FUNC_SIGNATURE, "KN_add_vars");
    if (!model.variables_initial_values.empty())
        knitro_check(
                KN_set_var_primal_init_values_all(
                        knitro.kc,
                        model.variables_initial_values.data()),
                FUNC_SIGNATURE, "KN_set_var_primal_init_values_all");
    std::vector<double> knitro_lower_bounds(model.number_of_variables());
    std::vector<double> knitro_upper_bounds(model.number_of_variables());
    for (int variable_id = 0; variable_id < model.number_of_variables(); ++variable_id) {
        knitro_lower_bounds[variable_id] = to_knitro(model.variables_lower_bounds[variable_id]);
        knitro_upper_bounds[variable_id] = to_knitro(model.variables_upper_bounds[variable_id]);
    }
    knitro_check(
            KN_set_var_lobnds_all(knitro.kc, knitro_lower_bounds.data()),
            FUNC_SIGNATURE, "KN_set_var_lobnds_all");
    knitro_check(
            KN_set_var_upbnds_all(knitro.kc, knitro_upper_bounds.data()),
            FUNC_SIGNATURE, "KN_set_var_upbnds_all");

    // Variable types (integer/binary).
    if (model.has_non_continuous_variables()) {
        std::vector<KNINT> variable_indices(model.number_of_variables());
        std::iota(variable_indices.begin(), variable_indices.end(), 0);
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
        knitro_check(
                KN_set_var_types(
                        knitro.kc,
                        model.number_of_variables(),
                        variable_indices.data(),
                        variable_types.data()),
                FUNC_SIGNATURE, "KN_set_var_types");
    }

    // Objective direction.
    knitro_check(
            KN_set_obj_goal(
                    knitro.kc,
                    model.objective_direction == ObjectiveDirection::Minimize?
                    KN_OBJGOAL_MINIMIZE:
                    KN_OBJGOAL_MAXIMIZE),
            FUNC_SIGNATURE, "KN_set_obj_goal");

    // Linear objective.
    if (!model.objective_coefficients.empty()) {
        std::vector<KNINT> variables(model.number_of_variables());
        std::iota(variables.begin(), variables.end(), 0);
        knitro_check(
                KN_add_obj_linear_struct(
                        knitro.kc,
                        model.number_of_variables(),
                        variables.data(),
                        model.objective_coefficients.data()),
                FUNC_SIGNATURE, "KN_add_obj_linear_struct");
    }

    // Quadratic objective.
    if (!model.objective_quadratic_elements_variables_1.empty()) {
        knitro_check(
                KN_add_obj_quadratic_struct(
                        knitro.kc,
                        (KNLONG)model.objective_quadratic_elements_variables_1.size(),
                        model.objective_quadratic_elements_variables_1.data(),
                        model.objective_quadratic_elements_variables_2.data(),
                        model.objective_quadratic_elements_coefficients.data()),
                FUNC_SIGNATURE, "KN_add_obj_quadratic_struct");
    }

    // Nonlinear objective.
    if (!model.objective_nonlinear_elements_operators.empty()) {
        KNExprId expr = to_knitro_expr(
                knitro.kc,
                model.objective_nonlinear_elements_operators,
                model.objective_nonlinear_elements_values,
                model.objective_nonlinear_elements_variables,
                model.objective_nonlinear_elements_parent,
                0,
                (int)model.objective_nonlinear_elements_operators.size());
        knitro_check(
                KN_add_obj_expr(knitro.kc, expr),
                FUNC_SIGNATURE, "KN_add_obj_expr");
    }

    // Constraints.
    if (model.number_of_constraints() > 0) {
        knitro_check(
                KN_add_cons(knitro.kc, model.number_of_constraints(), NULL),
                FUNC_SIGNATURE, "KN_add_cons");

        std::vector<double> knitro_constraints_lower_bounds(model.number_of_constraints());
        std::vector<double> knitro_constraints_upper_bounds(model.number_of_constraints());
        for (int constraint_id = 0; constraint_id < model.number_of_constraints(); ++constraint_id) {
            knitro_constraints_lower_bounds[constraint_id] = to_knitro(model.constraints_lower_bounds[constraint_id]);
            knitro_constraints_upper_bounds[constraint_id] = to_knitro(model.constraints_upper_bounds[constraint_id]);
        }
        knitro_check(
                KN_set_con_lobnds_all(knitro.kc, knitro_constraints_lower_bounds.data()),
                FUNC_SIGNATURE, "KN_set_con_lobnds_all");
        knitro_check(
                KN_set_con_upbnds_all(knitro.kc, knitro_constraints_upper_bounds.data()),
                FUNC_SIGNATURE, "KN_set_con_upbnds_all");

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
            knitro_check(
                    KN_add_con_linear_struct(
                            knitro.kc,
                            (KNLONG)model.number_of_elements(),
                            constraints.data(),
                            model.elements_variables.data(),
                            model.elements_coefficients.data()),
                    FUNC_SIGNATURE, "KN_add_con_linear_struct");
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
            knitro_check(
                    KN_add_con_quadratic_struct(
                            knitro.kc,
                            (KNLONG)model.quadratic_elements_variables_1.size(),
                            constraints.data(),
                            model.quadratic_elements_variables_1.data(),
                            model.quadratic_elements_variables_2.data(),
                            model.quadratic_elements_coefficients.data()),
                    FUNC_SIGNATURE, "KN_add_con_quadratic_struct");
        }

        // Nonlinear constraints.
        if (!model.nonlinear_elements_operators.empty()) {
            for (int constraint_id = 0;
                    constraint_id < model.number_of_constraints();
                    ++constraint_id) {
                int start = model.nonlinear_elements_constraints_starts[constraint_id];
                int end = model.nonlinear_constraint_end(constraint_id);
                if (end <= start)
                    continue;
                KNExprId expr = to_knitro_expr(
                        knitro.kc,
                        model.nonlinear_elements_operators,
                        model.nonlinear_elements_values,
                        model.nonlinear_elements_variables,
                        model.nonlinear_elements_parent,
                        start,
                        end);
                knitro_check(
                        KN_add_con_expr(knitro.kc, constraint_id, expr),
                        FUNC_SIGNATURE, "KN_add_con_expr");
            }
        }
    }

    // Resize reusable evaluation buffer.
    knitro.x_.resize(model.number_of_variables(), 0.0);

    // Black-box callbacks below share mutable state (knitro.x_ and a
    // per-callback output buffer) across calls, which is not safe under
    // Knitro's parallel MIP branch-and-bound (concurrent node subproblems
    // would call the same callback from multiple threads at once). Force
    // serial MIP solving whenever the model has a black-box function.
    if (model.has_black_box()) {
        knitro_check(
                KN_set_int_param(knitro.kc, KN_PARAM_MIP_NUMTHREADS, 1),
                FUNC_SIGNATURE, "KN_set_int_param");
    }

    // Black-box objective callback.
    if (model.objective_function) {
        auto output = std::make_shared<BlackBoxFunctionOutput>();
        std::unique_ptr<KnitroContext::EvalCallbackStruct> cb(new KnitroContext::EvalCallbackStruct());
        cb->eval = [&model, &knitro, output](
                KN_eval_request_ptr const eval_request,
                KN_eval_result_ptr const eval_result) {
            for (int i = 0; i < (int)knitro.x_.size(); ++i)
                knitro.x_[i] = eval_request->x[i];
            *output = model.objective_function(knitro.x_);
            *eval_result->obj = output->objective_value;
            return 0;
        };
        cb->grad = [output](
                KN_eval_request_ptr const,
                KN_eval_result_ptr const eval_result) {
            for (int i = 0; i < (int)output->gradient.size(); ++i)
                eval_result->objGrad[i] = output->gradient[i];
            return 0;
        };
        CB_context_ptr knitro_cb;
        knitro_check(
                KN_add_eval_callback(
                        knitro.kc, KNTRUE, 0, NULL,
                        eval_callback, &knitro_cb),
                FUNC_SIGNATURE, "KN_add_eval_callback");
        knitro_check(
                KN_set_cb_user_params(knitro.kc, knitro_cb, cb.get()),
                FUNC_SIGNATURE, "KN_set_cb_user_params");
        knitro_check(
                KN_set_cb_grad(
                        knitro.kc, knitro_cb, KN_DENSE, NULL, 0, NULL, NULL,
                        grad_callback),
                FUNC_SIGNATURE, "KN_set_cb_grad");
        knitro.eval_callbacks_.push_back(std::move(cb));
    }

    // Black-box constraint callbacks (one per black-box constraint).
    for (int constraint_id = 0; constraint_id < model.number_of_constraints(); ++constraint_id) {
        if (model.constraints_functions.empty()
                || !model.constraints_functions[constraint_id]) {
            continue;
        }
        auto output = std::make_shared<BlackBoxFunctionOutput>();
        std::unique_ptr<KnitroContext::EvalCallbackStruct> cb(new KnitroContext::EvalCallbackStruct());
        cb->eval = [&model, &knitro, output, constraint_id](
                KN_eval_request_ptr const eval_request,
                KN_eval_result_ptr const eval_result) {
            for (int i = 0; i < (int)knitro.x_.size(); ++i)
                knitro.x_[i] = eval_request->x[i];
            *output = model.constraints_functions[constraint_id](knitro.x_);
            eval_result->c[0] = output->objective_value;
            return 0;
        };
        cb->grad = [output](
                KN_eval_request_ptr const,
                KN_eval_result_ptr const eval_result) {
            for (int i = 0; i < (int)output->gradient.size(); ++i)
                eval_result->jac[i] = output->gradient[i];
            return 0;
        };
        CB_context_ptr knitro_cb;
        KNINT knitro_constraint_id = constraint_id;
        knitro_check(
                KN_add_eval_callback(
                        knitro.kc, KNFALSE, 1, &knitro_constraint_id,
                        eval_callback, &knitro_cb),
                FUNC_SIGNATURE, "KN_add_eval_callback");
        knitro_check(
                KN_set_cb_user_params(knitro.kc, knitro_cb, cb.get()),
                FUNC_SIGNATURE, "KN_set_cb_user_params");
        knitro_check(
                KN_set_cb_grad(
                        knitro.kc, knitro_cb, 0, NULL, KN_DENSE, NULL, NULL,
                        grad_callback),
                FUNC_SIGNATURE, "KN_set_cb_grad");
        knitro.eval_callbacks_.push_back(std::move(cb));
    }
}

void mathoptsolverscmake::solve(
        KnitroContext& knitro)
{
    KN_solve(knitro.kc);
}

double mathoptsolverscmake::get_solution_value(
        const KnitroContext& knitro)
{
    double obj = 0.0;
    knitro_check(
            KN_get_obj_value(knitro.kc, &obj),
            FUNC_SIGNATURE, "KN_get_obj_value");
    return obj;
}

std::vector<double> mathoptsolverscmake::get_solution(
        const KnitroContext& knitro)
{
    KNINT n = 0;
    knitro_check(
            KN_get_number_vars(knitro.kc, &n),
            FUNC_SIGNATURE, "KN_get_number_vars");
    std::vector<double> values(n, 0.0);
    knitro_check(
            KN_get_var_primal_values_all(knitro.kc, values.data()),
            FUNC_SIGNATURE, "KN_get_var_primal_values_all");
    return values;
}

int mathoptsolverscmake::get_number_of_nodes(
        const KnitroContext& knitro)
{
    int number_of_nodes = 0;
    knitro_check(
            KN_get_mip_number_nodes(knitro.kc, &number_of_nodes),
            FUNC_SIGNATURE, "KN_get_mip_number_nodes");
    return number_of_nodes;
}

bool mathoptsolverscmake::is_infeasible(
        const KnitroContext& knitro)
{
    KNINT n = 0;
    knitro_check(
            KN_get_number_vars(knitro.kc, &n),
            FUNC_SIGNATURE, "KN_get_number_vars");
    KNINT m = 0;
    knitro_check(
            KN_get_number_cons(knitro.kc, &m),
            FUNC_SIGNATURE, "KN_get_number_cons");

    int status = 0;
    double obj = 0.0;
    std::vector<double> x(n, 0.0);
    std::vector<double> lambda(n + m, 0.0);
    knitro_check(
            KN_get_solution(knitro.kc, &status, &obj, x.data(), lambda.data()),
            FUNC_SIGNATURE, "KN_get_solution");

    if (status <= KN_RC_INFEASIBLE && status > KN_RC_UNBOUNDED)
        return true;
    if (status == KN_RC_UNBOUNDED_OR_INFEAS)
        return true;
    if (status <= KN_RC_ITER_LIMIT_INFEAS && status >= KN_RC_MIP_NODE_LIMIT_INFEAS)
        return true;
    return false;
}

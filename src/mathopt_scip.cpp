#include "mathoptsolverscmake/mathopt_scip.hpp"

#include "scip/scipdefplugins.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace mathoptsolverscmake;

static void scip_check(SCIP_RETCODE retcode, const std::string& caller, const char* fn)
{
    if (retcode != SCIP_OKAY)
        throw std::runtime_error(
                caller + ": "
                + fn + " failed with return code "
                + std::to_string((int)retcode) + ".");
}

static double to_scip(SCIP* scip, double value)
{
    if (value >= SCIPinfinity(scip))
        return SCIPinfinity(scip);
    if (value <= -SCIPinfinity(scip))
        return -SCIPinfinity(scip);
    return value;
}

static double from_scip(SCIP* scip, double value)
{
    if (SCIPisInfinity(scip, value))
        return std::numeric_limits<double>::infinity();
    if (SCIPisInfinity(scip, -value))
        return -std::numeric_limits<double>::infinity();
    return value;
}

/**
 * Convert a nonlinear expression tree (stored as parallel element arrays,
 * see mathopt.hpp) into a SCIP expression and return its root. The caller
 * owns the returned expression and must release it.
 *
 * Nodes are stored in DFS pre-order (root first), so scanning from the last
 * element back to "start" is a post-order traversal in which every child
 * expression is built before its parent's (see to_knitro_expr in
 * mathopt_knitro.cpp).
 *
 * SCIP has no tan, division or variable-exponent power expressions:
 * - tan(a) is built as sin(a) * cos(a)^-1
 * - a / b as a * b^-1
 * - a ^ b requires a constant exponent b, or a constant base a > 0, in which
 *   case it is built as exp(b * log(a))
 */
static SCIP_EXPR* to_scip_expr(
        SCIP* scip,
        const std::vector<SCIP_VAR*>& scip_variables,
        const std::vector<char>& operators,
        const std::vector<double>& values,
        const std::vector<int>& variables,
        const std::vector<int>& number_of_children,
        int start,
        int end)
{
    // Pop the "number_of_operands" expressions on top of the stack, first
    // popped first.
    std::vector<SCIP_EXPR*> stk;
    stk.reserve(end - start);
    auto pop = [&stk](int number_of_operands)
    {
        std::vector<SCIP_EXPR*> operands(number_of_operands);
        for (int operand_pos = 0; operand_pos < number_of_operands; ++operand_pos) {
            operands[operand_pos] = stk.back();
            stk.pop_back();
        }
        return operands;
    };
    auto release = [scip](std::vector<SCIP_EXPR*>& exprs)
    {
        for (SCIP_EXPR*& expr: exprs)
            scip_check(SCIPreleaseExpr(scip, &expr), FUNC_SIGNATURE, "SCIPreleaseExpr");
    };

    for (int element_id = end - 1; element_id >= start; --element_id) {
        SCIP_EXPR* expr = nullptr;
        switch (operators[element_id]) {
        case 'k': {
            scip_check(
                    SCIPcreateExprValue(scip, &expr, values[element_id], NULL, NULL),
                    FUNC_SIGNATURE, "SCIPcreateExprValue");
            break;
        } case 'v': {
            scip_check(
                    SCIPcreateExprVar(scip, &expr, scip_variables[variables[element_id]], NULL, NULL),
                    FUNC_SIGNATURE, "SCIPcreateExprVar");
            break;
        } case 'n': {
            std::vector<SCIP_EXPR*> operands = pop(1);
            double coefficient = -1.0;
            scip_check(
                    SCIPcreateExprSum(scip, &expr, 1, operands.data(), &coefficient, 0.0, NULL, NULL),
                    FUNC_SIGNATURE, "SCIPcreateExprSum");
            release(operands);
            break;
        } case 'e': {
            std::vector<SCIP_EXPR*> operands = pop(1);
            scip_check(
                    SCIPcreateExprExp(scip, &expr, operands[0], NULL, NULL),
                    FUNC_SIGNATURE, "SCIPcreateExprExp");
            release(operands);
            break;
        } case 'l': {
            std::vector<SCIP_EXPR*> operands = pop(1);
            scip_check(
                    SCIPcreateExprLog(scip, &expr, operands[0], NULL, NULL),
                    FUNC_SIGNATURE, "SCIPcreateExprLog");
            release(operands);
            break;
        } case 'q': {
            std::vector<SCIP_EXPR*> operands = pop(1);
            scip_check(
                    SCIPcreateExprPow(scip, &expr, operands[0], 0.5, NULL, NULL),
                    FUNC_SIGNATURE, "SCIPcreateExprPow");
            release(operands);
            break;
        } case 's': {
            std::vector<SCIP_EXPR*> operands = pop(1);
            scip_check(
                    SCIPcreateExprSin(scip, &expr, operands[0], NULL, NULL),
                    FUNC_SIGNATURE, "SCIPcreateExprSin");
            release(operands);
            break;
        } case 'c': {
            std::vector<SCIP_EXPR*> operands = pop(1);
            scip_check(
                    SCIPcreateExprCos(scip, &expr, operands[0], NULL, NULL),
                    FUNC_SIGNATURE, "SCIPcreateExprCos");
            release(operands);
            break;
        } case 't': {
            std::vector<SCIP_EXPR*> operands = pop(1);
            std::vector<SCIP_EXPR*> factors(2, nullptr);
            SCIP_EXPR* cos_expr = nullptr;
            scip_check(
                    SCIPcreateExprSin(scip, &factors[0], operands[0], NULL, NULL),
                    FUNC_SIGNATURE, "SCIPcreateExprSin");
            scip_check(
                    SCIPcreateExprCos(scip, &cos_expr, operands[0], NULL, NULL),
                    FUNC_SIGNATURE, "SCIPcreateExprCos");
            scip_check(
                    SCIPcreateExprPow(scip, &factors[1], cos_expr, -1.0, NULL, NULL),
                    FUNC_SIGNATURE, "SCIPcreateExprPow");
            scip_check(SCIPreleaseExpr(scip, &cos_expr), FUNC_SIGNATURE, "SCIPreleaseExpr");
            scip_check(
                    SCIPcreateExprProduct(scip, &expr, 2, factors.data(), 1.0, NULL, NULL),
                    FUNC_SIGNATURE, "SCIPcreateExprProduct");
            release(factors);
            release(operands);
            break;
        } case '+': {
            std::vector<SCIP_EXPR*> operands = pop(number_of_children[element_id]);
            scip_check(
                    SCIPcreateExprSum(scip, &expr, (int)operands.size(), operands.data(), NULL, 0.0, NULL, NULL),
                    FUNC_SIGNATURE, "SCIPcreateExprSum");
            release(operands);
            break;
        } case '*': {
            std::vector<SCIP_EXPR*> operands = pop(number_of_children[element_id]);
            scip_check(
                    SCIPcreateExprProduct(scip, &expr, (int)operands.size(), operands.data(), 1.0, NULL, NULL),
                    FUNC_SIGNATURE, "SCIPcreateExprProduct");
            release(operands);
            break;
        } case '-': {
            std::vector<SCIP_EXPR*> operands = pop(2);
            double coefficients[2] = {1.0, -1.0};
            scip_check(
                    SCIPcreateExprSum(scip, &expr, 2, operands.data(), coefficients, 0.0, NULL, NULL),
                    FUNC_SIGNATURE, "SCIPcreateExprSum");
            release(operands);
            break;
        } case '/': {
            std::vector<SCIP_EXPR*> operands = pop(2);
            std::vector<SCIP_EXPR*> factors(2, nullptr);
            factors[0] = operands[0];
            scip_check(
                    SCIPcreateExprPow(scip, &factors[1], operands[1], -1.0, NULL, NULL),
                    FUNC_SIGNATURE, "SCIPcreateExprPow");
            scip_check(
                    SCIPcreateExprProduct(scip, &expr, 2, factors.data(), 1.0, NULL, NULL),
                    FUNC_SIGNATURE, "SCIPcreateExprProduct");
            scip_check(SCIPreleaseExpr(scip, &factors[1]), FUNC_SIGNATURE, "SCIPreleaseExpr");
            release(operands);
            break;
        } case 'p': {
            std::vector<SCIP_EXPR*> operands = pop(2);
            if (SCIPisExprValue(scip, operands[1])) {
                scip_check(
                        SCIPcreateExprPow(
                                scip, &expr, operands[0],
                                SCIPgetValueExprValue(operands[1]),
                                NULL, NULL),
                        FUNC_SIGNATURE, "SCIPcreateExprPow");
            } else if (SCIPisExprValue(scip, operands[0])
                    && SCIPgetValueExprValue(operands[0]) > 0) {
                // a ^ b = exp(log(a) * b).
                double log_base = std::log(SCIPgetValueExprValue(operands[0]));
                SCIP_EXPR* exponent_expr = nullptr;
                scip_check(
                        SCIPcreateExprSum(scip, &exponent_expr, 1, &operands[1], &log_base, 0.0, NULL, NULL),
                        FUNC_SIGNATURE, "SCIPcreateExprSum");
                scip_check(
                        SCIPcreateExprExp(scip, &expr, exponent_expr, NULL, NULL),
                        FUNC_SIGNATURE, "SCIPcreateExprExp");
                scip_check(SCIPreleaseExpr(scip, &exponent_expr), FUNC_SIGNATURE, "SCIPreleaseExpr");
            } else {
                release(operands);
                throw std::invalid_argument(
                        FUNC_SIGNATURE + ": "
                        "SCIP only supports powers with a constant exponent "
                        "or a constant positive base.");
            }
            release(operands);
            break;
        } default: {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "unknown nonlinear operator '"
                    + std::string(1, operators[element_id]) + "'.");
        }
        }
        stk.push_back(expr);
    }
    return stk.back();
}

ScipModel::ScipModel()
{
    scip_check(SCIPcreate(&scip), FUNC_SIGNATURE, "SCIPcreate");
    scip_check(SCIPincludeDefaultPlugins(scip), FUNC_SIGNATURE, "SCIPincludeDefaultPlugins");
    scip_check(SCIPcreateProbBasic(scip, "mathopt"), FUNC_SIGNATURE, "SCIPcreateProbBasic");
}

ScipModel::~ScipModel()
{
    if (scip != nullptr)
        SCIPfree(&scip);
}

void mathoptsolverscmake::reduce_printout(
        ScipModel& scip_model)
{
    SCIPsetMessagehdlrQuiet(scip_model.scip, TRUE);
}

void mathoptsolverscmake::set_time_limit(
        ScipModel& scip_model,
        double time_limit)
{
    scip_check(
            SCIPsetRealParam(scip_model.scip, "limits/time", time_limit),
            FUNC_SIGNATURE, "SCIPsetRealParam");
}

void mathoptsolverscmake::set_node_limit(
        ScipModel& scip_model,
        int node_limit)
{
    scip_check(
            SCIPsetLongintParam(scip_model.scip, "limits/nodes", node_limit),
            FUNC_SIGNATURE, "SCIPsetLongintParam");
}

void mathoptsolverscmake::set_log_file(
        ScipModel& scip_model,
        const std::string& log_file)
{
    SCIPsetMessagehdlrLogfile(scip_model.scip, log_file.c_str());
}

void mathoptsolverscmake::load(
        ScipModel& scip_model,
        const MathOptModel& model)
{
    if (model.has_black_box()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "SCIP does not support black-box functions.");
    }

    SCIP* scip = scip_model.scip;

    // Objective direction.
    scip_check(
            SCIPsetObjsense(
                    scip,
                    (model.objective_direction == ObjectiveDirection::Minimize)?
                    SCIP_OBJSENSE_MINIMIZE:
                    SCIP_OBJSENSE_MAXIMIZE),
            FUNC_SIGNATURE, "SCIPsetObjsense");

    // Variables, with the linear objective.
    scip_model.variables.resize(model.number_of_variables(), nullptr);
    for (int variable_id = 0; variable_id < model.number_of_variables(); ++variable_id) {
        double lower_bound = model.variables_lower_bounds[variable_id];
        double upper_bound = model.variables_upper_bounds[variable_id];
        SCIP_VARTYPE variable_type = SCIP_VARTYPE_CONTINUOUS;
        switch (model.variables_types[variable_id]) {
        case VariableType::Binary: {
            variable_type = SCIP_VARTYPE_BINARY;
            lower_bound = (std::max)(lower_bound, 0.0);
            upper_bound = (std::min)(upper_bound, 1.0);
            break;
        } case VariableType::Integer: {
            variable_type = SCIP_VARTYPE_INTEGER;
            break;
        } case VariableType::Continuous: {
            break;
        }
        }
        double objective_coefficient = model.objective_coefficients.empty()?
            0.0:
            model.objective_coefficients[variable_id];
        SCIP_VAR* var = nullptr;
        scip_check(
                SCIPcreateVarBasic(
                        scip,
                        &var,
                        model.variable_name(variable_id).c_str(),
                        to_scip(scip, lower_bound),
                        to_scip(scip, upper_bound),
                        objective_coefficient,
                        variable_type),
                FUNC_SIGNATURE, "SCIPcreateVarBasic");
        scip_check(SCIPaddVar(scip, var), FUNC_SIGNATURE, "SCIPaddVar");
        scip_model.variables[variable_id] = var;
        // The problem keeps its own reference to the variable.
        scip_check(SCIPreleaseVar(scip, &var), FUNC_SIGNATURE, "SCIPreleaseVar");
    }

    // Build the expression of the quadratic and nonlinear terms of a
    // constraint or of the objective. Return nullptr if there are none.
    auto build_nonlinear_expr = [scip, &scip_model](
            const std::vector<int>& quadratic_variables_1,
            const std::vector<int>& quadratic_variables_2,
            const std::vector<double>& quadratic_coefficients,
            int quadratic_start,
            int quadratic_end,
            const std::vector<char>& operators,
            const std::vector<double>& values,
            const std::vector<int>& variables,
            const std::vector<int>& number_of_children,
            int nonlinear_start,
            int nonlinear_end) -> SCIP_EXPR*
    {
        std::vector<SCIP_EXPR*> terms;
        if (quadratic_end > quadratic_start) {
            std::vector<SCIP_VAR*> vars_1;
            std::vector<SCIP_VAR*> vars_2;
            std::vector<double> coefficients;
            for (int element_id = quadratic_start; element_id < quadratic_end; ++element_id) {
                vars_1.push_back(scip_model.variables[quadratic_variables_1[element_id]]);
                vars_2.push_back(scip_model.variables[quadratic_variables_2[element_id]]);
                coefficients.push_back(quadratic_coefficients[element_id]);
            }
            SCIP_EXPR* expr = nullptr;
            scip_check(
                    SCIPcreateExprQuadratic(
                            scip, &expr,
                            0, NULL, NULL,
                            (int)coefficients.size(),
                            vars_1.data(),
                            vars_2.data(),
                            coefficients.data(),
                            NULL, NULL),
                    FUNC_SIGNATURE, "SCIPcreateExprQuadratic");
            terms.push_back(expr);
        }
        if (nonlinear_end > nonlinear_start) {
            terms.push_back(to_scip_expr(
                        scip,
                        scip_model.variables,
                        operators,
                        values,
                        variables,
                        number_of_children,
                        nonlinear_start,
                        nonlinear_end));
        }
        if (terms.empty())
            return nullptr;
        if (terms.size() == 1)
            return terms[0];
        SCIP_EXPR* expr = nullptr;
        scip_check(
                SCIPcreateExprSum(scip, &expr, (int)terms.size(), terms.data(), NULL, 0.0, NULL, NULL),
                FUNC_SIGNATURE, "SCIPcreateExprSum");
        for (SCIP_EXPR*& term: terms)
            scip_check(SCIPreleaseExpr(scip, &term), FUNC_SIGNATURE, "SCIPreleaseExpr");
        return expr;
    };

    // Add a constraint lhs <= linear terms + expr <= rhs, expr may be null.
    auto add_constraint = [scip](
            const std::string& name,
            std::vector<SCIP_VAR*>& linear_vars,
            std::vector<double>& linear_coefficients,
            SCIP_EXPR* expr,
            double lhs,
            double rhs)
    {
        SCIP_CONS* cons = nullptr;
        if (expr == nullptr) {
            scip_check(
                    SCIPcreateConsBasicLinear(
                            scip, &cons, name.c_str(),
                            (int)linear_vars.size(),
                            linear_vars.data(),
                            linear_coefficients.data(),
                            lhs, rhs),
                    FUNC_SIGNATURE, "SCIPcreateConsBasicLinear");
        } else {
            scip_check(
                    SCIPcreateConsBasicNonlinear(scip, &cons, name.c_str(), expr, lhs, rhs),
                    FUNC_SIGNATURE, "SCIPcreateConsBasicNonlinear");
            scip_check(SCIPreleaseExpr(scip, &expr), FUNC_SIGNATURE, "SCIPreleaseExpr");
            for (int pos = 0; pos < (int)linear_vars.size(); ++pos) {
                scip_check(
                        SCIPaddLinearVarNonlinear(scip, cons, linear_vars[pos], linear_coefficients[pos]),
                        FUNC_SIGNATURE, "SCIPaddLinearVarNonlinear");
            }
        }
        scip_check(SCIPaddCons(scip, cons), FUNC_SIGNATURE, "SCIPaddCons");
        scip_check(SCIPreleaseCons(scip, &cons), FUNC_SIGNATURE, "SCIPreleaseCons");
    };

    // Quadratic and nonlinear objective terms: minimize c x + z s.t.
    // f(x) - z <= 0 (or maximize c x + z s.t. f(x) - z >= 0).
    SCIP_EXPR* objective_expr = build_nonlinear_expr(
            model.objective_quadratic_elements_variables_1,
            model.objective_quadratic_elements_variables_2,
            model.objective_quadratic_elements_coefficients,
            0,
            (int)model.objective_quadratic_elements_variables_1.size(),
            model.objective_nonlinear_elements_operators,
            model.objective_nonlinear_elements_values,
            model.objective_nonlinear_elements_variables,
            model.objective_nonlinear_elements_number_of_children,
            0,
            (int)model.objective_nonlinear_elements_operators.size());
    if (objective_expr != nullptr) {
        SCIP_VAR* var = nullptr;
        scip_check(
                SCIPcreateVarBasic(
                        scip,
                        &var,
                        "objective",
                        -SCIPinfinity(scip),
                        SCIPinfinity(scip),
                        1.0,
                        SCIP_VARTYPE_CONTINUOUS),
                FUNC_SIGNATURE, "SCIPcreateVarBasic");
        scip_check(SCIPaddVar(scip, var), FUNC_SIGNATURE, "SCIPaddVar");
        scip_model.objective_variable = var;
        scip_check(SCIPreleaseVar(scip, &var), FUNC_SIGNATURE, "SCIPreleaseVar");

        std::vector<SCIP_VAR*> linear_vars = {scip_model.objective_variable};
        std::vector<double> linear_coefficients = {-1.0};
        bool minimize = (model.objective_direction == ObjectiveDirection::Minimize);
        add_constraint(
                "objective",
                linear_vars,
                linear_coefficients,
                objective_expr,
                minimize? -SCIPinfinity(scip): 0.0,
                minimize? 0.0: SCIPinfinity(scip));
    }

    // Constraints.
    for (int constraint_id = 0; constraint_id < model.number_of_constraints(); ++constraint_id) {
        std::vector<SCIP_VAR*> linear_vars;
        std::vector<double> linear_coefficients;
        if (!model.elements_variables.empty()) {
            for (int element_id = model.constraints_starts[constraint_id];
                    element_id < model.constraint_end(constraint_id);
                    ++element_id) {
                linear_vars.push_back(scip_model.variables[model.elements_variables[element_id]]);
                linear_coefficients.push_back(model.elements_coefficients[element_id]);
            }
        }

        int quadratic_start = 0;
        int quadratic_end = 0;
        if (!model.quadratic_elements_variables_1.empty()) {
            quadratic_start = model.quadratic_elements_constraints_starts[constraint_id];
            quadratic_end = model.quadratic_constraint_end(constraint_id);
        }
        int nonlinear_start = 0;
        int nonlinear_end = 0;
        if (!model.nonlinear_elements_operators.empty()) {
            nonlinear_start = model.nonlinear_elements_constraints_starts[constraint_id];
            nonlinear_end = model.nonlinear_constraint_end(constraint_id);
        }
        SCIP_EXPR* expr = build_nonlinear_expr(
                model.quadratic_elements_variables_1,
                model.quadratic_elements_variables_2,
                model.quadratic_elements_coefficients,
                quadratic_start,
                quadratic_end,
                model.nonlinear_elements_operators,
                model.nonlinear_elements_values,
                model.nonlinear_elements_variables,
                model.nonlinear_elements_number_of_children,
                nonlinear_start,
                nonlinear_end);

        add_constraint(
                model.constraint_name(constraint_id),
                linear_vars,
                linear_coefficients,
                expr,
                to_scip(scip, model.constraints_lower_bounds[constraint_id]),
                to_scip(scip, model.constraints_upper_bounds[constraint_id]));
    }

    // Initial solution.
    if (!model.variables_initial_values.empty()) {
        SCIP_SOL* sol = nullptr;
        scip_check(SCIPcreateOrigSol(scip, &sol, NULL), FUNC_SIGNATURE, "SCIPcreateOrigSol");
        for (int variable_id = 0; variable_id < model.number_of_variables(); ++variable_id) {
            scip_check(
                    SCIPsetSolVal(
                            scip, sol,
                            scip_model.variables[variable_id],
                            model.variables_initial_values[variable_id]),
                    FUNC_SIGNATURE, "SCIPsetSolVal");
        }
        if (scip_model.objective_variable != nullptr) {
            // The value of z is the value of the non-linear part of the
            // objective.
            double value = model.evaluate_objective(model.variables_initial_values);
            if (!model.objective_coefficients.empty()) {
                for (int variable_id = 0; variable_id < model.number_of_variables(); ++variable_id) {
                    value -= model.objective_coefficients[variable_id]
                        * model.variables_initial_values[variable_id];
                }
            }
            scip_check(
                    SCIPsetSolVal(scip, sol, scip_model.objective_variable, value),
                    FUNC_SIGNATURE, "SCIPsetSolVal");
        }
        SCIP_Bool stored = FALSE;
        scip_check(SCIPaddSolFree(scip, &sol, &stored), FUNC_SIGNATURE, "SCIPaddSolFree");
    }
}

void mathoptsolverscmake::solve(
        ScipModel& scip_model)
{
    scip_check(SCIPsolve(scip_model.scip), FUNC_SIGNATURE, "SCIPsolve");
}

double mathoptsolverscmake::get_solution_value(
        const ScipModel& scip_model)
{
    return from_scip(scip_model.scip, SCIPgetPrimalbound(scip_model.scip));
}

std::vector<double> mathoptsolverscmake::get_solution(
        const ScipModel& scip_model)
{
    SCIP_SOL* sol = SCIPgetBestSol(scip_model.scip);
    if (sol == nullptr)
        return {};
    std::vector<double> solution(scip_model.variables.size());
    for (int variable_id = 0; variable_id < (int)scip_model.variables.size(); ++variable_id) {
        solution[variable_id] = SCIPgetSolVal(
                scip_model.scip,
                sol,
                scip_model.variables[variable_id]);
    }
    return solution;
}

double mathoptsolverscmake::get_bound(
        const ScipModel& scip_model)
{
    return from_scip(scip_model.scip, SCIPgetDualbound(scip_model.scip));
}

int mathoptsolverscmake::get_number_of_nodes(
        const ScipModel& scip_model)
{
    if (SCIPgetStage(scip_model.scip) < SCIP_STAGE_TRANSFORMED)
        return 0;
    return (int)SCIPgetNNodes(scip_model.scip);
}

bool mathoptsolverscmake::is_infeasible(
        const ScipModel& scip_model)
{
    SCIP_STATUS status = SCIPgetStatus(scip_model.scip);
    return status == SCIP_STATUS_INFEASIBLE
        || status == SCIP_STATUS_INFORUNBD;
}

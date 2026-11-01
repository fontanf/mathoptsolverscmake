#include "mathoptsolverscmake/mathopt_cplex.hpp"

#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace mathoptsolverscmake;

static void cplex_check(CPXCENVptr env, int status, const std::string& caller, const char* fn)
{
    if (status != 0) {
        char buffer[CPXMESSAGEBUFSIZE];
        const char* message = (env != nullptr)?
            CPXgeterrorstring(env, status, buffer):
            nullptr;
        throw std::runtime_error(
                caller + ": "
                + fn + " failed with status "
                + std::to_string(status)
                + ((message != nullptr)? std::string(": ") + message: std::string(""))
                + ".");
    }
}

static double to_cplex(double value)
{
    if (value >= CPX_INFBOUND)
        return CPX_INFBOUND;
    if (value <= -CPX_INFBOUND)
        return -CPX_INFBOUND;
    return value;
}

static bool is_mip(const CplexModel& cplex_model)
{
    switch (CPXgetprobtype(cplex_model.env, cplex_model.lp)) {
    case CPXPROB_LP:
    case CPXPROB_QP:
    case CPXPROB_QCP:
        return false;
    default:
        return true;
    }
}

CplexModel::CplexModel()
{
    int status = 0;
    env = CPXopenCPLEX(&status);
    cplex_check(nullptr, status, FUNC_SIGNATURE, "CPXopenCPLEX");
    // The output is disabled by default in the callable library.
    cplex_check(
            env,
            CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON),
            FUNC_SIGNATURE, "CPXsetintparam");
    lp = CPXcreateprob(env, &status, "mathopt");
    cplex_check(env, status, FUNC_SIGNATURE, "CPXcreateprob");
}

CplexModel::~CplexModel()
{
    if (lp != nullptr)
        CPXfreeprob(env, &lp);
    if (env != nullptr)
        CPXcloseCPLEX(&env);
}

void mathoptsolverscmake::reduce_printout(
        CplexModel& cplex_model)
{
    cplex_check(
            cplex_model.env,
            CPXsetintparam(cplex_model.env, CPXPARAM_ScreenOutput, CPX_OFF),
            FUNC_SIGNATURE, "CPXsetintparam");
}

void mathoptsolverscmake::set_time_limit(
        CplexModel& cplex_model,
        double time_limit)
{
    cplex_check(
            cplex_model.env,
            CPXsetdblparam(cplex_model.env, CPXPARAM_TimeLimit, time_limit),
            FUNC_SIGNATURE, "CPXsetdblparam");
}

void mathoptsolverscmake::set_node_limit(
        CplexModel& cplex_model,
        int node_limit)
{
    cplex_check(
            cplex_model.env,
            CPXsetlongparam(cplex_model.env, CPXPARAM_MIP_Limits_Nodes, node_limit),
            FUNC_SIGNATURE, "CPXsetlongparam");
}

void mathoptsolverscmake::set_log_file(
        CplexModel& cplex_model,
        const std::string& log_file)
{
    cplex_check(
            cplex_model.env,
            CPXsetlogfilename(cplex_model.env, log_file.c_str(), "w"),
            FUNC_SIGNATURE, "CPXsetlogfilename");
}

void mathoptsolverscmake::load(
        CplexModel& cplex_model,
        const MathOptModel& model)
{
    if (model.has_black_box()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "CPLEX does not support black-box functions.");
    }
    if (model.has_nonlinear()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "CPLEX does not support nonlinear expression trees.");
    }

    CPXENVptr env = cplex_model.env;
    CPXLPptr lp = cplex_model.lp;

    cplex_check(
            env,
            CPXchgobjsen(
                    env, lp,
                    (model.objective_direction == ObjectiveDirection::Minimize)?
                    CPX_MIN:
                    CPX_MAX),
            FUNC_SIGNATURE, "CPXchgobjsen");

    // Variables, with the linear objective.
    int number_of_variables = model.number_of_variables();
    {
        std::vector<double> objective(number_of_variables, 0.0);
        std::vector<double> lower_bounds(number_of_variables);
        std::vector<double> upper_bounds(number_of_variables);
        std::vector<char> types(number_of_variables);
        std::vector<std::string> names(number_of_variables);
        std::vector<char*> names_ptrs(number_of_variables);
        for (int variable_id = 0; variable_id < number_of_variables; ++variable_id) {
            if (!model.objective_coefficients.empty())
                objective[variable_id] = model.objective_coefficients[variable_id];
            lower_bounds[variable_id] = to_cplex(model.variables_lower_bounds[variable_id]);
            upper_bounds[variable_id] = to_cplex(model.variables_upper_bounds[variable_id]);
            switch (model.variables_types[variable_id]) {
            case VariableType::Continuous: types[variable_id] = CPX_CONTINUOUS; break;
            case VariableType::Binary: types[variable_id] = CPX_BINARY; break;
            case VariableType::Integer: types[variable_id] = CPX_INTEGER; break;
            }
            names[variable_id] = model.variable_name(variable_id);
            names_ptrs[variable_id] = const_cast<char*>(names[variable_id].c_str());
        }
        // Without variable types, the problem remains an LP.
        cplex_check(
                env,
                CPXnewcols(
                        env, lp,
                        number_of_variables,
                        objective.data(),
                        lower_bounds.data(),
                        upper_bounds.data(),
                        model.has_non_continuous_variables()? types.data(): NULL,
                        names_ptrs.data()),
                FUNC_SIGNATURE, "CPXnewcols");
    }

    // Quadratic objective. CPLEX minimizes c x + 1/2 x'Qx with Q symmetric,
    // given column by column.
    if (!model.objective_quadratic_elements_variables_1.empty()) {
        std::vector<std::map<int, double>> q_columns(number_of_variables);
        for (int element_id = 0;
                element_id < (int)model.objective_quadratic_elements_variables_1.size();
                ++element_id) {
            int variable_id_1 = model.objective_quadratic_elements_variables_1[element_id];
            int variable_id_2 = model.objective_quadratic_elements_variables_2[element_id];
            double coefficient = model.objective_quadratic_elements_coefficients[element_id];
            if (variable_id_1 == variable_id_2) {
                q_columns[variable_id_1][variable_id_1] += 2 * coefficient;
            } else {
                q_columns[variable_id_1][variable_id_2] += coefficient;
                q_columns[variable_id_2][variable_id_1] += coefficient;
            }
        }
        std::vector<int> q_starts(number_of_variables);
        std::vector<int> q_counts(number_of_variables);
        std::vector<int> q_rows;
        std::vector<double> q_values;
        for (int variable_id = 0; variable_id < number_of_variables; ++variable_id) {
            q_starts[variable_id] = (int)q_rows.size();
            q_counts[variable_id] = (int)q_columns[variable_id].size();
            for (const auto& entry: q_columns[variable_id]) {
                q_rows.push_back(entry.first);
                q_values.push_back(entry.second);
            }
        }
        cplex_check(
                env,
                CPXcopyquad(
                        env, lp,
                        q_starts.data(),
                        q_counts.data(),
                        q_rows.data(),
                        q_values.data()),
                FUNC_SIGNATURE, "CPXcopyquad");
        cplex_check(
                env,
                CPXsetintparam(
                        env,
                        CPXPARAM_OptimalityTarget,
                        CPX_OPTIMALITYTARGET_OPTIMALGLOBAL),
                FUNC_SIGNATURE, "CPXsetintparam");
    }

    // Constraints. Linear constraints are added all at once, quadratic
    // constraints one by one.
    std::vector<double> rows_rhs;
    std::vector<char> rows_senses;
    std::vector<double> rows_ranges;
    std::vector<int> rows_starts;
    std::vector<int> rows_variables;
    std::vector<double> rows_coefficients;
    std::vector<std::string> rows_names;
    for (int constraint_id = 0; constraint_id < model.number_of_constraints(); ++constraint_id) {
        double lower_bound = model.constraints_lower_bounds[constraint_id];
        double upper_bound = model.constraints_upper_bounds[constraint_id];
        bool has_lower_bound = (lower_bound > -CPX_INFBOUND);
        bool has_upper_bound = (upper_bound < CPX_INFBOUND);
        if (!has_lower_bound && !has_upper_bound)
            continue;

        int linear_start = 0;
        int linear_end = 0;
        if (!model.elements_variables.empty()) {
            linear_start = model.constraints_starts[constraint_id];
            linear_end = model.constraint_end(constraint_id);
        }
        int quadratic_start = 0;
        int quadratic_end = 0;
        if (!model.quadratic_elements_variables_1.empty()) {
            quadratic_start = model.quadratic_elements_constraints_starts[constraint_id];
            quadratic_end = model.quadratic_constraint_end(constraint_id);
        }

        if (quadratic_end == quadratic_start) {
            rows_starts.push_back((int)rows_variables.size());
            for (int element_id = linear_start; element_id < linear_end; ++element_id) {
                rows_variables.push_back(model.elements_variables[element_id]);
                rows_coefficients.push_back(model.elements_coefficients[element_id]);
            }
            // A range row is rhs <= a x <= rhs + range.
            double range = 0.0;
            if (lower_bound == upper_bound) {
                rows_senses.push_back('E');
                rows_rhs.push_back(lower_bound);
            } else if (has_lower_bound && has_upper_bound) {
                rows_senses.push_back('R');
                rows_rhs.push_back(lower_bound);
                range = upper_bound - lower_bound;
            } else if (has_lower_bound) {
                rows_senses.push_back('G');
                rows_rhs.push_back(lower_bound);
            } else {
                rows_senses.push_back('L');
                rows_rhs.push_back(upper_bound);
            }
            rows_ranges.push_back(range);
            rows_names.push_back(model.constraint_name(constraint_id));
            continue;
        }

        // Quadratic constraints have a single sense: equality and range
        // constraints are split in two. Unlike in the objective, there is no
        // 1/2 factor.
        std::vector<int> linear_variables(
                model.elements_variables.begin() + linear_start,
                model.elements_variables.begin() + linear_end);
        std::vector<double> linear_coefficients(
                model.elements_coefficients.begin() + linear_start,
                model.elements_coefficients.begin() + linear_end);
        std::vector<int> quadratic_variables_1(
                model.quadratic_elements_variables_1.begin() + quadratic_start,
                model.quadratic_elements_variables_1.begin() + quadratic_end);
        std::vector<int> quadratic_variables_2(
                model.quadratic_elements_variables_2.begin() + quadratic_start,
                model.quadratic_elements_variables_2.begin() + quadratic_end);
        std::vector<double> quadratic_coefficients(
                model.quadratic_elements_coefficients.begin() + quadratic_start,
                model.quadratic_elements_coefficients.begin() + quadratic_end);
        std::string name = model.constraint_name(constraint_id);
        auto add_quadratic_constraint = [&](char sense, double rhs, const std::string& qc_name)
        {
            cplex_check(
                    env,
                    CPXaddqconstr(
                            env, lp,
                            (int)linear_variables.size(),
                            (int)quadratic_coefficients.size(),
                            rhs,
                            sense,
                            linear_variables.data(),
                            linear_coefficients.data(),
                            quadratic_variables_1.data(),
                            quadratic_variables_2.data(),
                            quadratic_coefficients.data(),
                            qc_name.c_str()),
                    FUNC_SIGNATURE, "CPXaddqconstr");
        };
        if (has_lower_bound && has_upper_bound) {
            add_quadratic_constraint('G', lower_bound, name + "_lower");
            add_quadratic_constraint('L', upper_bound, name + "_upper");
        } else {
            add_quadratic_constraint(
                    has_lower_bound? 'G': 'L',
                    has_lower_bound? lower_bound: upper_bound,
                    name);
        }
    }
    if (!rows_rhs.empty()) {
        std::vector<char*> rows_names_ptrs;
        for (const std::string& name: rows_names)
            rows_names_ptrs.push_back(const_cast<char*>(name.c_str()));
        cplex_check(
                env,
                CPXaddrows(
                        env, lp,
                        0,
                        (int)rows_rhs.size(),
                        (int)rows_variables.size(),
                        rows_rhs.data(),
                        rows_senses.data(),
                        rows_starts.data(),
                        rows_variables.data(),
                        rows_coefficients.data(),
                        NULL,
                        rows_names_ptrs.data()),
                FUNC_SIGNATURE, "CPXaddrows");
        std::vector<int> range_rows;
        std::vector<double> range_values;
        for (int row_id = 0; row_id < (int)rows_senses.size(); ++row_id) {
            if (rows_senses[row_id] == 'R') {
                range_rows.push_back(row_id);
                range_values.push_back(rows_ranges[row_id]);
            }
        }
        if (!range_rows.empty()) {
            cplex_check(
                    env,
                    CPXchgrngval(
                            env, lp,
                            (int)range_rows.size(),
                            range_rows.data(),
                            range_values.data()),
                    FUNC_SIGNATURE, "CPXchgrngval");
        }
    }

    // Initial solution.
    if (!model.variables_initial_values.empty()) {
        if (is_mip(cplex_model)) {
            std::vector<int> begins = {0};
            std::vector<int> variables(number_of_variables);
            for (int variable_id = 0; variable_id < number_of_variables; ++variable_id)
                variables[variable_id] = variable_id;
            std::vector<int> effort_levels = {CPX_MIPSTART_AUTO};
            cplex_check(
                    env,
                    CPXaddmipstarts(
                            env, lp,
                            1,
                            number_of_variables,
                            begins.data(),
                            variables.data(),
                            model.variables_initial_values.data(),
                            effort_levels.data(),
                            NULL),
                    FUNC_SIGNATURE, "CPXaddmipstarts");
        } else {
            cplex_check(
                    env,
                    CPXcopystart(
                            env, lp,
                            NULL, NULL,
                            model.variables_initial_values.data(),
                            NULL, NULL, NULL),
                    FUNC_SIGNATURE, "CPXcopystart");
        }
    }
}

void mathoptsolverscmake::solve(
        CplexModel& cplex_model)
{
    CPXENVptr env = cplex_model.env;
    CPXLPptr lp = cplex_model.lp;
    switch (CPXgetprobtype(env, lp)) {
    case CPXPROB_LP: {
        cplex_check(env, CPXlpopt(env, lp), FUNC_SIGNATURE, "CPXlpopt");
        break;
    } case CPXPROB_QP: {
        cplex_check(env, CPXqpopt(env, lp), FUNC_SIGNATURE, "CPXqpopt");
        break;
    } case CPXPROB_QCP: {
        cplex_check(env, CPXbaropt(env, lp), FUNC_SIGNATURE, "CPXbaropt");
        break;
    } default: {
        cplex_check(env, CPXmipopt(env, lp), FUNC_SIGNATURE, "CPXmipopt");
        break;
    }
    }
}

static bool has_solution(const CplexModel& cplex_model)
{
    int primal_feasible = 0;
    cplex_check(
            cplex_model.env,
            CPXsolninfo(cplex_model.env, cplex_model.lp, NULL, NULL, &primal_feasible, NULL),
            FUNC_SIGNATURE, "CPXsolninfo");
    return primal_feasible != 0;
}

static bool is_minimization(const CplexModel& cplex_model)
{
    return CPXgetobjsen(cplex_model.env, cplex_model.lp) == CPX_MIN;
}

double mathoptsolverscmake::get_solution_value(
        const CplexModel& cplex_model)
{
    if (!has_solution(cplex_model)) {
        return is_minimization(cplex_model)?
            std::numeric_limits<double>::infinity():
            -std::numeric_limits<double>::infinity();
    }
    double objective_value = 0.0;
    cplex_check(
            cplex_model.env,
            CPXgetobjval(cplex_model.env, cplex_model.lp, &objective_value),
            FUNC_SIGNATURE, "CPXgetobjval");
    return objective_value;
}

std::vector<double> mathoptsolverscmake::get_solution(
        const CplexModel& cplex_model)
{
    if (!has_solution(cplex_model))
        return {};
    int number_of_variables = CPXgetnumcols(cplex_model.env, cplex_model.lp);
    std::vector<double> solution(number_of_variables);
    cplex_check(
            cplex_model.env,
            CPXgetx(cplex_model.env, cplex_model.lp, solution.data(), 0, number_of_variables - 1),
            FUNC_SIGNATURE, "CPXgetx");
    return solution;
}

double mathoptsolverscmake::get_bound(
        const CplexModel& cplex_model)
{
    double bound = 0.0;
    if (is_mip(cplex_model)
            && CPXgetbestobjval(cplex_model.env, cplex_model.lp, &bound) == 0) {
        return bound;
    }
    // For continuous models, the bound is the objective value if the model
    // has been solved to optimality.
    if (CPXgetstat(cplex_model.env, cplex_model.lp) == CPX_STAT_OPTIMAL)
        return get_solution_value(cplex_model);
    return is_minimization(cplex_model)?
        -std::numeric_limits<double>::infinity():
        std::numeric_limits<double>::infinity();
}

int mathoptsolverscmake::get_number_of_nodes(
        const CplexModel& cplex_model)
{
    if (!is_mip(cplex_model))
        return 0;
    return (int)CPXgetnodecnt(cplex_model.env, cplex_model.lp);
}

bool mathoptsolverscmake::is_infeasible(
        const CplexModel& cplex_model)
{
    int status = CPXgetstat(cplex_model.env, cplex_model.lp);
    return status == CPX_STAT_INFEASIBLE
        || status == CPX_STAT_INForUNBD
        || status == CPXMIP_INFEASIBLE
        || status == CPXMIP_INForUNBD;
}

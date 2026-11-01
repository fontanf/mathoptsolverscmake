#include "mathoptsolverscmake/mathopt_gurobi.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace mathoptsolverscmake;

static void gurobi_check(GRBenv* env, int error, const std::string& caller, const char* fn)
{
    if (error != 0) {
        throw std::runtime_error(
                caller + ": "
                + fn + " failed with error code "
                + std::to_string(error) + ": "
                + ((env != nullptr)? GRBgeterrormsg(env): "") + ".");
    }
}

static double to_gurobi(double value)
{
    if (value >= GRB_INFINITY)
        return GRB_INFINITY;
    if (value <= -GRB_INFINITY)
        return -GRB_INFINITY;
    return value;
}

/** Environment in which to set parameters: the one of the model once loaded. */
static GRBenv* parameters_env(const GurobiModel& gurobi_model)
{
    if (gurobi_model.model != nullptr)
        return GRBgetenv(gurobi_model.model);
    return gurobi_model.env;
}

/**
 * Add a nonlinear expression tree (stored as parallel element arrays, see
 * mathopt.hpp) as a general constraint resvar = f(x).
 *
 * Gurobi takes the nodes of the tree in an order in which every node comes
 * after its parent, each node given by its parent's index. The DFS pre-order
 * of MathOptModel is such an order, so only the parents have to be computed.
 * Gurobi supports all the operators of MathOptModel.
 */
static void add_nonlinear_constraint(
        GRBmodel* model,
        const std::string& name,
        int resvar,
        const std::vector<char>& operators,
        const std::vector<double>& values,
        const std::vector<int>& variables,
        const std::vector<int>& number_of_children,
        int start,
        int end)
{
    std::vector<int> opcodes;
    std::vector<double> data;
    std::vector<int> parents;
    // Nodes whose children have not all been visited yet, with their number
    // of remaining children.
    std::vector<std::pair<int, int>> stack;
    for (int element_id = start; element_id < end; ++element_id) {
        int parent = -1;
        if (!stack.empty()) {
            parent = stack.back().first;
            if (--stack.back().second == 0)
                stack.pop_back();
        }
        int node_id = element_id - start;
        parents.push_back(parent);
        if (number_of_children[element_id] > 0)
            stack.push_back({node_id, number_of_children[element_id]});

        int opcode = -1;
        double node_data = -1.0;
        switch (operators[element_id]) {
        case 'k': opcode = GRB_OPCODE_CONSTANT; node_data = values[element_id]; break;
        case 'v': opcode = GRB_OPCODE_VARIABLE; node_data = variables[element_id]; break;
        case '+': opcode = GRB_OPCODE_PLUS; break;
        case '-': opcode = GRB_OPCODE_MINUS; break;
        case '*': opcode = GRB_OPCODE_MULTIPLY; break;
        case '/': opcode = GRB_OPCODE_DIVIDE; break;
        case 'n': opcode = GRB_OPCODE_UMINUS; break;
        case 'e': opcode = GRB_OPCODE_EXP; break;
        case 'l': opcode = GRB_OPCODE_LOG; break;
        case 'q': opcode = GRB_OPCODE_SQRT; break;
        case 's': opcode = GRB_OPCODE_SIN; break;
        case 'c': opcode = GRB_OPCODE_COS; break;
        case 't': opcode = GRB_OPCODE_TAN; break;
        case 'p': opcode = GRB_OPCODE_POW; break;
        default: {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "unknown nonlinear operator '"
                    + std::string(1, operators[element_id]) + "'.");
        }
        }
        opcodes.push_back(opcode);
        data.push_back(node_data);
    }
    gurobi_check(
            GRBgetenv(model),
            GRBaddgenconstrNL(
                    model,
                    name.c_str(),
                    resvar,
                    (int)opcodes.size(),
                    opcodes.data(),
                    data.data(),
                    parents.data()),
            FUNC_SIGNATURE, "GRBaddgenconstrNL");
}

GurobiModel::GurobiModel()
{
    gurobi_check(nullptr, GRBemptyenv(&env), FUNC_SIGNATURE, "GRBemptyenv");
}

GurobiModel::~GurobiModel()
{
    if (model != nullptr)
        GRBfreemodel(model);
    if (env != nullptr)
        GRBfreeenv(env);
}

void mathoptsolverscmake::reduce_printout(
        GurobiModel& gurobi_model)
{
    GRBenv* env = parameters_env(gurobi_model);
    gurobi_check(
            env,
            GRBsetintparam(env, GRB_INT_PAR_LOGTOCONSOLE, 0),
            FUNC_SIGNATURE, "GRBsetintparam");
}

void mathoptsolverscmake::set_time_limit(
        GurobiModel& gurobi_model,
        double time_limit)
{
    GRBenv* env = parameters_env(gurobi_model);
    gurobi_check(
            env,
            GRBsetdblparam(env, GRB_DBL_PAR_TIMELIMIT, time_limit),
            FUNC_SIGNATURE, "GRBsetdblparam");
}

void mathoptsolverscmake::set_node_limit(
        GurobiModel& gurobi_model,
        int node_limit)
{
    GRBenv* env = parameters_env(gurobi_model);
    gurobi_check(
            env,
            GRBsetdblparam(env, GRB_DBL_PAR_NODELIMIT, node_limit),
            FUNC_SIGNATURE, "GRBsetdblparam");
}

void mathoptsolverscmake::set_log_file(
        GurobiModel& gurobi_model,
        const std::string& log_file)
{
    GRBenv* env = parameters_env(gurobi_model);
    gurobi_check(
            env,
            GRBsetstrparam(env, GRB_STR_PAR_LOGFILE, log_file.c_str()),
            FUNC_SIGNATURE, "GRBsetstrparam");
}

void mathoptsolverscmake::load(
        GurobiModel& gurobi_model,
        const MathOptModel& model)
{
    if (model.has_black_box()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "Gurobi does not support black-box functions.");
    }
    if (gurobi_model.model != nullptr) {
        throw std::logic_error(
                FUNC_SIGNATURE + ": "
                "a model has already been loaded.");
    }

    GRBenv* env = gurobi_model.env;
    gurobi_check(env, GRBstartenv(env), FUNC_SIGNATURE, "GRBstartenv");
    gurobi_check(
            env,
            GRBnewmodel(env, &gurobi_model.model, "mathopt", 0, NULL, NULL, NULL, NULL, NULL),
            FUNC_SIGNATURE, "GRBnewmodel");
    GRBmodel* grb_model = gurobi_model.model;
    env = GRBgetenv(grb_model);
    gurobi_model.number_of_variables = model.number_of_variables();

    gurobi_check(
            env,
            GRBsetintattr(
                    grb_model,
                    GRB_INT_ATTR_MODELSENSE,
                    (model.objective_direction == ObjectiveDirection::Minimize)?
                    GRB_MINIMIZE:
                    GRB_MAXIMIZE),
            FUNC_SIGNATURE, "GRBsetintattr");

    // Variables of the model, with the linear objective.
    {
        int number_of_variables = model.number_of_variables();
        std::vector<double> objective(number_of_variables, 0.0);
        std::vector<double> lower_bounds(number_of_variables);
        std::vector<double> upper_bounds(number_of_variables);
        std::vector<char> types(number_of_variables);
        std::vector<std::string> names(number_of_variables);
        std::vector<char*> names_ptrs(number_of_variables);
        for (int variable_id = 0; variable_id < number_of_variables; ++variable_id) {
            if (!model.objective_coefficients.empty())
                objective[variable_id] = model.objective_coefficients[variable_id];
            lower_bounds[variable_id] = to_gurobi(model.variables_lower_bounds[variable_id]);
            upper_bounds[variable_id] = to_gurobi(model.variables_upper_bounds[variable_id]);
            switch (model.variables_types[variable_id]) {
            case VariableType::Continuous: types[variable_id] = GRB_CONTINUOUS; break;
            case VariableType::Binary: types[variable_id] = GRB_BINARY; break;
            case VariableType::Integer: types[variable_id] = GRB_INTEGER; break;
            }
            names[variable_id] = model.variable_name(variable_id);
            names_ptrs[variable_id] = const_cast<char*>(names[variable_id].c_str());
        }
        gurobi_check(
                env,
                GRBaddvars(
                        grb_model,
                        number_of_variables,
                        0, NULL, NULL, NULL,
                        objective.data(),
                        lower_bounds.data(),
                        upper_bounds.data(),
                        types.data(),
                        names_ptrs.data()),
                FUNC_SIGNATURE, "GRBaddvars");
    }

    // Auxiliary variables of the nonlinear parts of the objective and of the
    // constraints. They are added before the constraints since
    // GRBaddrangeconstr also adds variables.
    int next_variable_id = model.number_of_variables();
    auto add_auxiliary_variable = [env, grb_model, &next_variable_id](
            const std::string& name,
            double objective_coefficient)
    {
        gurobi_check(
                env,
                GRBaddvar(
                        grb_model,
                        0, NULL, NULL,
                        objective_coefficient,
                        -GRB_INFINITY,
                        GRB_INFINITY,
                        GRB_CONTINUOUS,
                        name.c_str()),
                FUNC_SIGNATURE, "GRBaddvar");
        return next_variable_id++;
    };
    if (!model.objective_nonlinear_elements_operators.empty()) {
        int objective_variable_id = add_auxiliary_variable("objective", 1.0);
        add_nonlinear_constraint(
                grb_model,
                "objective",
                objective_variable_id,
                model.objective_nonlinear_elements_operators,
                model.objective_nonlinear_elements_values,
                model.objective_nonlinear_elements_variables,
                model.objective_nonlinear_elements_number_of_children,
                0,
                (int)model.objective_nonlinear_elements_operators.size());
    }
    std::vector<int> nonlinear_variables(model.number_of_constraints(), -1);
    if (!model.nonlinear_elements_operators.empty()) {
        for (int constraint_id = 0; constraint_id < model.number_of_constraints(); ++constraint_id) {
            int nonlinear_start = model.nonlinear_elements_constraints_starts[constraint_id];
            int nonlinear_end = model.nonlinear_constraint_end(constraint_id);
            if (nonlinear_end == nonlinear_start)
                continue;
            std::string name = model.constraint_name(constraint_id) + "_nonlinear";
            nonlinear_variables[constraint_id] = add_auxiliary_variable(name, 0.0);
            add_nonlinear_constraint(
                    grb_model,
                    name,
                    nonlinear_variables[constraint_id],
                    model.nonlinear_elements_operators,
                    model.nonlinear_elements_values,
                    model.nonlinear_elements_variables,
                    model.nonlinear_elements_number_of_children,
                    nonlinear_start,
                    nonlinear_end);
        }
    }

    // Quadratic objective.
    if (!model.objective_quadratic_elements_variables_1.empty()) {
        std::vector<int> rows = model.objective_quadratic_elements_variables_1;
        std::vector<int> columns = model.objective_quadratic_elements_variables_2;
        std::vector<double> coefficients = model.objective_quadratic_elements_coefficients;
        gurobi_check(
                env,
                GRBaddqpterms(
                        grb_model,
                        (int)coefficients.size(),
                        rows.data(),
                        columns.data(),
                        coefficients.data()),
                FUNC_SIGNATURE, "GRBaddqpterms");
    }

    // Constraints.
    for (int constraint_id = 0; constraint_id < model.number_of_constraints(); ++constraint_id) {
        double lower_bound = model.constraints_lower_bounds[constraint_id];
        double upper_bound = model.constraints_upper_bounds[constraint_id];
        bool has_lower_bound = (lower_bound > -GRB_INFINITY);
        bool has_upper_bound = (upper_bound < GRB_INFINITY);
        if (!has_lower_bound && !has_upper_bound)
            continue;
        std::string name = model.constraint_name(constraint_id);

        std::vector<int> linear_variables;
        std::vector<double> linear_coefficients;
        if (!model.elements_variables.empty()) {
            for (int element_id = model.constraints_starts[constraint_id];
                    element_id < model.constraint_end(constraint_id);
                    ++element_id) {
                linear_variables.push_back(model.elements_variables[element_id]);
                linear_coefficients.push_back(model.elements_coefficients[element_id]);
            }
        }
        if (nonlinear_variables[constraint_id] != -1) {
            linear_variables.push_back(nonlinear_variables[constraint_id]);
            linear_coefficients.push_back(1.0);
        }

        std::vector<int> quadratic_variables_1;
        std::vector<int> quadratic_variables_2;
        std::vector<double> quadratic_coefficients;
        if (!model.quadratic_elements_variables_1.empty()) {
            for (int element_id = model.quadratic_elements_constraints_starts[constraint_id];
                    element_id < model.quadratic_constraint_end(constraint_id);
                    ++element_id) {
                quadratic_variables_1.push_back(model.quadratic_elements_variables_1[element_id]);
                quadratic_variables_2.push_back(model.quadratic_elements_variables_2[element_id]);
                quadratic_coefficients.push_back(model.quadratic_elements_coefficients[element_id]);
            }
        }

        if (quadratic_coefficients.empty()) {
            if (lower_bound == upper_bound) {
                gurobi_check(
                        env,
                        GRBaddconstr(
                                grb_model,
                                (int)linear_variables.size(),
                                linear_variables.data(),
                                linear_coefficients.data(),
                                GRB_EQUAL,
                                lower_bound,
                                name.c_str()),
                        FUNC_SIGNATURE, "GRBaddconstr");
            } else if (has_lower_bound && has_upper_bound) {
                gurobi_check(
                        env,
                        GRBaddrangeconstr(
                                grb_model,
                                (int)linear_variables.size(),
                                linear_variables.data(),
                                linear_coefficients.data(),
                                lower_bound,
                                upper_bound,
                                name.c_str()),
                        FUNC_SIGNATURE, "GRBaddrangeconstr");
            } else {
                gurobi_check(
                        env,
                        GRBaddconstr(
                                grb_model,
                                (int)linear_variables.size(),
                                linear_variables.data(),
                                linear_coefficients.data(),
                                has_lower_bound? GRB_GREATER_EQUAL: GRB_LESS_EQUAL,
                                has_lower_bound? lower_bound: upper_bound,
                                name.c_str()),
                        FUNC_SIGNATURE, "GRBaddconstr");
            }
            continue;
        }

        // Quadratic constraints have a single sense: a range constraint is
        // split in two.
        auto add_quadratic_constraint = [&](char sense, double rhs, const std::string& qc_name)
        {
            gurobi_check(
                    env,
                    GRBaddqconstr(
                            grb_model,
                            (int)linear_variables.size(),
                            linear_variables.data(),
                            linear_coefficients.data(),
                            (int)quadratic_coefficients.size(),
                            quadratic_variables_1.data(),
                            quadratic_variables_2.data(),
                            quadratic_coefficients.data(),
                            sense,
                            rhs,
                            qc_name.c_str()),
                    FUNC_SIGNATURE, "GRBaddqconstr");
        };
        if (lower_bound == upper_bound) {
            add_quadratic_constraint(GRB_EQUAL, lower_bound, name);
        } else if (has_lower_bound && has_upper_bound) {
            add_quadratic_constraint(GRB_GREATER_EQUAL, lower_bound, name + "_lower");
            add_quadratic_constraint(GRB_LESS_EQUAL, upper_bound, name + "_upper");
        } else {
            add_quadratic_constraint(
                    has_lower_bound? GRB_GREATER_EQUAL: GRB_LESS_EQUAL,
                    has_lower_bound? lower_bound: upper_bound,
                    name);
        }
    }

    // Initial solution. The values of the auxiliary variables are left for
    // Gurobi to complete.
    if (!model.variables_initial_values.empty()) {
        std::vector<double> values = model.variables_initial_values;
        gurobi_check(
                env,
                GRBsetdblattrarray(
                        grb_model,
                        GRB_DBL_ATTR_START,
                        0,
                        (int)values.size(),
                        values.data()),
                FUNC_SIGNATURE, "GRBsetdblattrarray");
    }

    gurobi_check(env, GRBupdatemodel(grb_model), FUNC_SIGNATURE, "GRBupdatemodel");
}

void mathoptsolverscmake::solve(
        GurobiModel& gurobi_model)
{
    gurobi_check(
            GRBgetenv(gurobi_model.model),
            GRBoptimize(gurobi_model.model),
            FUNC_SIGNATURE, "GRBoptimize");
}

static int number_of_solutions(const GurobiModel& gurobi_model)
{
    int solution_count = 0;
    gurobi_check(
            GRBgetenv(gurobi_model.model),
            GRBgetintattr(gurobi_model.model, GRB_INT_ATTR_SOLCOUNT, &solution_count),
            FUNC_SIGNATURE, "GRBgetintattr");
    return solution_count;
}

static bool is_minimization(const GurobiModel& gurobi_model)
{
    int model_sense = GRB_MINIMIZE;
    gurobi_check(
            GRBgetenv(gurobi_model.model),
            GRBgetintattr(gurobi_model.model, GRB_INT_ATTR_MODELSENSE, &model_sense),
            FUNC_SIGNATURE, "GRBgetintattr");
    return model_sense == GRB_MINIMIZE;
}

double mathoptsolverscmake::get_solution_value(
        const GurobiModel& gurobi_model)
{
    if (number_of_solutions(gurobi_model) == 0) {
        return is_minimization(gurobi_model)?
            std::numeric_limits<double>::infinity():
            -std::numeric_limits<double>::infinity();
    }
    double objective_value = 0.0;
    gurobi_check(
            GRBgetenv(gurobi_model.model),
            GRBgetdblattr(gurobi_model.model, GRB_DBL_ATTR_OBJVAL, &objective_value),
            FUNC_SIGNATURE, "GRBgetdblattr");
    return objective_value;
}

std::vector<double> mathoptsolverscmake::get_solution(
        const GurobiModel& gurobi_model)
{
    if (number_of_solutions(gurobi_model) == 0)
        return {};
    std::vector<double> solution(gurobi_model.number_of_variables);
    gurobi_check(
            GRBgetenv(gurobi_model.model),
            GRBgetdblattrarray(
                    gurobi_model.model,
                    GRB_DBL_ATTR_X,
                    0,
                    gurobi_model.number_of_variables,
                    solution.data()),
            FUNC_SIGNATURE, "GRBgetdblattrarray");
    return solution;
}

double mathoptsolverscmake::get_bound(
        const GurobiModel& gurobi_model)
{
    double bound = 0.0;
    if (GRBgetdblattr(gurobi_model.model, GRB_DBL_ATTR_OBJBOUND, &bound) == 0)
        return bound;
    // The bound is not available for continuous models: it is the objective
    // value if the model has been solved to optimality.
    int status = 0;
    gurobi_check(
            GRBgetenv(gurobi_model.model),
            GRBgetintattr(gurobi_model.model, GRB_INT_ATTR_STATUS, &status),
            FUNC_SIGNATURE, "GRBgetintattr");
    if (status == GRB_OPTIMAL)
        return get_solution_value(gurobi_model);
    return is_minimization(gurobi_model)?
        -std::numeric_limits<double>::infinity():
        std::numeric_limits<double>::infinity();
}

int mathoptsolverscmake::get_number_of_nodes(
        const GurobiModel& gurobi_model)
{
    double number_of_nodes = 0.0;
    // Not available for continuous models.
    if (GRBgetdblattr(gurobi_model.model, GRB_DBL_ATTR_NODECOUNT, &number_of_nodes) != 0)
        return 0;
    return (int)number_of_nodes;
}

bool mathoptsolverscmake::is_infeasible(
        const GurobiModel& gurobi_model)
{
    int status = 0;
    gurobi_check(
            GRBgetenv(gurobi_model.model),
            GRBgetintattr(gurobi_model.model, GRB_INT_ATTR_STATUS, &status),
            FUNC_SIGNATURE, "GRBgetintattr");
    return status == GRB_INFEASIBLE
        || status == GRB_INF_OR_UNBD;
}

#include "mathoptsolverscmake/mathopt_xpress.hpp"

#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

using namespace mathoptsolverscmake;

void mathoptsolverscmake::load(
        XPRSprob& xpress_model,
        const MathOptModel& model)
{
    if (!model.is_milp()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "model is not an MILP; Xpress only supports MILP models.");
    }

    int status = 0;

    // Add variables.
    status = XPRSaddcols(
            xpress_model,
            model.number_of_variables(),
            0,  // ncoeffs
            model.objective_coefficients.data(),
            NULL,  // starts
            NULL,  // constraintind
            NULL,  // constraintcoef
            model.variables_lower_bounds.data(),
            model.variables_upper_bounds.data());
    if (status) {
        throw std::runtime_error(
                "mathoptsolverscmake::load: "
                "error loading variables; "
                "status: " + std::to_string(status) + ".");
    }

    // Add constraints.
    std::vector<char> constraints_types(model.number_of_constraints());
    std::vector<double> rhs(model.number_of_constraints());
    std::vector<double> rng(model.number_of_constraints());
    for (int constraint_id = 0;
            constraint_id < model.number_of_constraints();
            ++constraint_id) {
        if (model.constraints_lower_bounds[constraint_id] == -std::numeric_limits<double>::infinity()) {
            constraints_types[constraint_id] = 'L';
            rhs[constraint_id] = model.constraints_upper_bounds[constraint_id];
        } else if (model.constraints_upper_bounds[constraint_id] == std::numeric_limits<double>::infinity()) {
            constraints_types[constraint_id] = 'G';
            rhs[constraint_id] = model.constraints_lower_bounds[constraint_id];
        } else {
            constraints_types[constraint_id] = 'R';
            rhs[constraint_id] = model.constraints_upper_bounds[constraint_id];
            rng[constraint_id] = model.constraints_upper_bounds[constraint_id] - model.constraints_lower_bounds[constraint_id];
        }
    }
    status = XPRSaddrows(
            xpress_model,
            model.number_of_constraints(),
            model.number_of_elements(),
            constraints_types.data(),
            rhs.data(),
            rng.data(),
            model.constraints_starts.data(),
            model.elements_variables.data(),
            model.elements_coefficients.data());
    if (status) {
        throw std::runtime_error(
                "mathoptsolverscmake::load: "
                "error loading constraints; "
                "status: " + std::to_string(status) + ".");
    }

    // Set objective sense.
    if (model.objective_direction == ObjectiveDirection::Minimize) {
        XPRSchgobjsense(xpress_model, XPRS_OBJ_MINIMIZE);
    } else {
        XPRSchgobjsense(xpress_model, XPRS_OBJ_MAXIMIZE);
    }

    // Set variable types
    std::vector<int> variables_indices(model.number_of_variables(), 0);
    std::iota(variables_indices.begin(), variables_indices.end(), 0);
    std::vector<char> variables_types(model.number_of_variables(), 'A');
    for (int variable_id = 0;
            variable_id < model.number_of_variables();
            ++variable_id) {
        switch (model.variables_types[variable_id]) {
        case VariableType::Continuous:
            variables_types[variable_id] = 'B';
            break;
        case VariableType::Binary:
            variables_types[variable_id] = 'B';
            break;
        case VariableType::Integer:
            variables_types[variable_id] = 'I';
            break;
        }
    }
    status = XPRSchgcoltype(
            xpress_model,
            variables_indices.size(),
            variables_indices.data(),
            variables_types.data());
    if (status) {
        throw std::runtime_error(
                "mathoptsolverscmake::load: "
                "error setting variable types; "
                "status: " + std::to_string(status) + ".");
    }
}

void mathoptsolverscmake::set_time_limit(
        XPRSprob& xpress_model,
        double time_limit)
{
    XPRSsetdblcontrol(xpress_model, XPRS_TIMELIMIT, time_limit);
}

void mathoptsolverscmake::set_log_file(
        XPRSprob& xpress_model,
        const std::string& log_file)
{
    XPRSsetlogfile(xpress_model, log_file.c_str());
}

void mathoptsolverscmake::write_mps(
        const XPRSprob& xpress_model,
        const std::string& mps_file)
{
    XPRSwriteprob(xpress_model, mps_file.c_str(), "ot");
}

void mathoptsolverscmake::solve(
        XPRSprob& xpress_model)
{
    int status = XPRSmipoptimize(xpress_model, "");
    if (status) {
        throw std::runtime_error(
                "mathoptsolverscmake::solve: "
                "error solving model; "
                "status: " + std::to_string(status) + ".");
    }
}

double mathoptsolverscmake::get_solution_value(
        const XPRSprob& xpress_model)
{
    double objective_value = 0;
    XPRSgetdblattrib(xpress_model, XPRS_MIPOBJVAL, &objective_value);
    return objective_value;
}

std::vector<double> mathoptsolverscmake::get_solution(
        const XPRSprob& xpress_model)
{
    int number_of_variables = 0;
    XPRSgetintattrib(xpress_model, XPRS_ORIGINALCOLS, &number_of_variables);
    std::vector<double> solution(number_of_variables);
    int status = XPRSgetmipsol(
            xpress_model,
            solution.data(),
            NULL);
    if (status) {
        throw std::runtime_error(
                "mathoptsolverscmake::get_solution: "
                "error retrieving solution; "
                "status: " + std::to_string(status) + ".");
    }
    return solution;
}

double mathoptsolverscmake::get_bound(
        const XPRSprob& xpress_model)
{
    double bound = 0.0;
    XPRSgetdblattrib(xpress_model, XPRS_BESTBOUND, &bound);
    return bound;
}

int mathoptsolverscmake::get_number_of_nodes(
        const XPRSprob& xpress_model)
{
    int number_of_nodes = 0;
    XPRSgetintattrib(xpress_model, XPRS_NODES, &number_of_nodes);
    return number_of_nodes;
}

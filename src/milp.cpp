#include "mathoptsolverscmake/milp.hpp"

#include <cstddef>
#include <stdexcept>
#include <numeric>

using namespace mathoptsolverscmake;

#ifdef CBC_FOUND

void mathoptsolverscmake::load(
        CbcModel& cbc_model,
        const MilpModel& model)
{
    std::vector<int> constraints_lengths(model.number_of_constraints(), 0);
    for (int constraint_id = 0;
            constraint_id < model.number_of_constraints() - 1;
            ++constraint_id) {
        constraints_lengths[constraint_id] = model.constraints_starts[constraint_id + 1] - model.constraints_starts[constraint_id];
    }
    constraints_lengths[model.number_of_constraints() - 1] = model.number_of_elements() - model.constraints_starts[model.number_of_constraints() - 1];
    CoinPackedMatrix matrix(
            false,
            model.number_of_variables(),
            model.number_of_constraints(),
            model.number_of_elements(),
            model.elements_coefficients.data(),
            model.elements_variables.data(),
            model.constraints_starts.data(),
            constraints_lengths.data());
    cbc_model.solver()->loadProblem(
            matrix,
            model.variables_lower_bounds.data(),
            model.variables_upper_bounds.data(),
            model.objective_coefficients.data(),
            model.constraints_lower_bounds.data(),
            model.constraints_upper_bounds.data());

    // Set objective direction.
    if (model.objective_direction == ObjectiveDirection::Minimize) {
        cbc_model.setObjSense(1);
    } else {
        cbc_model.setObjSense(-1);
    }

    // Set variable types.
    std::vector<int> integer_variables;
    for (int variable_id = 0;
            variable_id < model.number_of_variables();
            ++variable_id) {
        if (model.variables_types[variable_id] != VariableType::Continuous)
            integer_variables.push_back(variable_id);
    }
    cbc_model.solver()->setInteger(integer_variables.data(), integer_variables.size());
}

void mathoptsolverscmake::reduce_printout(
        CbcModel& cbc_model)
{
    cbc_model.setLogLevel(0);
    cbc_model.messageHandler()->setLogLevel(0);
    cbc_model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
}

void mathoptsolverscmake::set_time_limit(
        CbcModel& cbc_model,
        double time_limit)
{
    cbc_model.setMaximumSeconds(time_limit);
}

void mathoptsolverscmake::solve(
        CbcModel& cbc_model)
{
    cbc_model.branchAndBound();
}

double mathoptsolverscmake::get_solution_value(
        const CbcModel& cbc_model)
{
    if (cbc_model.bestSolution() == nullptr)
        return cbc_model.getObjSense() * std::numeric_limits<double>::infinity();
    return cbc_model.getObjValue();
}

std::vector<double> mathoptsolverscmake::get_solution(
        const CbcModel& cbc_model)
{
    return std::vector<double>(
            cbc_model.bestSolution(),
            cbc_model.bestSolution() + cbc_model.getNumCols());
}

double mathoptsolverscmake::get_bound(
        const CbcModel& cbc_model)
{
    return cbc_model.getBestPossibleObjValue();
}

int mathoptsolverscmake::get_number_of_nodes(
        const CbcModel& cbc_model)
{
    return cbc_model.getNodeCount();
}

#endif

#ifdef HIGHS_FOUND

void mathoptsolverscmake::load(
        Highs& highs_model,
        const MilpModel& model)
{
    HighsStatus highs_status;
    highs_status = highs_model.addCols(
            model.number_of_variables(),
            model.objective_coefficients.data(),
            model.variables_lower_bounds.data(),
            model.variables_upper_bounds.data(),
            0,
            NULL,
            NULL,
            NULL);
    highs_status = highs_model.addRows(
            model.number_of_constraints(),
            model.constraints_lower_bounds.data(),
            model.constraints_upper_bounds.data(),
            model.elements_coefficients.size(),
            model.constraints_starts.data(),
            model.elements_variables.data(),
            model.elements_coefficients.data());

    // Set objective sense.
    if (model.objective_direction == ObjectiveDirection::Minimize) {
        highs_model.changeObjectiveSense(ObjSense::kMinimize);
    } else {
        highs_model.changeObjectiveSense(ObjSense::kMaximize);
    }

    // Set variable types.
    std::vector<int> variables_indices(model.number_of_variables(), 0);
    std::iota(variables_indices.begin(), variables_indices.end(), 0);
    std::vector<HighsVarType> variables_types(model.number_of_variables(), HighsVarType::kContinuous);
    for (int variable_id = 0;
            variable_id < model.number_of_variables();
            ++variable_id) {
        switch (model.variables_types[variable_id]) {
        case VariableType::Continuous:
            variables_types[variable_id] = HighsVarType::kContinuous;
            break;
        case VariableType::Binary:
            variables_types[variable_id] = HighsVarType::kInteger;
            break;
        case VariableType::Integer:
            variables_types[variable_id] = HighsVarType::kInteger;
            break;
        }
    }
    highs_model.changeColsIntegrality(
            variables_indices.size(),
            variables_indices.data(),
            variables_types.data());
}

void mathoptsolverscmake::reduce_printout(
        Highs& highs_model)
{
    HighsStatus return_status = highs_model.setOptionValue(
            "log_to_console",
            false);
}

void mathoptsolverscmake::set_time_limit(
        Highs& highs_model,
        double time_limit)
{
    HighsStatus highs_status = highs_model.setOptionValue(
            "time_limit",
            time_limit);
}

void mathoptsolverscmake::set_log_file(
        Highs& highs_model,
        const std::string& log_file)
{
    HighsStatus return_status = highs_model.setOptionValue(
            "log_file",
            log_file);
}

void mathoptsolverscmake::solve(
        Highs& highs_model)
{
    HighsStatus return_status = highs_model.run();
}

std::vector<double> mathoptsolverscmake::get_solution(
        const Highs& highs_model)
{
    return highs_model.getSolution().col_value;
}

double mathoptsolverscmake::get_bound(
        const Highs& highs_model)
{
    return highs_model.getInfo().mip_dual_bound;
}

#endif

#ifdef XPRESS_FOUND

void mathoptsolverscmake::load(
        XPRSprob& xpress_model,
        const MilpModel& model)
{
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

int get_number_of_nodes(
        const XPRSprob& xpress_model)
{
    int number_of_nodes = 0;
    XPRSgetintattrib(xpress_model, XPRS_NODES, &number_of_nodes);
    return number_of_nodes;
}

#endif

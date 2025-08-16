#include "mathoptsolverscmake/milp.hpp"

#include <cstddef>
#include <stdexcept>
#include <numeric>

using namespace mathoptsolverscmake;

#ifdef CBC_FOUND

void mathoptsolverscmake::cbc::load(
        OsiCbcSolverInterface& cbc_model,
        const MilpModel& model)
{
    cbc_model.loadProblem(
            model.number_of_variables(),
            model.number_of_constraints(),
            model.constraints_starts.data(),
            model.elements_variables.data(),
            model.elements_coefficients.data(),
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
    cbc_model.setInteger(integer_variables.data(), integer_variables.size());
}

#endif

#ifdef HIGHS_FOUND

void mathoptsolverscmake::highs::load(
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

#endif

#ifdef XPRESS_FOUND

void mathoptsolverscmake::xpress::load(
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
                "knapsackwithconflictssolver::xpress_milp: "
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
                "knapsackwithconflictssolver::xpress_milp: "
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
                "knapsackwithconflictssolver::xpress_milp: "
                "error setting variable types; "
                "status: " + std::to_string(status) + ".");
    }

}

#endif

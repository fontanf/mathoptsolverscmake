#include "mathoptsolverscmake/mathopt_highs.hpp"

#include <numeric>
#include <stdexcept>
#include <vector>

using namespace mathoptsolverscmake;

void mathoptsolverscmake::load(
        Highs& highs_model,
        const MathOptModel& model)
{
    if (!model.is_milp()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "model is not an MILP; HiGHS only supports MILP models.");
    }

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
        case VariableType::SemiContinuous:
            variables_types[variable_id] = HighsVarType::kSemiContinuous;
            break;
        case VariableType::SemiInteger:
            variables_types[variable_id] = HighsVarType::kSemiInteger;
            break;
        }
    }
    highs_model.changeColsIntegrality(
            variables_indices.size(),
            variables_indices.data(),
            variables_types.data());

    // Set variable names.
    if (!model.variables_names.empty()) {
        for (int variable_id = 0;
                variable_id < model.number_of_variables();
                ++variable_id) {
            highs_model.passColName(variable_id, model.variables_names[variable_id]);
        }
    }

    // Set constraint names.
    if (!model.constraints_names.empty()) {
        for (int constraint_id = 0;
                constraint_id < model.number_of_constraints();
                ++constraint_id) {
            highs_model.passRowName(constraint_id, model.constraints_names[constraint_id]);
        }
    }

    // Set initial values.
    if (!model.variables_initial_values.empty()) {
        HighsSolution highs_solution;
        highs_solution.col_value = model.variables_initial_values;
        highs_solution.value_valid = true;
        highs_model.setSolution(highs_solution);
    }

    // Set tolerances. A tolerance of 0 (the default of the model) means the
    // default of HiGHS; HiGHS rejects 0.
    if (model.feasibility_tolerance > 0) {
        highs_model.setOptionValue(
                "primal_feasibility_tolerance",
                model.feasibility_tolerance);
    }
    if (model.integrality_tolerance > 0) {
        highs_model.setOptionValue(
                "mip_feasibility_tolerance",
                model.integrality_tolerance);
    }
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

void mathoptsolverscmake::set_node_limit(
        Highs& highs_model,
        int node_limit)
{
    HighsStatus highs_status = highs_model.setOptionValue(
            "mip_max_nodes",
            node_limit);
}

void mathoptsolverscmake::set_log_file(
        Highs& highs_model,
        const std::string& log_file)
{
    HighsStatus return_status = highs_model.setOptionValue(
            "log_file",
            log_file);
}

void mathoptsolverscmake::write_mps(
        Highs& highs_model,
        const std::string& mps_file)
{
    HighsStatus return_status = highs_model.writeModel(mps_file);
}

void mathoptsolverscmake::solve(
        Highs& highs_model)
{
    // HiGHS keeps one task scheduler per calling thread, created by the first
    // solve in that thread with the number of threads of its 'threads'
    // option, and kept until the thread exits. A later solve in the same
    // thread with another nonzero 'threads' value returns an error without
    // solving. And when a thread holding a scheduler with worker threads
    // exits, HiGHS detaches the workers, which crashes the process on
    // glibc < 2.43 (https://github.com/ERGO-Code/HiGHS/issues/3826).
    //
    // So, unless the solve uses a single thread (its scheduler has no worker
    // threads), reset the scheduler of the calling thread before the solve,
    // so that it doesn't conflict with a scheduler left by a previous solve,
    // and after it, so that no worker threads are left behind; 'true' waits
    // for the workers to end, joining them rather than detaching them. The
    // solution and the information of the model remain available.
    //
    // This must not be called while another HiGHS solve is in progress in
    // the same thread (for example, from a HiGHS callback).
    //
    // 'resetGlobalScheduler' rather than 'resetThreadScheduler', which
    // replaces it, since the latter doesn't exist in HiGHS 1.13.
    HighsInt threads = 0;
    highs_model.getOptionValue("threads", threads);
    bool reset_scheduler = (threads != 1);
    if (reset_scheduler)
        Highs::resetGlobalScheduler(true);
    HighsStatus return_status = highs_model.run();
    if (reset_scheduler)
        Highs::resetGlobalScheduler(true);
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

int mathoptsolverscmake::get_number_of_nodes(
        const Highs& highs_model)
{
    return highs_model.getInfo().mip_node_count;
}

MathOptModel mathoptsolverscmake::to_mathopt(
        Highs& highs_model)
{
    HighsLp lp = highs_model.getLp();
    lp.ensureRowwise();

    MathOptModel model;

    // Objective direction.
    model.objective_direction = (lp.sense_ == ObjSense::kMinimize)?
        ObjectiveDirection::Minimize:
        ObjectiveDirection::Maximize;

    // Variables.
    model.variables_lower_bounds = lp.col_lower_;
    model.variables_upper_bounds = lp.col_upper_;
    model.objective_coefficients = lp.col_cost_;
    if (!lp.col_names_.empty())
        model.variables_names = lp.col_names_;

    // Variable types.
    model.variables_types.resize(lp.num_col_, VariableType::Continuous);
    for (int variable_id = 0; variable_id < lp.num_col_; ++variable_id) {
        if (lp.integrality_.empty())
            break;
        if (lp.integrality_[variable_id] == HighsVarType::kInteger
                || lp.integrality_[variable_id] == HighsVarType::kImplicitInteger) {
            if (lp.col_lower_[variable_id] == 0.0
                    && lp.col_upper_[variable_id] == 1.0) {
                model.variables_types[variable_id] = VariableType::Binary;
            } else {
                model.variables_types[variable_id] = VariableType::Integer;
            }
        }
    }

    // Constraints.
    model.constraints_lower_bounds = lp.row_lower_;
    model.constraints_upper_bounds = lp.row_upper_;
    if (!lp.row_names_.empty())
        model.constraints_names = lp.row_names_;

    // Constraint matrix (row-wise CSR).
    model.constraints_starts = std::vector<int>(
            lp.a_matrix_.start_.begin(),
            lp.a_matrix_.start_.end());
    model.elements_variables = std::vector<int>(
            lp.a_matrix_.index_.begin(),
            lp.a_matrix_.index_.end());
    model.elements_coefficients = lp.a_matrix_.value_;

    return model;
}

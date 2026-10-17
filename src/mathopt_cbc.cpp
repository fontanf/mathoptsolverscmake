#include "mathoptsolverscmake/mathopt_cbc.hpp"

#include <cstdio>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace mathoptsolverscmake;

namespace
{

/**
 * Message handler which writes to a file that it owns.
 */
class FileMessageHandler: public CoinMessageHandler
{

public:

    FileMessageHandler(const std::string& log_file):
        file_(std::fopen(log_file.c_str(), "w"))
    {
        if (file_ == nullptr) {
            throw std::runtime_error(
                    FUNC_SIGNATURE + ": "
                    "unable to open file \"" + log_file + "\".");
        }
        setFilePointer(file_);
    }

    FileMessageHandler(const FileMessageHandler&) = delete;

    FileMessageHandler& operator=(const FileMessageHandler&) = delete;

    virtual ~FileMessageHandler() { std::fclose(file_); }

private:

    FILE* file_;

};

}

mathoptsolverscmake::CbcModel::CbcModel():
    model(OsiCbcSolverInterface())
{
}

void mathoptsolverscmake::load(
        CbcModel& cbc_model,
        const MathOptModel& model)
{
    if (!model.is_milp()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "model is not an MILP; Cbc only supports MILP models.");
    }

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
    cbc_model.model.solver()->loadProblem(
            matrix,
            model.variables_lower_bounds.data(),
            model.variables_upper_bounds.data(),
            model.objective_coefficients.data(),
            model.constraints_lower_bounds.data(),
            model.constraints_upper_bounds.data());

    // Set objective direction.
    if (model.objective_direction == ObjectiveDirection::Minimize) {
        cbc_model.model.setObjSense(1);
    } else {
        cbc_model.model.setObjSense(-1);
    }

    // Set variable types.
    std::vector<int> integer_variables;
    for (int variable_id = 0;
            variable_id < model.number_of_variables();
            ++variable_id) {
        if (model.variables_types[variable_id] != VariableType::Continuous)
            integer_variables.push_back(variable_id);
    }
    cbc_model.model.solver()->setInteger(integer_variables.data(), integer_variables.size());
}

void mathoptsolverscmake::reduce_printout(
        CbcModel& cbc_model)
{
    cbc_model.model.setLogLevel(0);
    cbc_model.model.messageHandler()->setLogLevel(0);
    cbc_model.model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
}

void mathoptsolverscmake::set_time_limit(
        CbcModel& cbc_model,
        double time_limit)
{
    cbc_model.model.setMaximumSeconds(time_limit);
}

void mathoptsolverscmake::set_log_file(
        CbcModel& cbc_model,
        const std::string& log_file)
{
    std::unique_ptr<CoinMessageHandler> log_handler(
            new FileMessageHandler(log_file));
    log_handler->setLogLevel(1);
    cbc_model.model.passInMessageHandler(log_handler.get());
    // The previous handler, if any, isn't used by the model anymore.
    cbc_model.log_handler = std::move(log_handler);
}

void mathoptsolverscmake::solve(
        CbcModel& cbc_model)
{
    cbc_model.model.branchAndBound();
}

double mathoptsolverscmake::get_solution_value(
        const CbcModel& cbc_model)
{
    return get_solution_value(cbc_model.model);
}

std::vector<double> mathoptsolverscmake::get_solution(
        const CbcModel& cbc_model)
{
    return get_solution(cbc_model.model);
}

double mathoptsolverscmake::get_bound(
        const CbcModel& cbc_model)
{
    return get_bound(cbc_model.model);
}

int mathoptsolverscmake::get_number_of_nodes(
        const CbcModel& cbc_model)
{
    return get_number_of_nodes(cbc_model.model);
}

double mathoptsolverscmake::get_solution_value(
        const ::CbcModel& cbc_model)
{
    if (cbc_model.bestSolution() == nullptr)
        return cbc_model.getObjSense() * std::numeric_limits<double>::infinity();
    return cbc_model.getObjValue();
}

std::vector<double> mathoptsolverscmake::get_solution(
        const ::CbcModel& cbc_model)
{
    if (cbc_model.bestSolution() == nullptr)
        return {};
    return std::vector<double>(
            cbc_model.bestSolution(),
            cbc_model.bestSolution() + cbc_model.getNumCols());
}

double mathoptsolverscmake::get_bound(
        const ::CbcModel& cbc_model)
{
    return cbc_model.getBestPossibleObjValue();
}

int mathoptsolverscmake::get_number_of_nodes(
        const ::CbcModel& cbc_model)
{
    return cbc_model.getNodeCount();
}

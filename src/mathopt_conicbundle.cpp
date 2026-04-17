#include "mathoptsolverscmake/mathopt_conicbundle.hpp"

#include <stdexcept>
#include <vector>

using namespace mathoptsolverscmake;

class ConicBundleFunction: public ConicBundle::FunctionOracle
{

public:

    ConicBundleFunction(
            const MathOptModel& model):
        model_(model) { }

    int evaluate(
            const double* x,
            double /* relprec */,
            double& objective_value,
            std::vector<ConicBundle::Minorant*>& minorants,
            ConicBundle::PrimalExtender*&)
    {
        std::vector<double> x_vec(x, x + this->model_.number_of_variables());
        BlackBoxFunctionOutput sp_output = this->model_.objective_function(x_vec);
        if (this->model_.objective_direction == ObjectiveDirection::Minimize) {
            objective_value = sp_output.objective_value;
            minorants.push_back(new ConicBundle::Minorant(objective_value, sp_output.gradient));
        } else {
            objective_value = -sp_output.objective_value;
            std::vector<double> gradient(this->model_.number_of_variables(), 0);
            for (int variable_id = 0;
                    variable_id < this->model_.number_of_variables();
                    ++variable_id) {
                gradient[variable_id] = -sp_output.gradient[variable_id];
            }
            minorants.push_back(new ConicBundle::Minorant(objective_value, gradient));
        }
        return 0;
    }

private:

    const MathOptModel& model_;

};

void mathoptsolverscmake::solve(
        const MathOptModel& model,
        ConicBundle::CBSolver& solver)
{
    if (!model.is_box_constrained()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "model is not box-constrained; ConicBundle only supports box-constrained models.");
    }

    ConicBundleFunction cb_function(model);
    solver.init_problem(
            model.number_of_variables(),
            &model.variables_lower_bounds,
            &model.variables_upper_bounds);
    solver.add_function(cb_function);
    solver.solve();
}

double mathoptsolverscmake::get_solution_value(
        const MathOptModel& model,
        const ConicBundle::CBSolver& solver)
{
    return (model.objective_direction == ObjectiveDirection::Minimize)?
        solver.get_objval():
        -solver.get_objval();
}

std::vector<double> mathoptsolverscmake::get_solution(
        const MathOptModel& model,
        const ConicBundle::CBSolver& solver)
{
    std::vector<double> solution(model.number_of_variables(), 0.0);
    ConicBundle::DVector variables;
    solver.get_center(variables);

    if (model.objective_direction == ObjectiveDirection::Minimize) {
        for (int variable_id = 0;
                variable_id < model.number_of_variables();
                ++variable_id) {
            solution[variable_id] = variables[variable_id];
        }
    } else {
        for (int variable_id = 0;
                variable_id < model.number_of_variables();
                ++variable_id) {
            solution[variable_id] = -variables[variable_id];
        }
    }
    return solution;
}

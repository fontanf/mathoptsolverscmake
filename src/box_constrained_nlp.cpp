#include "mathoptsolverscmake/box_constrained_nlp.hpp"

#ifdef DLIB_FOUND
#include "dlib/optimization.h"
#endif

using namespace mathoptsolverscmake;

#ifdef KNITRO_FOUND

inline double to_knitro(double value)
{
    if (value == std::numeric_limits<double>::infinity())
        return KN_INFINITY;
    if (value == -std::numeric_limits<double>::infinity())
        return -KN_INFINITY;
    return value;
}

void mathoptsolverscmake::solve(
        const BoxConstrainedNlpModel& model,
        knitrocpp::Context& knitro_context)
{
    knitro_context.add_vars(model.number_of_variables());
    if (!model.variables_initial_values.empty())
        knitro_context.set_var_primal_init_values(model.variables_initial_values);
    std::vector<double> knitro_lower_bounds(model.number_of_variables(), 0.0);
    std::vector<double> knitro_upper_bounds(model.number_of_variables(), 0.0);
    for (int variable_id = 0; variable_id < model.number_of_variables(); ++variable_id) {
        knitro_lower_bounds[variable_id] = to_knitro(model.variables_lower_bounds[variable_id]);
        knitro_upper_bounds[variable_id] = to_knitro(model.variables_upper_bounds[variable_id]);
    }
    knitro_context.set_var_lobnds(knitro_lower_bounds);
    knitro_context.set_var_upbnds(knitro_upper_bounds);

    if (model.objective_direction == ObjectiveDirection::Minimize) {
        knitro_context.set_obj_goal(KN_OBJGOAL_MINIMIZE);
    } else {
        knitro_context.set_obj_goal(KN_OBJGOAL_MAXIMIZE);
    }
    std::vector<double> x(model.number_of_variables(), 0.0);
    BoxConstrainedNlpFunctionOutput f_output;
    CB_context* callback_context = knitro_context.add_eval_callback(
            true,  // evaluate objective?
            {},  // constraints
            [&model, &x, &f_output](
                const knitrocpp::Context&,
                CB_context*,
                KN_eval_request_ptr const eval_request,
                KN_eval_result_ptr const eval_result)
            {
                for (int variable_id = 0; variable_id < model.number_of_variables(); ++variable_id)
                    x[variable_id] = eval_request->x[variable_id];
                f_output = model.objective_function(x);
                *eval_result->obj = f_output.objective_value;
                return 0;
            });
    knitro_context.set_cb_grad(
            callback_context,
            nullptr,
            nullptr,
            nullptr,
            [&model, &f_output](
                const knitrocpp::Context&,
                CB_context*,
                KN_eval_request_ptr const eval_request,
                KN_eval_result_ptr const eval_result)
            {
                for (int variable_id = 0; variable_id < model.number_of_variables(); ++variable_id)
                    eval_result->objGrad[variable_id] = f_output.gradient[variable_id];
                return 0;
            });

    // Solve
    int knitro_return_status = knitro_context.solve();
}

double mathoptsolverscmake::get_solution_value(
        const knitrocpp::Context& knitro_context)
{
    return knitro_context.get_obj_value();
}

std::vector<double> mathoptsolverscmake::get_solution(
        const knitrocpp::Context& knitro_context)
{
    return knitro_context.get_var_primal_values();
}

#endif

#ifdef DLIB_FOUND

typedef dlib::matrix<double,0,1> column_vector;

class BoxConstrainedNlpDlibFunction
{

public:

    BoxConstrainedNlpDlibFunction(
            const BoxConstrainedNlpModel& model):
        model_(model),
        grad_(model.number_of_variables()) { }

    double f(const column_vector& x)
    {
        std::vector<double> variables(model_.number_of_variables());
        for (int variable_id = 0;
                variable_id < model_.number_of_variables();
                ++variable_id) {
            variables[variable_id] = x(variable_id);
        }
        BoxConstrainedNlpFunctionOutput output = model_.objective_function(variables);

        for (int variable_id = 0;
                variable_id < model_.number_of_variables();
                ++variable_id) {
            grad_(variable_id) = output.gradient[variable_id];
        }
        return output.objective_value;
    }

    const column_vector der(const column_vector& x) const { (void)x; return grad_; }

private:

    const BoxConstrainedNlpModel& model_;

    column_vector grad_;

};

BoxConstrainedNlpDlibOutput mathoptsolverscmake::solve_dlib(
        const BoxConstrainedNlpModel& model)
{
    BoxConstrainedNlpDlibOutput output;

    column_vector dlib_variables(model.number_of_variables());
    if (!model.variables_initial_values.empty()) {
        for (int variable_id = 0;
                variable_id < model.number_of_variables();
                ++variable_id) {
            dlib_variables(variable_id) = model.variables_initial_values[variable_id];
        }
    }

    column_vector dlib_variables_lower_bounds(model.number_of_variables());
    column_vector dlib_variables_upper_bounds(model.number_of_variables());
    for (int variable_id = 0;
            variable_id < model.number_of_variables();
            ++variable_id) {
        dlib_variables_lower_bounds(variable_id) = model.variables_lower_bounds[variable_id];
        dlib_variables_upper_bounds(variable_id) = model.variables_upper_bounds[variable_id];
    }

    // Solve
    BoxConstrainedNlpDlibFunction dlib_function(model);
    auto f = [&dlib_function](const column_vector& x) { return dlib_function.f(x); };
    auto def = [&dlib_function](const column_vector& x) { return dlib_function.der(x); };
    //auto stop_strategy = dlib::objective_delta_stop_strategy(0.0001);
    auto stop_strategy = dlib::objective_delta_stop_strategy(0.0001).be_verbose();
    output.objective_value = (model.objective_direction == ObjectiveDirection::Minimize)?
        dlib::find_min_box_constrained(
            dlib::lbfgs_search_strategy(256),
            stop_strategy,
            f,
            def,
            dlib_variables,
            dlib_variables_lower_bounds,
            dlib_variables_upper_bounds):
        dlib::find_max_box_constrained(
            dlib::lbfgs_search_strategy(256),
            stop_strategy,
            f,
            def,
            dlib_variables,
            dlib_variables_lower_bounds,
            dlib_variables_upper_bounds);

    // Compute output parameters
    output.solution = std::vector<double>(model.number_of_variables());
    for (int variable_id = 0;
            variable_id < model.number_of_variables();
            ++variable_id) {
        output.solution[variable_id] = dlib_variables(variable_id);
    }
    return output;
}

#endif

#ifdef CONICBUNDLE_FOUND

class BoxConstrainedNlpConicBundleFunction: public ConicBundle::FunctionOracle
{

public:

    BoxConstrainedNlpConicBundleFunction(
            const BoxConstrainedNlpModel& model):
        model_(model) { }

    int evaluate(
            const double* x,
            double /* relprec */,
            double& objective_value,
            std::vector<ConicBundle::Minorant*>& minorants,
            ConicBundle::PrimalExtender*&)
    {
        std::vector<double> x_vec(x, x + this->model_.number_of_variables());
        BoxConstrainedNlpFunctionOutput sp_output = this->model_.objective_function(x_vec);
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

    const BoxConstrainedNlpModel& model_;

};

void mathoptsolverscmake::solve(
        const BoxConstrainedNlpModel& model,
        ConicBundle::CBSolver& solver)
{
    BoxConstrainedNlpConicBundleFunction cb_function(model);
    solver.init_problem(
            model.number_of_variables(),
            &model.variables_lower_bounds,
            &model.variables_upper_bounds);
    solver.add_function(cb_function);
    solver.solve();
}

double mathoptsolverscmake::get_solution_value(
        const BoxConstrainedNlpModel& model,
        const ConicBundle::CBSolver& solver)
{
    return (model.objective_direction == ObjectiveDirection::Minimize)?
        solver.get_objval():
        -solver.get_objval();
}

std::vector<double> mathoptsolverscmake::get_solution(
        const BoxConstrainedNlpModel& model,
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

#endif

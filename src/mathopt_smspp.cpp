#include "mathoptsolverscmake/mathopt_smspp.hpp"

#include "AbstractBlock.h"
#include "BundleSolver.h"
#include "FRealObjective.h"
#include "OneVarConstraint.h"

#include <algorithm>
#include <limits>
#include <list>
#include <memory>
#include <stdexcept>
#include <vector>

using namespace mathoptsolverscmake;
using namespace SMSpp_di_unipi_it;

namespace
{

/**
 * A C05Function wrapping a MathOptModel's black-box objective (a callback
 * returning a value and a gradient) as an SMS++ nonlinear-oracle Function.
 *
 * The active Variable are exactly model.objective_function's input vector,
 * in order; they never change (remove_variable() is unsupported), and the
 * Function only ever exposes the single "diagonal" linearization (the
 * gradient) computed by the last compute() call, i.e. it does not use the
 * "global pool" machinery (any name other than the default, meaning "the
 * last computed one", is rejected).
 */
class MathOptModelC05Function: public C05Function
{

public:

    MathOptModelC05Function(
            const MathOptModel& model,
            std::vector<ColVariable*>&& variables):
        model_(model),
        v_vars(std::move(variables)),
        gradient_(v_vars.size(), 0.0) { }

    /*
     * ThinVarDepInterface: the "active" Variable are exactly v_vars, fixed
     * once and for all at construction time.
     */

    class v_iterator: public ThinVarDepInterface::v_iterator
    {

    public:

        v_iterator(std::vector<ColVariable*>::iterator it): it_(it) { }

        v_iterator* clone() override { return new v_iterator(it_); }

        void operator++() override { ++it_; }
        reference operator*() const override { return *(*it_); }
        pointer operator->() const override { return *it_; }

        bool operator==(const ThinVarDepInterface::v_iterator& rhs) const override
        {
            auto other = dynamic_cast<const v_iterator*>(&rhs);
            return other? it_ == other->it_: false;
        }

        bool operator!=(const ThinVarDepInterface::v_iterator& rhs) const override
        {
            auto other = dynamic_cast<const v_iterator*>(&rhs);
            return other? it_ != other->it_: true;
        }

    private:

        std::vector<ColVariable*>::iterator it_;

    };

    class v_const_iterator: public ThinVarDepInterface::v_const_iterator
    {

    public:

        v_const_iterator(std::vector<ColVariable*>::const_iterator it): it_(it) { }

        v_const_iterator* clone() override { return new v_const_iterator(it_); }

        void operator++() override { ++it_; }
        reference operator*() const override { return *(*it_); }
        pointer operator->() const override { return *it_; }

        bool operator==(const ThinVarDepInterface::v_const_iterator& rhs) const override
        {
            auto other = dynamic_cast<const v_const_iterator*>(&rhs);
            return other? it_ == other->it_: false;
        }

        bool operator!=(const ThinVarDepInterface::v_const_iterator& rhs) const override
        {
            auto other = dynamic_cast<const v_const_iterator*>(&rhs);
            return other? it_ != other->it_: true;
        }

    private:

        std::vector<ColVariable*>::const_iterator it_;

    };

    Index get_num_active_var() const override final
    {
        return (Index)v_vars.size();
    }

    Index is_active(const Variable* var) const override final
    {
        auto it = std::find(v_vars.begin(), v_vars.end(), var);
        return it != v_vars.end()?
            (Index)std::distance(v_vars.begin(), it):
            Inf<Index>();
    }

    Variable* get_active_var(const Index i) const override final
    {
        return v_vars[i];
    }

    v_iterator* v_begin() override final
    {
        return new v_iterator(v_vars.begin());
    }

    v_const_iterator* v_begin() const override final
    {
        return new v_const_iterator(v_vars.begin());
    }

    v_iterator* v_end() override final
    {
        return new v_iterator(v_vars.end());
    }

    v_const_iterator* v_end() const override final
    {
        return new v_const_iterator(v_vars.end());
    }

    void remove_variable(Index, c_ModParam = eModBlck) override final
    {
        throw std::logic_error(
                FUNC_SIGNATURE + ": removing variables is not supported.");
    }

    /*
     * Function / C05Function.
     */

    int compute(bool /* changedvars */ = true) override
    {
        std::vector<double> x(v_vars.size());
        for (std::size_t i = 0; i < v_vars.size(); ++i)
            x[i] = v_vars[i]->get_value();

        BlackBoxFunctionOutput output = this->model_.objective_function(x);
        if (this->model_.objective_direction == ObjectiveDirection::Minimize) {
            value_ = output.objective_value;
            gradient_ = output.gradient;
        } else {
            value_ = -output.objective_value;
            for (std::size_t i = 0; i < v_vars.size(); ++i)
                gradient_[i] = -output.gradient[i];
        }
        return kOK;
    }

    FunctionValue get_value() override { return value_; }

    // BundleSolver::set_Block() hard-requires the objective's Function to
    // report is_convex() or is_concave() (it throws otherwise), and further
    // requires the Objective's sense to be Minimize for a convex function or
    // Maximize for a concave one. Since compute() above always normalizes to
    // a Minimize (negating value and gradient when the model itself
    // maximizes, exactly like this codebase's ConicBundle wrapper does), this
    // Function must always report itself convex here to satisfy that check.
    // This is not a claim this wrapper can verify: BundleSolver is only
    // correct if model.objective_function truly is convex when minimizing
    // (or concave when maximizing) -- the primary intended use, per this
    // library's Lagrangian-relaxation-oriented box-constrained NLP modeler,
    // where that always holds by construction.

    bool is_convex() override { return true; }

    bool is_continuously_differentiable() const override { return true; }

    // BundleSolver drives its bundle method by keeping a "global pool" of
    // past linearizations (set via intGPMaxSz, a param on every C05Function);
    // the base class's set_par(int) throws on any attempt to change it, so it
    // has to be overridden here with a real (if simple, dense-vector-backed)
    // pool, sized to whatever BundleSolver asks for.

    void set_par(idx_type par, int value) override
    {
        if (par == intGPMaxSz) {
            gp_size_ = (value > 0)? (Index)value: 0;
            pool_gradient_.resize(gp_size_);
            pool_alpha_.resize(gp_size_, 0.0);
            pool_valid_.assign(gp_size_, false);
            return;
        }
        if (par == intLPMaxSz) {
            // Stored (so get_int_par() below reports it back accurately) but
            // otherwise unused: compute_new_linearization() always returns
            // false (the gradient is the only linearization this Function
            // ever produces), so the local pool never holds more than one
            // entry regardless of the capacity BundleSolver asks for.
            lp_size_ = value;
            return;
        }
        C05Function::set_par(par, value);
    }

    // The base class's get_int_par() (called by, e.g., C05Function's default
    // delete_linearizations() to learn the pool size) just returns the
    // compile-time default unless overridden -- it does not know that
    // set_par() above actually applies the value, so without this override
    // it would report a stale 0-sized pool no matter what set_par() had set,
    // and any pool_gradient_/pool_alpha_ index would look "out of range".

    int get_int_par(idx_type par) const override
    {
        if (par == intGPMaxSz)
            return (int)gp_size_;
        if (par == intLPMaxSz)
            return lp_size_;
        return C05Function::get_int_par(par);
    }

    using C05Function::set_par;

    void store_linearization(Index name, ModParam /* issueMod */ = eModBlck) override
    {
        check_pool_name(name);
        pool_gradient_[name] = gradient_;
        pool_alpha_[name] = linearization_constant();
        pool_valid_[name] = true;
    }

    bool is_linearization_there(Index name) const override
    {
        return (name < gp_size_) && pool_valid_[name];
    }

    void delete_linearization(Index name, ModParam /* issueMod */ = eModBlck) override
    {
        if (name < gp_size_)
            pool_valid_[name] = false;
    }

    void get_linearization_coefficients(
            FunctionValue* g,
            Range range = INFRange,
            Index name = Inf<Index>()) override
    {
        const std::vector<double>& grad = linearization_gradient(name);
        range.second = std::min(range.second, get_num_active_var());
        for (Index i = range.first; i < range.second; ++i)
            *(g++) = grad[i];
    }

    void get_linearization_coefficients(
            FunctionValue* g,
            c_Subset& subset,
            bool /* ordered */ = false,
            Index name = Inf<Index>()) override
    {
        const std::vector<double>& grad = linearization_gradient(name);
        for (Index i: subset)
            *(g++) = grad[i];
    }

    FunctionValue get_linearization_constant(Index name = Inf<Index>()) override
    {
        if (name == Inf<Index>())
            return linearization_constant();
        check_pool_name(name);
        if (!pool_valid_[name]) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": linearization not in the global pool.");
        }
        return pool_alpha_[name];
    }

private:

    // The constant term alpha of the tangent-plane linearization L(x) =
    // alpha + gradient . x through the point of the last compute() call:
    // since the gradient is exact, alpha = f(x) - gradient . x.
    FunctionValue linearization_constant() const
    {
        FunctionValue alpha = value_;
        for (std::size_t i = 0; i < v_vars.size(); ++i)
            alpha -= v_vars[i]->get_value() * gradient_[i];
        return alpha;
    }

    const std::vector<double>& linearization_gradient(Index name) const
    {
        if (name == Inf<Index>())
            return gradient_;
        check_pool_name(name);
        if (!pool_valid_[name]) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": linearization not in the global pool.");
        }
        return pool_gradient_[name];
    }

    void check_pool_name(Index name) const
    {
        if (name >= gp_size_) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": linearization name out of range.");
        }
    }

    Index gp_size_ = 0;
    int lp_size_ = 1;
    std::vector<std::vector<double>> pool_gradient_;
    std::vector<double> pool_alpha_;
    std::vector<bool> pool_valid_;

    const MathOptModel& model_;

    std::vector<ColVariable*> v_vars;

    double value_ = 0;
    std::vector<double> gradient_;

};

}

SmsppOutput mathoptsolverscmake::solve_smspp(
        const MathOptModel& model)
{
    if (!model.is_box_constrained()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "model is not box-constrained; SMS++ only supports box-constrained models.");
    }
    if (!model.objective_function) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "model.objective_function must be set; SMS++'s BundleSolver requires a "
                "black-box objective providing a value and a gradient.");
    }
    // BundleSolver::FakeFiOracle::GetUC() only accepts a variable lower
    // bound that is exactly 0 or -infinity (it follows the Lagrangian-
    // multiplier convention lambda >= 0 or lambda free), throwing
    // std::logic_error("finite lhs different from zero not allowed")
    // otherwise; checked eagerly here for a clearer error message.
    for (int variable_id = 0; variable_id < model.number_of_variables(); ++variable_id) {
        double lb = model.variables_lower_bounds[variable_id];
        if (lb != -std::numeric_limits<double>::infinity() && lb != 0.0) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "variable " + std::to_string(variable_id) + " has lower bound " +
                    std::to_string(lb) + "; SMS++'s BundleSolver only accepts a lower "
                    "bound of 0 or -infinity (the Lagrangian-multiplier convention).");
        }
    }

    // AbstractBlock deletes, in its own destructor, every Variable/Constraint
    // container and every Objective it was given, so they must all be
    // heap-allocated (mirroring the SMS++ tests, e.g. tests/PolyhedralFunction).
    std::unique_ptr<AbstractBlock> block(new AbstractBlock());

    auto variables = new std::vector<ColVariable>(model.number_of_variables());
    block->add_static_variable(*variables, "x");

    std::vector<ColVariable*> variable_pointers(model.number_of_variables());
    for (int variable_id = 0; variable_id < model.number_of_variables(); ++variable_id) {
        variable_pointers[variable_id] = &(*variables)[variable_id];
        if (!model.variables_initial_values.empty())
            (*variables)[variable_id].set_value(model.variables_initial_values[variable_id]);
    }

    auto bounds = new std::list<BoxConstraint>();
    for (int variable_id = 0; variable_id < model.number_of_variables(); ++variable_id) {
        bounds->emplace_back(
                block.get(),
                &(*variables)[variable_id],
                model.variables_lower_bounds[variable_id],
                model.variables_upper_bounds[variable_id]);
    }
    block->add_dynamic_constraint(*bounds, "xbnd");

    auto objective = new FRealObjective();
    objective->set_function(new MathOptModelC05Function(model, std::move(variable_pointers)));
    objective->set_sense(Objective::eMin, eNoMod);
    block->set_objective(objective);

    // AbstractBlock::~AbstractBlock() deletes the Variable containers before
    // the Objective (see its implementation), but the Objective's Function
    // here holds pointers into those same Variable: were it left to that
    // destructor, ~FRealObjective() would try to unregister itself from
    // already-freed ColVariable while unwinding, a use-after-free. So the
    // Objective (and, through it, the Function) is always explicitly
    // detached and deleted below -- on every exit path, including on
    // exceptions -- strictly before block's destructor gets a chance to run.
    // Block::set_objective() unconditionally dereferences the new Objective
    // (newOF->set_Block(this), no null check), so the old one can't simply be
    // swapped for nullptr; an empty placeholder (no Function, so nothing of
    // ours left for ~AbstractBlock() to touch) stands in instead.
    auto detach_objective = [&block, &objective]()
    {
        delete objective;
        block->set_objective(new FRealObjective(), eNoMod);
    };

    auto solver = new BundleSolver();
    try {
        block->register_Solver(solver);
    } catch (...) {
        delete solver;  // set_Block() threw, so it was never actually registered
        detach_objective();
        throw;
    }

    int return_code;
    try {
        return_code = solver->compute(false);
    } catch (...) {
        block->unregister_Solver(solver, true);
        detach_objective();
        throw;
    }
    bool has_solution =
        ((return_code >= Solver::kOK) && (return_code < Solver::kError))
        || (return_code == Solver::kLowPrecision);
    if (!has_solution) {
        block->unregister_Solver(solver, true);
        detach_objective();
        throw std::runtime_error(
                FUNC_SIGNATURE + ": "
                "SMS++'s BundleSolver did not find a solution (return code "
                + std::to_string(return_code) + ").");
    }

    solver->get_var_solution();

    SmsppOutput output;
    output.solution.resize(model.number_of_variables());
    for (int variable_id = 0; variable_id < model.number_of_variables(); ++variable_id)
        output.solution[variable_id] = (*variables)[variable_id].get_value();

    double smspp_value = solver->get_ub();
    output.objective_value = (model.objective_direction == ObjectiveDirection::Minimize)?
        smspp_value:
        -smspp_value;

    block->unregister_Solver(solver, true);
    detach_objective();
    return output;
}

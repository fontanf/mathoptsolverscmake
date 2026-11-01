#include "mathoptsolverscmake/mathopt_ipopt.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace mathoptsolverscmake;

namespace
{

/** Value used for infinite bounds (Ipopt's default "nlp_upper_bound_inf" is 1e19). */
const double ipopt_infinity = 2e19;

double to_ipopt(double value)
{
    if (value >= ipopt_infinity)
        return ipopt_infinity;
    if (value <= -ipopt_infinity)
        return -ipopt_infinity;
    return value;
}

/**
 * A nonlinear expression tree (stored as parallel element arrays, see
 * mathopt.hpp), with its exact first and second derivatives.
 *
 * Nodes are stored in DFS pre-order, so children always have larger indices
 * than their parent: scanning nodes backward evaluates children before their
 * parent, scanning them forward propagates derivatives from parents to
 * children.
 *
 * - The gradient is computed by reverse mode differentiation: the adjoint
 *   of a node is the derivative of the root with respect to it.
 * - The Hessian is computed by forward-over-reverse differentiation: for
 *   each variable j, a forward pass computes the derivative of each node
 *   with respect to x_j ("tangent"), then a reverse pass computes the
 *   derivative of the adjoints with respect to x_j, which, at the variable
 *   nodes, give column j of the Hessian.
 */
class ExpressionTree
{

public:

    ExpressionTree(
            const std::vector<char>& operators,
            const std::vector<double>& values,
            const std::vector<int>& variables,
            const std::vector<int>& number_of_children,
            int start,
            int end):
        operators_(operators.begin() + start, operators.begin() + end),
        values_(values.begin() + start, values.begin() + end),
        node_variables_(end - start, -1),
        children_starts_(end - start, 0),
        number_of_children_(number_of_children.begin() + start, number_of_children.begin() + end),
        has_variables_(end - start, false),
        second_partials_starts_(end - start, -1),
        node_values_(end - start, 0.0),
        adjoints_(end - start, 0.0),
        tangents_(end - start, 0.0),
        adjoint_tangents_(end - start, 0.0)
    {
        int number_of_nodes = end - start;

        // Local variables.
        std::map<int, int> local_ids;
        for (int node_id = 0; node_id < number_of_nodes; ++node_id) {
            if (operators_[node_id] != 'v')
                continue;
            int variable_id = variables[start + node_id];
            auto it = local_ids.find(variable_id);
            if (it == local_ids.end()) {
                it = local_ids.insert({variable_id, (int)variables_.size()}).first;
                variables_.push_back(variable_id);
            }
            node_variables_[node_id] = it->second;
        }

        // Children, with the same stack-based scan as the evaluation: the
        // first child popped is the first (left) operand.
        std::vector<int> stk;
        for (int node_id = number_of_nodes - 1; node_id >= 0; --node_id) {
            children_starts_[node_id] = (int)children_.size();
            for (int child_pos = 0; child_pos < number_of_children_[node_id]; ++child_pos) {
                int child_id = stk.back();
                stk.pop_back();
                children_.push_back(child_id);
                if (has_variables_[child_id])
                    has_variables_[node_id] = true;
            }
            if (operators_[node_id] == 'v')
                has_variables_[node_id] = true;
            stk.push_back(node_id);
        }

        // Partial derivatives storage.
        int number_of_partials = 0;
        int number_of_second_partials = 0;
        for (int node_id = 0; node_id < number_of_nodes; ++node_id) {
            number_of_partials += number_of_children_[node_id];
            switch (operators_[node_id]) {
            case '+': case '-': case 'n': case 'k': case 'v':
                break;
            default:
                second_partials_starts_[node_id] = number_of_second_partials;
                number_of_second_partials += number_of_children_[node_id] * number_of_children_[node_id];
            }
        }
        partials_.resize(number_of_partials, 0.0);
        second_partials_.resize(number_of_second_partials, 0.0);
        gradient_.resize(variables_.size(), 0.0);
    }

    /** Get the (global) ids of the variables of the tree, by local id. */
    const std::vector<int>& variables() const { return variables_; }

    /**
     * Evaluate the tree, its partial derivatives and its gradient at x
     * (indexed by global variable ids).
     */
    void evaluate(const double* x)
    {
        int number_of_nodes = (int)operators_.size();

        // Values.
        for (int node_id = number_of_nodes - 1; node_id >= 0; --node_id) {
            const int* children = children_.data() + children_starts_[node_id];
            double& value = node_values_[node_id];
            switch (operators_[node_id]) {
            case 'k': value = values_[node_id]; break;
            case 'v': value = x[variables_[node_variables_[node_id]]]; break;
            case 'n': value = -node_values_[children[0]]; break;
            case 'e': value = std::exp(node_values_[children[0]]); break;
            case 'l': value = std::log(node_values_[children[0]]); break;
            case 'q': value = std::sqrt(node_values_[children[0]]); break;
            case 's': value = std::sin(node_values_[children[0]]); break;
            case 'c': value = std::cos(node_values_[children[0]]); break;
            case 't': value = std::tan(node_values_[children[0]]); break;
            case '+': {
                value = 0.0;
                for (int child_pos = 0; child_pos < number_of_children_[node_id]; ++child_pos)
                    value += node_values_[children[child_pos]];
                break;
            } case '*': {
                value = 1.0;
                for (int child_pos = 0; child_pos < number_of_children_[node_id]; ++child_pos)
                    value *= node_values_[children[child_pos]];
                break;
            }
            case '-': value = node_values_[children[0]] - node_values_[children[1]]; break;
            case '/': value = node_values_[children[0]] / node_values_[children[1]]; break;
            case 'p': value = std::pow(node_values_[children[0]], node_values_[children[1]]); break;
            default:
                throw std::invalid_argument(
                        FUNC_SIGNATURE + ": "
                        "unknown nonlinear operator '"
                        + std::string(1, operators_[node_id]) + "'.");
            }
        }

        // Partial derivatives.
        for (int node_id = 0; node_id < number_of_nodes; ++node_id)
            compute_partials(node_id);

        // Adjoints and gradient.
        std::fill(adjoints_.begin(), adjoints_.end(), 0.0);
        std::fill(gradient_.begin(), gradient_.end(), 0.0);
        adjoints_[0] = 1.0;
        for (int node_id = 0; node_id < number_of_nodes; ++node_id) {
            if (!has_variables_[node_id])
                continue;
            if (operators_[node_id] == 'v') {
                gradient_[node_variables_[node_id]] += adjoints_[node_id];
                continue;
            }
            const int* children = children_.data() + children_starts_[node_id];
            const double* partials = partials_.data() + children_starts_[node_id];
            for (int child_pos = 0; child_pos < number_of_children_[node_id]; ++child_pos)
                adjoints_[children[child_pos]] += adjoints_[node_id] * partials[child_pos];
        }
    }

    /** Get the value of the tree (requires evaluate()). */
    double value() const { return node_values_[0]; }

    /** Get the gradient of the tree, by local variable id (requires evaluate()). */
    const std::vector<double>& gradient() const { return gradient_; }

    /**
     * Add "weight" times the lower triangle of the Hessian of the tree to
     * "hessian", where the entry of local variables a >= b is at index
     * a * (a + 1) / 2 + b (requires evaluate()).
     */
    void add_hessian(
            double weight,
            std::vector<double>& hessian)
    {
        int number_of_nodes = (int)operators_.size();
        for (int local_variable_id = 0;
                local_variable_id < (int)variables_.size();
                ++local_variable_id) {

            // Tangents.
            for (int node_id = number_of_nodes - 1; node_id >= 0; --node_id) {
                double& tangent = tangents_[node_id];
                tangent = 0.0;
                if (!has_variables_[node_id])
                    continue;
                if (operators_[node_id] == 'v') {
                    tangent = (node_variables_[node_id] == local_variable_id)? 1.0: 0.0;
                    continue;
                }
                const int* children = children_.data() + children_starts_[node_id];
                const double* partials = partials_.data() + children_starts_[node_id];
                for (int child_pos = 0; child_pos < number_of_children_[node_id]; ++child_pos)
                    tangent += partials[child_pos] * tangents_[children[child_pos]];
            }

            // Tangents of the adjoints.
            std::fill(adjoint_tangents_.begin(), adjoint_tangents_.end(), 0.0);
            for (int node_id = 0; node_id < number_of_nodes; ++node_id) {
                if (!has_variables_[node_id])
                    continue;
                if (operators_[node_id] == 'v') {
                    int variable_a = node_variables_[node_id];
                    if (variable_a >= local_variable_id) {
                        hessian[variable_a * (variable_a + 1) / 2 + local_variable_id]
                            += weight * adjoint_tangents_[node_id];
                    }
                    continue;
                }
                const int* children = children_.data() + children_starts_[node_id];
                const double* partials = partials_.data() + children_starts_[node_id];
                int number_of_children = number_of_children_[node_id];
                for (int child_pos = 0; child_pos < number_of_children; ++child_pos) {
                    double adjoint_tangent = adjoint_tangents_[node_id] * partials[child_pos];
                    if (second_partials_starts_[node_id] != -1) {
                        const double* second_partials = second_partials_.data()
                            + second_partials_starts_[node_id]
                            + child_pos * number_of_children;
                        for (int child_pos_2 = 0; child_pos_2 < number_of_children; ++child_pos_2) {
                            adjoint_tangent += adjoints_[node_id]
                                * second_partials[child_pos_2]
                                * tangents_[children[child_pos_2]];
                        }
                    }
                    adjoint_tangents_[children[child_pos]] += adjoint_tangent;
                }
            }
        }
    }

private:

    /**
     * Compute the first and second partial derivatives of a node with
     * respect to its children. Derivatives with respect to children which
     * don't depend on any variable are not needed and set to 0 (this
     * avoids, for example, computing log(a) for a ^ b with b constant).
     */
    void compute_partials(int node_id)
    {
        int number_of_children = number_of_children_[node_id];
        if (number_of_children == 0)
            return;
        const int* children = children_.data() + children_starts_[node_id];
        double* p = partials_.data() + children_starts_[node_id];
        double* pp = (second_partials_starts_[node_id] == -1)?
            nullptr:
            second_partials_.data() + second_partials_starts_[node_id];
        if (pp != nullptr)
            std::fill(pp, pp + number_of_children * number_of_children, 0.0);
        double f = node_values_[node_id];
        double a = node_values_[children[0]];
        switch (operators_[node_id]) {
        case 'n': p[0] = -1.0; break;
        case 'e': p[0] = f; pp[0] = f; break;
        case 'l': p[0] = 1.0 / a; pp[0] = -1.0 / (a * a); break;
        case 'q': p[0] = 0.5 / f; pp[0] = -0.25 / (a * f); break;
        case 's': p[0] = std::cos(a); pp[0] = -f; break;
        case 'c': p[0] = -std::sin(a); pp[0] = -f; break;
        case 't': p[0] = 1.0 + f * f; pp[0] = 2.0 * f * (1.0 + f * f); break;
        case '+': {
            for (int child_pos = 0; child_pos < number_of_children; ++child_pos)
                p[child_pos] = 1.0;
            break;
        } case '-': {
            p[0] = 1.0;
            p[1] = -1.0;
            break;
        } case '*': {
            // Products of all factors but one or two, without dividing (some
            // factors may be 0).
            for (int child_pos = 0; child_pos < number_of_children; ++child_pos) {
                double product = 1.0;
                for (int factor_pos = 0; factor_pos < number_of_children; ++factor_pos)
                    if (factor_pos != child_pos)
                        product *= node_values_[children[factor_pos]];
                p[child_pos] = product;
                for (int child_pos_2 = 0; child_pos_2 < number_of_children; ++child_pos_2) {
                    if (child_pos_2 == child_pos)
                        continue;
                    double product_2 = 1.0;
                    for (int factor_pos = 0; factor_pos < number_of_children; ++factor_pos)
                        if (factor_pos != child_pos && factor_pos != child_pos_2)
                            product_2 *= node_values_[children[factor_pos]];
                    pp[child_pos * number_of_children + child_pos_2] = product_2;
                }
            }
            break;
        } case '/': {
            double b = node_values_[children[1]];
            p[0] = 1.0 / b;
            p[1] = -a / (b * b);
            pp[0 * 2 + 1] = -1.0 / (b * b);
            pp[1 * 2 + 0] = -1.0 / (b * b);
            pp[1 * 2 + 1] = 2.0 * a / (b * b * b);
            break;
        } case 'p': {
            double b = node_values_[children[1]];
            p[0] = (b == 0.0)? 0.0: b * std::pow(a, b - 1.0);
            pp[0 * 2 + 0] = (b == 0.0 || b == 1.0)? 0.0: b * (b - 1.0) * std::pow(a, b - 2.0);
            if (has_variables_[children[1]]) {
                double log_a = std::log(a);
                p[1] = f * log_a;
                pp[0 * 2 + 1] = std::pow(a, b - 1.0) * (1.0 + b * log_a);
                pp[1 * 2 + 0] = pp[0 * 2 + 1];
                pp[1 * 2 + 1] = f * log_a * log_a;
            } else {
                p[1] = 0.0;
            }
            break;
        }
        }

        for (int child_pos = 0; child_pos < number_of_children; ++child_pos) {
            if (has_variables_[children[child_pos]])
                continue;
            p[child_pos] = 0.0;
            if (pp == nullptr)
                continue;
            for (int child_pos_2 = 0; child_pos_2 < number_of_children; ++child_pos_2) {
                pp[child_pos * number_of_children + child_pos_2] = 0.0;
                pp[child_pos_2 * number_of_children + child_pos] = 0.0;
            }
        }
    }

    /*
     * Structure.
     */

    std::vector<char> operators_;
    std::vector<double> values_;
    /** Local variable id of each 'v' node, -1 for other nodes. */
    std::vector<int> node_variables_;
    std::vector<int> children_starts_;
    std::vector<int> number_of_children_;
    std::vector<int> children_;
    std::vector<bool> has_variables_;
    std::vector<int> second_partials_starts_;
    std::vector<int> variables_;

    /*
     * Evaluation buffers.
     */

    std::vector<double> node_values_;
    /** Partial derivatives of each node with respect to its children. */
    std::vector<double> partials_;
    /** Second partial derivatives, row-major children x children matrices. */
    std::vector<double> second_partials_;
    std::vector<double> adjoints_;
    std::vector<double> gradient_;
    std::vector<double> tangents_;
    std::vector<double> adjoint_tangents_;

};

class MathOptTnlp: public Ipopt::TNLP
{

public:

    MathOptTnlp(
            const MathOptModel& model,
            IpoptModel& ipopt_model);

    bool exact_hessian() const { return exact_hessian_; }

    bool get_nlp_info(
            Ipopt::Index& n,
            Ipopt::Index& m,
            Ipopt::Index& nnz_jac_g,
            Ipopt::Index& nnz_h_lag,
            IndexStyleEnum& index_style) override;

    bool get_bounds_info(
            Ipopt::Index n,
            Ipopt::Number* x_l,
            Ipopt::Number* x_u,
            Ipopt::Index m,
            Ipopt::Number* g_l,
            Ipopt::Number* g_u) override;

    bool get_starting_point(
            Ipopt::Index n,
            bool init_x,
            Ipopt::Number* x,
            bool init_z,
            Ipopt::Number* z_L,
            Ipopt::Number* z_U,
            Ipopt::Index m,
            bool init_lambda,
            Ipopt::Number* lambda) override;

    bool eval_f(
            Ipopt::Index n,
            const Ipopt::Number* x,
            bool new_x,
            Ipopt::Number& obj_value) override;

    bool eval_grad_f(
            Ipopt::Index n,
            const Ipopt::Number* x,
            bool new_x,
            Ipopt::Number* grad_f) override;

    bool eval_g(
            Ipopt::Index n,
            const Ipopt::Number* x,
            bool new_x,
            Ipopt::Index m,
            Ipopt::Number* g) override;

    bool eval_jac_g(
            Ipopt::Index n,
            const Ipopt::Number* x,
            bool new_x,
            Ipopt::Index m,
            Ipopt::Index nele_jac,
            Ipopt::Index* iRow,
            Ipopt::Index* jCol,
            Ipopt::Number* values) override;

    bool eval_h(
            Ipopt::Index n,
            const Ipopt::Number* x,
            bool new_x,
            Ipopt::Number obj_factor,
            Ipopt::Index m,
            const Ipopt::Number* lambda,
            bool new_lambda,
            Ipopt::Index nele_hess,
            Ipopt::Index* iRow,
            Ipopt::Index* jCol,
            Ipopt::Number* values) override;

    void finalize_solution(
            Ipopt::SolverReturn status,
            Ipopt::Index n,
            const Ipopt::Number* x,
            const Ipopt::Number* z_L,
            const Ipopt::Number* z_U,
            Ipopt::Index m,
            const Ipopt::Number* g,
            const Ipopt::Number* lambda,
            Ipopt::Number obj_value,
            const Ipopt::IpoptData* ip_data,
            Ipopt::IpoptCalculatedQuantities* ip_cq) override;

private:

    /**
     * Update the current point; trees and black-box functions are then
     * evaluated lazily, at most once per point.
     */
    void update_x(const Ipopt::Number* x, bool new_x);

    void evaluate_trees();

    const BlackBoxFunctionOutput& objective_black_box_output();

    const BlackBoxFunctionOutput& constraint_black_box_output(int constraint_id);

    /** Get the position of the Hessian entry (i, j), adding it if needed. */
    int hessian_position(int variable_id_1, int variable_id_2);

    /** Copy of the model (only black-box functions are used after construction). */
    MathOptModel model_;

    IpoptModel& ipopt_model_;

    /** 1 for a minimization problem, -1 for a maximization problem. */
    double sign_;

    bool exact_hessian_;

    std::unique_ptr<ExpressionTree> objective_tree_;

    /** Nonlinear expression tree of each constraint, null if none. */
    std::vector<std::unique_ptr<ExpressionTree>> constraints_trees_;

    /*
     * Jacobian structure.
     */

    std::vector<Ipopt::Index> jacobian_rows_;
    std::vector<Ipopt::Index> jacobian_columns_;
    std::vector<int> linear_jacobian_positions_;
    std::vector<int> quadratic_jacobian_positions_1_;
    std::vector<int> quadratic_jacobian_positions_2_;
    /** For each constraint, the position of each local variable of its tree. */
    std::vector<std::vector<int>> tree_jacobian_positions_;
    /** For each black-box constraint, the position of each variable. */
    std::vector<std::vector<int>> black_box_jacobian_positions_;

    /*
     * Hessian structure.
     */

    std::map<std::pair<int, int>, int> hessian_positions_;
    std::vector<Ipopt::Index> hessian_rows_;
    std::vector<Ipopt::Index> hessian_columns_;
    std::vector<int> objective_quadratic_hessian_positions_;
    std::vector<int> quadratic_hessian_positions_;
    /** Position of each entry of the local lower triangle of the objective tree Hessian. */
    std::vector<int> objective_tree_hessian_positions_;
    /** Same, for each constraint tree. */
    std::vector<std::vector<int>> constraints_trees_hessian_positions_;
    std::vector<double> tree_hessian_buffer_;

    /*
     * Evaluation cache.
     */

    std::vector<double> x_;
    bool trees_evaluated_ = false;
    bool objective_black_box_evaluated_ = false;
    BlackBoxFunctionOutput objective_black_box_output_;
    std::vector<bool> constraints_black_box_evaluated_;
    std::vector<BlackBoxFunctionOutput> constraints_black_box_outputs_;

};

MathOptTnlp::MathOptTnlp(
        const MathOptModel& model,
        IpoptModel& ipopt_model):
    model_(model),
    ipopt_model_(ipopt_model),
    sign_((model.objective_direction == ObjectiveDirection::Minimize)? 1.0: -1.0),
    exact_hessian_(!model.has_black_box()),
    constraints_trees_(model.number_of_constraints()),
    tree_jacobian_positions_(model.number_of_constraints()),
    black_box_jacobian_positions_(model.number_of_constraints()),
    constraints_trees_hessian_positions_(model.number_of_constraints()),
    x_(model.number_of_variables(), 0.0),
    constraints_black_box_evaluated_(model.number_of_constraints(), false),
    constraints_black_box_outputs_(model.number_of_constraints())
{
    int number_of_variables = model.number_of_variables();

    // Trees.
    if (!model.objective_nonlinear_elements_operators.empty()) {
        objective_tree_.reset(new ExpressionTree(
                    model.objective_nonlinear_elements_operators,
                    model.objective_nonlinear_elements_values,
                    model.objective_nonlinear_elements_variables,
                    model.objective_nonlinear_elements_number_of_children,
                    0,
                    (int)model.objective_nonlinear_elements_operators.size()));
    }
    if (!model.nonlinear_elements_operators.empty()) {
        for (int constraint_id = 0; constraint_id < model.number_of_constraints(); ++constraint_id) {
            int start = model.nonlinear_elements_constraints_starts[constraint_id];
            int end = model.nonlinear_constraint_end(constraint_id);
            if (end <= start)
                continue;
            constraints_trees_[constraint_id].reset(new ExpressionTree(
                        model.nonlinear_elements_operators,
                        model.nonlinear_elements_values,
                        model.nonlinear_elements_variables,
                        model.nonlinear_elements_number_of_children,
                        start,
                        end));
        }
    }

    // Jacobian structure: for each constraint, the union of the variables
    // of its terms, sorted.
    linear_jacobian_positions_.resize(model.elements_variables.size(), -1);
    quadratic_jacobian_positions_1_.resize(model.quadratic_elements_variables_1.size(), -1);
    quadratic_jacobian_positions_2_.resize(model.quadratic_elements_variables_1.size(), -1);
    for (int constraint_id = 0; constraint_id < model.number_of_constraints(); ++constraint_id) {
        std::map<int, int> positions;
        bool black_box = !model.constraints_functions.empty()
            && model.constraints_functions[constraint_id];
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

        for (int element_id = linear_start; element_id < linear_end; ++element_id)
            positions[model.elements_variables[element_id]] = -1;
        for (int element_id = quadratic_start; element_id < quadratic_end; ++element_id) {
            positions[model.quadratic_elements_variables_1[element_id]] = -1;
            positions[model.quadratic_elements_variables_2[element_id]] = -1;
        }
        if (constraints_trees_[constraint_id] != nullptr)
            for (int variable_id: constraints_trees_[constraint_id]->variables())
                positions[variable_id] = -1;
        if (black_box)
            for (int variable_id = 0; variable_id < number_of_variables; ++variable_id)
                positions[variable_id] = -1;

        for (auto& position: positions) {
            position.second = (int)jacobian_rows_.size();
            jacobian_rows_.push_back(constraint_id);
            jacobian_columns_.push_back(position.first);
        }

        for (int element_id = linear_start; element_id < linear_end; ++element_id)
            linear_jacobian_positions_[element_id] = positions[model.elements_variables[element_id]];
        for (int element_id = quadratic_start; element_id < quadratic_end; ++element_id) {
            quadratic_jacobian_positions_1_[element_id] = positions[model.quadratic_elements_variables_1[element_id]];
            quadratic_jacobian_positions_2_[element_id] = positions[model.quadratic_elements_variables_2[element_id]];
        }
        if (constraints_trees_[constraint_id] != nullptr)
            for (int variable_id: constraints_trees_[constraint_id]->variables())
                tree_jacobian_positions_[constraint_id].push_back(positions[variable_id]);
        if (black_box)
            for (int variable_id = 0; variable_id < number_of_variables; ++variable_id)
                black_box_jacobian_positions_[constraint_id].push_back(positions[variable_id]);
    }

    // Hessian structure.
    if (exact_hessian_) {
        for (int element_id = 0;
                element_id < (int)model.objective_quadratic_elements_variables_1.size();
                ++element_id) {
            objective_quadratic_hessian_positions_.push_back(hessian_position(
                        model.objective_quadratic_elements_variables_1[element_id],
                        model.objective_quadratic_elements_variables_2[element_id]));
        }
        for (int element_id = 0;
                element_id < (int)model.quadratic_elements_variables_1.size();
                ++element_id) {
            quadratic_hessian_positions_.push_back(hessian_position(
                        model.quadratic_elements_variables_1[element_id],
                        model.quadratic_elements_variables_2[element_id]));
        }
        auto tree_hessian_positions = [this](const ExpressionTree& tree)
        {
            const std::vector<int>& variables = tree.variables();
            std::vector<int> positions;
            for (int variable_a = 0; variable_a < (int)variables.size(); ++variable_a)
                for (int variable_b = 0; variable_b <= variable_a; ++variable_b)
                    positions.push_back(hessian_position(variables[variable_a], variables[variable_b]));
            return positions;
        };
        if (objective_tree_ != nullptr)
            objective_tree_hessian_positions_ = tree_hessian_positions(*objective_tree_);
        for (int constraint_id = 0; constraint_id < model.number_of_constraints(); ++constraint_id) {
            if (constraints_trees_[constraint_id] != nullptr) {
                constraints_trees_hessian_positions_[constraint_id]
                    = tree_hessian_positions(*constraints_trees_[constraint_id]);
            }
        }
    }
}

int MathOptTnlp::hessian_position(int variable_id_1, int variable_id_2)
{
    // Ipopt expects the lower triangle.
    std::pair<int, int> entry(
            (std::max)(variable_id_1, variable_id_2),
            (std::min)(variable_id_1, variable_id_2));
    auto it = hessian_positions_.find(entry);
    if (it != hessian_positions_.end())
        return it->second;
    int position = (int)hessian_rows_.size();
    hessian_positions_[entry] = position;
    hessian_rows_.push_back(entry.first);
    hessian_columns_.push_back(entry.second);
    return position;
}

void MathOptTnlp::update_x(const Ipopt::Number* x, bool new_x)
{
    if (!new_x)
        return;
    std::copy(x, x + x_.size(), x_.begin());
    trees_evaluated_ = false;
    objective_black_box_evaluated_ = false;
    std::fill(constraints_black_box_evaluated_.begin(), constraints_black_box_evaluated_.end(), false);
}

void MathOptTnlp::evaluate_trees()
{
    if (trees_evaluated_)
        return;
    if (objective_tree_ != nullptr)
        objective_tree_->evaluate(x_.data());
    for (const auto& tree: constraints_trees_)
        if (tree != nullptr)
            tree->evaluate(x_.data());
    trees_evaluated_ = true;
}

const BlackBoxFunctionOutput& MathOptTnlp::objective_black_box_output()
{
    if (!objective_black_box_evaluated_) {
        objective_black_box_output_ = model_.objective_function(x_);
        objective_black_box_evaluated_ = true;
    }
    return objective_black_box_output_;
}

const BlackBoxFunctionOutput& MathOptTnlp::constraint_black_box_output(int constraint_id)
{
    if (!constraints_black_box_evaluated_[constraint_id]) {
        constraints_black_box_outputs_[constraint_id] = model_.constraints_functions[constraint_id](x_);
        constraints_black_box_evaluated_[constraint_id] = true;
    }
    return constraints_black_box_outputs_[constraint_id];
}

bool MathOptTnlp::get_nlp_info(
        Ipopt::Index& n,
        Ipopt::Index& m,
        Ipopt::Index& nnz_jac_g,
        Ipopt::Index& nnz_h_lag,
        IndexStyleEnum& index_style)
{
    n = model_.number_of_variables();
    m = model_.number_of_constraints();
    nnz_jac_g = (Ipopt::Index)jacobian_rows_.size();
    nnz_h_lag = (Ipopt::Index)hessian_rows_.size();
    index_style = C_STYLE;
    return true;
}

bool MathOptTnlp::get_bounds_info(
        Ipopt::Index n,
        Ipopt::Number* x_l,
        Ipopt::Number* x_u,
        Ipopt::Index m,
        Ipopt::Number* g_l,
        Ipopt::Number* g_u)
{
    for (int variable_id = 0; variable_id < n; ++variable_id) {
        x_l[variable_id] = to_ipopt(model_.variables_lower_bounds[variable_id]);
        x_u[variable_id] = to_ipopt(model_.variables_upper_bounds[variable_id]);
    }
    for (int constraint_id = 0; constraint_id < m; ++constraint_id) {
        g_l[constraint_id] = to_ipopt(model_.constraints_lower_bounds[constraint_id]);
        g_u[constraint_id] = to_ipopt(model_.constraints_upper_bounds[constraint_id]);
    }
    return true;
}

bool MathOptTnlp::get_starting_point(
        Ipopt::Index n,
        bool init_x,
        Ipopt::Number* x,
        bool init_z,
        Ipopt::Number*,
        Ipopt::Number*,
        Ipopt::Index,
        bool init_lambda,
        Ipopt::Number*)
{
    if (init_z || init_lambda)
        return false;
    if (!init_x)
        return true;
    for (int variable_id = 0; variable_id < n; ++variable_id) {
        if (!model_.variables_initial_values.empty()) {
            x[variable_id] = model_.variables_initial_values[variable_id];
        } else {
            // 0, projected on the bounds.
            x[variable_id] = (std::min)(
                    (std::max)(0.0, model_.variables_lower_bounds[variable_id]),
                    model_.variables_upper_bounds[variable_id]);
        }
    }
    return true;
}

bool MathOptTnlp::eval_f(
        Ipopt::Index,
        const Ipopt::Number* x,
        bool new_x,
        Ipopt::Number& obj_value)
{
    update_x(x, new_x);
    double value = 0.0;
    if (!model_.objective_coefficients.empty())
        for (int variable_id = 0; variable_id < model_.number_of_variables(); ++variable_id)
            value += model_.objective_coefficients[variable_id] * x_[variable_id];
    for (int element_id = 0;
            element_id < (int)model_.objective_quadratic_elements_variables_1.size();
            ++element_id) {
        value += model_.objective_quadratic_elements_coefficients[element_id]
            * x_[model_.objective_quadratic_elements_variables_1[element_id]]
            * x_[model_.objective_quadratic_elements_variables_2[element_id]];
    }
    if (objective_tree_ != nullptr) {
        evaluate_trees();
        value += objective_tree_->value();
    }
    if (model_.objective_function)
        value += objective_black_box_output().objective_value;
    obj_value = sign_ * value;
    return true;
}

bool MathOptTnlp::eval_grad_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool new_x,
        Ipopt::Number* grad_f)
{
    update_x(x, new_x);
    for (int variable_id = 0; variable_id < n; ++variable_id) {
        grad_f[variable_id] = model_.objective_coefficients.empty()?
            0.0:
            model_.objective_coefficients[variable_id];
    }
    for (int element_id = 0;
            element_id < (int)model_.objective_quadratic_elements_variables_1.size();
            ++element_id) {
        int variable_id_1 = model_.objective_quadratic_elements_variables_1[element_id];
        int variable_id_2 = model_.objective_quadratic_elements_variables_2[element_id];
        double coefficient = model_.objective_quadratic_elements_coefficients[element_id];
        grad_f[variable_id_1] += coefficient * x_[variable_id_2];
        grad_f[variable_id_2] += coefficient * x_[variable_id_1];
    }
    if (objective_tree_ != nullptr) {
        evaluate_trees();
        const std::vector<int>& variables = objective_tree_->variables();
        const std::vector<double>& gradient = objective_tree_->gradient();
        for (int local_variable_id = 0; local_variable_id < (int)variables.size(); ++local_variable_id)
            grad_f[variables[local_variable_id]] += gradient[local_variable_id];
    }
    if (model_.objective_function) {
        const BlackBoxFunctionOutput& output = objective_black_box_output();
        for (int variable_id = 0; variable_id < n; ++variable_id)
            grad_f[variable_id] += output.gradient[variable_id];
    }
    for (int variable_id = 0; variable_id < n; ++variable_id)
        grad_f[variable_id] *= sign_;
    return true;
}

bool MathOptTnlp::eval_g(
        Ipopt::Index,
        const Ipopt::Number* x,
        bool new_x,
        Ipopt::Index m,
        Ipopt::Number* g)
{
    update_x(x, new_x);
    evaluate_trees();
    for (int constraint_id = 0; constraint_id < m; ++constraint_id) {
        double value = 0.0;
        if (!model_.elements_variables.empty()) {
            for (int element_id = model_.constraints_starts[constraint_id];
                    element_id < model_.constraint_end(constraint_id);
                    ++element_id) {
                value += model_.elements_coefficients[element_id]
                    * x_[model_.elements_variables[element_id]];
            }
        }
        if (!model_.quadratic_elements_variables_1.empty()) {
            for (int element_id = model_.quadratic_elements_constraints_starts[constraint_id];
                    element_id < model_.quadratic_constraint_end(constraint_id);
                    ++element_id) {
                value += model_.quadratic_elements_coefficients[element_id]
                    * x_[model_.quadratic_elements_variables_1[element_id]]
                    * x_[model_.quadratic_elements_variables_2[element_id]];
            }
        }
        if (constraints_trees_[constraint_id] != nullptr)
            value += constraints_trees_[constraint_id]->value();
        if (!black_box_jacobian_positions_[constraint_id].empty())
            value += constraint_black_box_output(constraint_id).objective_value;
        g[constraint_id] = value;
    }
    return true;
}

bool MathOptTnlp::eval_jac_g(
        Ipopt::Index,
        const Ipopt::Number* x,
        bool new_x,
        Ipopt::Index m,
        Ipopt::Index nele_jac,
        Ipopt::Index* iRow,
        Ipopt::Index* jCol,
        Ipopt::Number* values)
{
    if (values == NULL) {
        std::copy(jacobian_rows_.begin(), jacobian_rows_.end(), iRow);
        std::copy(jacobian_columns_.begin(), jacobian_columns_.end(), jCol);
        return true;
    }

    update_x(x, new_x);
    evaluate_trees();
    std::fill(values, values + nele_jac, 0.0);
    for (int element_id = 0; element_id < (int)model_.elements_variables.size(); ++element_id)
        values[linear_jacobian_positions_[element_id]] += model_.elements_coefficients[element_id];
    for (int element_id = 0;
            element_id < (int)model_.quadratic_elements_variables_1.size();
            ++element_id) {
        double coefficient = model_.quadratic_elements_coefficients[element_id];
        values[quadratic_jacobian_positions_1_[element_id]]
            += coefficient * x_[model_.quadratic_elements_variables_2[element_id]];
        values[quadratic_jacobian_positions_2_[element_id]]
            += coefficient * x_[model_.quadratic_elements_variables_1[element_id]];
    }
    for (int constraint_id = 0; constraint_id < m; ++constraint_id) {
        if (constraints_trees_[constraint_id] != nullptr) {
            const std::vector<double>& gradient = constraints_trees_[constraint_id]->gradient();
            const std::vector<int>& positions = tree_jacobian_positions_[constraint_id];
            for (int local_variable_id = 0; local_variable_id < (int)positions.size(); ++local_variable_id)
                values[positions[local_variable_id]] += gradient[local_variable_id];
        }
        if (!black_box_jacobian_positions_[constraint_id].empty()) {
            const BlackBoxFunctionOutput& output = constraint_black_box_output(constraint_id);
            const std::vector<int>& positions = black_box_jacobian_positions_[constraint_id];
            for (int variable_id = 0; variable_id < (int)positions.size(); ++variable_id)
                values[positions[variable_id]] += output.gradient[variable_id];
        }
    }
    return true;
}

bool MathOptTnlp::eval_h(
        Ipopt::Index,
        const Ipopt::Number* x,
        bool new_x,
        Ipopt::Number obj_factor,
        Ipopt::Index m,
        const Ipopt::Number* lambda,
        bool,
        Ipopt::Index nele_hess,
        Ipopt::Index* iRow,
        Ipopt::Index* jCol,
        Ipopt::Number* values)
{
    if (!exact_hessian_)
        return false;

    if (values == NULL) {
        std::copy(hessian_rows_.begin(), hessian_rows_.end(), iRow);
        std::copy(hessian_columns_.begin(), hessian_columns_.end(), jCol);
        return true;
    }

    update_x(x, new_x);
    evaluate_trees();
    std::fill(values, values + nele_hess, 0.0);

    // Quadratic terms: the Hessian of c x_i x_j has c at (i, j) and (j, i),
    // or 2 c at (i, i) if i = j.
    double objective_factor = sign_ * obj_factor;
    for (int element_id = 0;
            element_id < (int)model_.objective_quadratic_elements_variables_1.size();
            ++element_id) {
        double coefficient = model_.objective_quadratic_elements_coefficients[element_id];
        if (model_.objective_quadratic_elements_variables_1[element_id]
                == model_.objective_quadratic_elements_variables_2[element_id]) {
            coefficient *= 2;
        }
        values[objective_quadratic_hessian_positions_[element_id]] += objective_factor * coefficient;
    }
    if (!model_.quadratic_elements_variables_1.empty()) {
        for (int constraint_id = 0; constraint_id < m; ++constraint_id) {
            for (int element_id = model_.quadratic_elements_constraints_starts[constraint_id];
                    element_id < model_.quadratic_constraint_end(constraint_id);
                    ++element_id) {
                double coefficient = model_.quadratic_elements_coefficients[element_id];
                if (model_.quadratic_elements_variables_1[element_id]
                        == model_.quadratic_elements_variables_2[element_id]) {
                    coefficient *= 2;
                }
                values[quadratic_hessian_positions_[element_id]] += lambda[constraint_id] * coefficient;
            }
        }
    }

    // Trees.
    auto add_tree_hessian = [this, values](
            ExpressionTree& tree,
            const std::vector<int>& positions,
            double weight)
    {
        if (weight == 0.0)
            return;
        tree_hessian_buffer_.assign(positions.size(), 0.0);
        tree.add_hessian(weight, tree_hessian_buffer_);
        for (int pos = 0; pos < (int)positions.size(); ++pos)
            values[positions[pos]] += tree_hessian_buffer_[pos];
    };
    if (objective_tree_ != nullptr)
        add_tree_hessian(*objective_tree_, objective_tree_hessian_positions_, objective_factor);
    for (int constraint_id = 0; constraint_id < m; ++constraint_id) {
        if (constraints_trees_[constraint_id] != nullptr) {
            add_tree_hessian(
                    *constraints_trees_[constraint_id],
                    constraints_trees_hessian_positions_[constraint_id],
                    lambda[constraint_id]);
        }
    }
    return true;
}

void MathOptTnlp::finalize_solution(
        Ipopt::SolverReturn,
        Ipopt::Index n,
        const Ipopt::Number* x,
        const Ipopt::Number*,
        const Ipopt::Number*,
        Ipopt::Index,
        const Ipopt::Number*,
        const Ipopt::Number*,
        Ipopt::Number obj_value,
        const Ipopt::IpoptData*,
        Ipopt::IpoptCalculatedQuantities*)
{
    ipopt_model_.solution.assign(x, x + n);
    ipopt_model_.objective_value = sign_ * obj_value;
}

}

IpoptModel::IpoptModel():
    application(IpoptApplicationFactory())
{
    Ipopt::ApplicationReturnStatus status = application->Initialize();
    if (status != Ipopt::Solve_Succeeded) {
        throw std::runtime_error(
                FUNC_SIGNATURE + ": "
                "IpoptApplication::Initialize failed with status "
                + std::to_string((int)status) + ".");
    }
}

void mathoptsolverscmake::load(
        IpoptModel& ipopt_model,
        const MathOptModel& model)
{
    if (model.has_non_continuous_variables()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "Ipopt does not support integer variables.");
    }
    MathOptTnlp* problem = new MathOptTnlp(model, ipopt_model);
    ipopt_model.problem = problem;
    if (!problem->exact_hessian()) {
        ipopt_model.application->Options()->SetStringValue(
                "hessian_approximation", "limited-memory");
    }
}

void mathoptsolverscmake::reduce_printout(
        IpoptModel& ipopt_model)
{
    ipopt_model.application->Options()->SetIntegerValue("print_level", 0);
    ipopt_model.application->Options()->SetStringValue("sb", "yes");
}

void mathoptsolverscmake::set_time_limit(
        IpoptModel& ipopt_model,
        double time_limit)
{
    ipopt_model.application->Options()->SetNumericValue("max_wall_time", time_limit);
}

void mathoptsolverscmake::set_log_file(
        IpoptModel& ipopt_model,
        const std::string& log_file)
{
    if (!ipopt_model.application->OpenOutputFile(log_file, Ipopt::J_ITERSUMMARY)) {
        throw std::runtime_error(
                FUNC_SIGNATURE + ": "
                "unable to open file \"" + log_file + "\".");
    }
}

void mathoptsolverscmake::solve(
        IpoptModel& ipopt_model)
{
    if (!Ipopt::IsValid(ipopt_model.problem)) {
        throw std::logic_error(
                FUNC_SIGNATURE + ": "
                "load() must be called before solve().");
    }
    ipopt_model.status = ipopt_model.application->OptimizeTNLP(ipopt_model.problem);
}

double mathoptsolverscmake::get_solution_value(
        const IpoptModel& ipopt_model)
{
    return ipopt_model.objective_value;
}

std::vector<double> mathoptsolverscmake::get_solution(
        const IpoptModel& ipopt_model)
{
    return ipopt_model.solution;
}

bool mathoptsolverscmake::is_infeasible(
        const IpoptModel& ipopt_model)
{
    return ipopt_model.status == Ipopt::Infeasible_Problem_Detected;
}

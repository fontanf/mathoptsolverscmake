#include "mathoptsolverscmake/mathopt_build.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>

using namespace mathoptsolverscmake;

namespace
{

void invalid(const std::string& message)
{
    throw std::invalid_argument(message);
}

/** Check that a model array of size 'current' can grow by 'added' entries. */
void check_capacity(
        std::size_t current,
        std::size_t added,
        const char* what)
{
    if (added > (std::size_t)std::numeric_limits<int>::max() - current)
        invalid(std::string("too many ") + what + ": MathOptModel indexes them with int.");
}

/** Check that 'array' is given if the batch has entries. */
void check_given(
        const void* array,
        std::size_t size,
        const char* name)
{
    if (array == nullptr && size > 0)
        invalid(std::string(name) + " is required.");
}

/** Check the CSR 'starts' of 'count' rows over 'number_of_elements' elements. */
void check_starts(
        const int* starts,
        std::size_t count,
        std::size_t number_of_elements,
        const char* name,
        const char* elements_name)
{
    if (starts == nullptr) {
        if (number_of_elements > 0)
            invalid(std::string(name) + " is required with " + elements_name + ".");
        return;
    }
    for (std::size_t i = 0; i < count; ++i) {
        int start = starts[i];
        if ((i == 0 && start != 0) || (i > 0 && start < starts[i - 1])
                || start < 0 || (std::size_t)start > number_of_elements) {
            invalid(std::string(name) + "[" + std::to_string(i) + "] = " + std::to_string(start)
                    + ": starts must begin at 0, be non-decreasing and not exceed the size of "
                    + elements_name + " (" + std::to_string(number_of_elements) + ").");
        }
    }
}

/** End of row 'i' of CSR 'starts' with 'count' rows over 'size' elements. */
std::size_t row_end(
        const int* starts,
        std::size_t count,
        std::size_t size,
        std::size_t i)
{
    return (i + 1 < count)? (std::size_t)starts[i + 1]: size;
}

void check_variables(
        const int* variables,
        std::size_t size,
        int number_of_variables,
        const char* name)
{
    for (std::size_t i = 0; i < size; ++i) {
        if (variables[i] < 0 || variables[i] >= number_of_variables) {
            invalid(std::string(name) + "[" + std::to_string(i) + "] = " + std::to_string(variables[i])
                    + " is not a variable index (the model has " + std::to_string(number_of_variables)
                    + " variables).");
        }
    }
}

/** The number of children an expression-tree operator accepts. */
bool operator_arity(
        char op,
        int& min,
        int& max)
{
    switch (op) {
    case 'v': case 'd': case 'k':
        min = max = 0; return true;
    case 'n': case 'e': case 'l': case 'g': case 'q':
    case 's': case 'c': case 't': case 'h': case 'a':
        min = max = 1; return true;
    case '-': case '/': case 'p':
        min = max = 2; return true;
    case '+': case '*':
        min = 2; max = std::numeric_limits<int>::max(); return true;
    case 'm': case 'M':
        min = 1; max = std::numeric_limits<int>::max(); return true;
    default:
        return false;
    }
}

void check_trees_arrays(
        const ExpressionTreesBatch& trees,
        const std::string& prefix)
{
    check_given(trees.operators, trees.number_of_nodes, (prefix + "operators").c_str());
    check_given(trees.number_of_children, trees.number_of_nodes,
            (prefix + "number_of_children").c_str());
}

/**
 * Check that the nodes [begin, end) of 'trees' form exactly one well-formed
 * pre-order tree (or nothing, if 'allow_empty'), whose variables are below
 * 'number_of_variables' and defined variables below
 * 'number_of_defined_variables'.
 */
void check_tree(
        const ExpressionTreesBatch& trees,
        std::size_t begin,
        std::size_t end,
        int number_of_variables,
        int number_of_defined_variables,
        bool allow_empty,
        const std::string& what)
{
    if (begin == end) {
        if (!allow_empty)
            invalid(what + " has an empty expression tree.");
        return;
    }
    long long pending = 1;
    for (std::size_t i = begin; i < end; ++i) {
        if (pending == 0)
            invalid(what + ": nodes " + std::to_string(i) + " to " + std::to_string(end - 1)
                    + " are not part of its expression tree (a single root is expected).");
        --pending;
        char op = trees.operators[i];
        int number_of_children = trees.number_of_children[i];
        int min = 0;
        int max = 0;
        if (!operator_arity(op, min, max))
            invalid(what + ": node " + std::to_string(i) + " has unknown operator code "
                    + std::to_string((int)(unsigned char)op) + ".");
        if (number_of_children < min || number_of_children > max)
            invalid(what + ": node " + std::to_string(i) + " (operator '" + std::string(1, op)
                    + "') has " + std::to_string(number_of_children) + " children.");
        if (op == 'v' || op == 'd') {
            int index = (trees.variables != nullptr)? trees.variables[i]: -1;
            int limit = (op == 'v')? number_of_variables: number_of_defined_variables;
            if (index < 0 || index >= limit)
                invalid(what + ": node " + std::to_string(i) + " references "
                        + (op == 'v'? "variable ": "defined variable ") + std::to_string(index)
                        + ", which is out of range [0, " + std::to_string(limit) + ").");
        }
        pending += number_of_children;
    }
    if (pending != 0)
        invalid(what + ": its expression tree is missing " + std::to_string(pending) + " node(s).");
}

/** Append all the nodes of 'trees' to the given model arrays. */
void append_trees(
        const ExpressionTreesBatch& trees,
        std::vector<char>& operators,
        std::vector<double>& values,
        std::vector<int>& variables,
        std::vector<int>& number_of_children)
{
    std::size_t n = trees.number_of_nodes;
    if (n == 0)
        return;
    operators.insert(operators.end(), trees.operators, trees.operators + n);
    number_of_children.insert(number_of_children.end(),
            trees.number_of_children, trees.number_of_children + n);
    if (trees.values != nullptr)
        values.insert(values.end(), trees.values, trees.values + n);
    else
        values.resize(values.size() + n, 0.0);
    if (trees.variables != nullptr)
        variables.insert(variables.end(), trees.variables, trees.variables + n);
    else
        variables.resize(variables.size() + n, -1);
}

/** Append 'count' names, giving the existing entries default names if needed. */
void append_names(
        std::vector<std::string>& names,
        int current_count,
        std::size_t count,
        const std::vector<std::string>* new_names,
        const char* prefix)
{
    if (new_names == nullptr && names.empty())
        return;
    for (int i = (int)names.size(); i < current_count; ++i)
        names.push_back(prefix + std::to_string(i));
    for (std::size_t i = 0; i < count; ++i) {
        if (new_names != nullptr)
            names.push_back((*new_names)[i]);
        else
            names.push_back(prefix + std::to_string(current_count + i));
    }
}

/**
 * Make room for 'size' entries, keeping the geometric growth of std::vector:
 * reserving exactly 'size' on every batch would reallocate (and copy) the
 * whole vector each time, making a model built one entry at a time quadratic.
 */
template <typename T>
void grow(
        std::vector<T>& vector,
        std::size_t size)
{
    if (size > vector.capacity())
        vector.reserve(std::max(size, 2 * vector.capacity()));
}

template <typename T>
void append(
        std::vector<T>& vector,
        const T* data,
        std::size_t size)
{
    if (size > 0)
        vector.insert(vector.end(), data, data + size);
}

}

int mathoptsolverscmake::add_variables(
        MathOptModel& model,
        const VariablesBatch& batch)
{
    std::size_t n = batch.number;
    int first = model.number_of_variables();
    check_capacity(first, n, "variables");
    if (batch.names != nullptr && batch.names->size() != n)
        invalid("names has " + std::to_string(batch.names->size())
                + " entries, but the batch has " + std::to_string(n) + " variables.");
    for (std::size_t i = 0; i < n; ++i) {
        int type = (batch.types != nullptr)? batch.types[i]: (int)batch.type;
        if (type < 0 || type > (int)VariableType::Integer)
            invalid("types[" + std::to_string(i) + "] = " + std::to_string(type)
                    + " is not a VariableType value.");
    }

    grow(model.variables_lower_bounds, first + n);
    grow(model.variables_upper_bounds, first + n);
    grow(model.variables_types, first + n);
    grow(model.objective_coefficients, first + n);
    for (std::size_t i = 0; i < n; ++i) {
        model.variables_lower_bounds.push_back(
                (batch.lower_bounds != nullptr)? batch.lower_bounds[i]: batch.lower_bound);
        model.variables_upper_bounds.push_back(
                (batch.upper_bounds != nullptr)? batch.upper_bounds[i]: batch.upper_bound);
        model.variables_types.push_back(
                (batch.types != nullptr)? static_cast<VariableType>(batch.types[i]): batch.type);
        model.objective_coefficients.push_back(
                (batch.objective_coefficients != nullptr)?
                batch.objective_coefficients[i]: batch.objective_coefficient);
    }
    if (batch.initial_values != nullptr || !model.variables_initial_values.empty()) {
        model.variables_initial_values.resize(first, 0.0);
        if (batch.initial_values != nullptr)
            append(model.variables_initial_values, batch.initial_values, n);
        else
            model.variables_initial_values.resize(first + n, 0.0);
    }
    append_names(model.variables_names, first, n, batch.names, "x");
    return first;
}

int mathoptsolverscmake::add_constraints(
        MathOptModel& model,
        const ConstraintsBatch& batch)
{
    std::size_t m = batch.number;
    const ExpressionTreesBatch& trees = batch.nonlinear;
    int n = model.number_of_variables();
    int first = model.number_of_constraints();

    // Validate everything before modifying the model.
    check_capacity(first, m, "constraints");
    check_capacity(model.number_of_elements(), batch.number_of_elements, "linear terms");
    check_capacity(model.quadratic_elements_variables_1.size(),
            batch.number_of_quadratic_elements, "quadratic terms");
    check_capacity(model.nonlinear_elements_operators.size(), trees.number_of_nodes,
            "nonlinear nodes");
    if (batch.names != nullptr && batch.names->size() != m)
        invalid("names has " + std::to_string(batch.names->size())
                + " entries, but the batch has " + std::to_string(m) + " constraints.");

    check_given(batch.variables, batch.number_of_elements, "variables");
    check_given(batch.coefficients, batch.number_of_elements, "coefficients");
    check_starts(batch.starts, m, batch.number_of_elements, "starts", "variables");
    check_variables(batch.variables, batch.number_of_elements, n, "variables");

    check_given(batch.quadratic_variables_1, batch.number_of_quadratic_elements,
            "quadratic_variables_1");
    check_given(batch.quadratic_variables_2, batch.number_of_quadratic_elements,
            "quadratic_variables_2");
    check_given(batch.quadratic_coefficients, batch.number_of_quadratic_elements,
            "quadratic_coefficients");
    check_starts(batch.quadratic_starts, m, batch.number_of_quadratic_elements,
            "quadratic_starts", "quadratic_variables_1");
    check_variables(batch.quadratic_variables_1, batch.number_of_quadratic_elements, n,
            "quadratic_variables_1");
    check_variables(batch.quadratic_variables_2, batch.number_of_quadratic_elements, n,
            "quadratic_variables_2");

    check_trees_arrays(trees, "nonlinear_");
    check_starts(batch.nonlinear_starts, m, trees.number_of_nodes,
            "nonlinear_starts", "nonlinear_operators");
    if (trees.number_of_nodes > 0) {
        for (std::size_t i = 0; i < m; ++i) {
            check_tree(trees,
                    batch.nonlinear_starts[i],
                    row_end(batch.nonlinear_starts, m, trees.number_of_nodes, i),
                    n, model.number_of_defined_variables(), true,
                    "constraint " + std::to_string(i) + " of the batch");
        }
    }

    // Append.
    int linear_offset = model.number_of_elements();
    int quadratic_offset = (int)model.quadratic_elements_variables_1.size();
    int nonlinear_offset = (int)model.nonlinear_elements_operators.size();
    grow(model.constraints_lower_bounds, first + m);
    grow(model.constraints_upper_bounds, first + m);
    grow(model.constraints_starts, first + m);
    for (std::size_t i = 0; i < m; ++i) {
        model.constraints_lower_bounds.push_back(
                (batch.lower_bounds != nullptr)? batch.lower_bounds[i]: batch.lower_bound);
        model.constraints_upper_bounds.push_back(
                (batch.upper_bounds != nullptr)? batch.upper_bounds[i]: batch.upper_bound);
        model.constraints_starts.push_back(
                linear_offset + ((batch.starts != nullptr)? batch.starts[i]: 0));
    }
    append(model.elements_variables, batch.variables, batch.number_of_elements);
    append(model.elements_coefficients, batch.coefficients, batch.number_of_elements);

    model.quadratic_elements_constraints_starts.resize(first, quadratic_offset);
    grow(model.quadratic_elements_constraints_starts, first + m);
    for (std::size_t i = 0; i < m; ++i) {
        model.quadratic_elements_constraints_starts.push_back(
                quadratic_offset + ((batch.quadratic_starts != nullptr)? batch.quadratic_starts[i]: 0));
    }
    append(model.quadratic_elements_variables_1, batch.quadratic_variables_1,
            batch.number_of_quadratic_elements);
    append(model.quadratic_elements_variables_2, batch.quadratic_variables_2,
            batch.number_of_quadratic_elements);
    append(model.quadratic_elements_coefficients, batch.quadratic_coefficients,
            batch.number_of_quadratic_elements);

    // The nonlinear starts are only stored once the model has expression trees.
    if (trees.number_of_nodes > 0 || !model.nonlinear_elements_constraints_starts.empty()) {
        model.nonlinear_elements_constraints_starts.resize(first, nonlinear_offset);
        grow(model.nonlinear_elements_constraints_starts, first + m);
        for (std::size_t i = 0; i < m; ++i) {
            model.nonlinear_elements_constraints_starts.push_back(
                    nonlinear_offset + ((batch.nonlinear_starts != nullptr)? batch.nonlinear_starts[i]: 0));
        }
        append_trees(trees,
                model.nonlinear_elements_operators,
                model.nonlinear_elements_values,
                model.nonlinear_elements_variables,
                model.nonlinear_elements_number_of_children);
    }

    if (!model.constraints_functions.empty())
        model.constraints_functions.resize(first + m);
    append_names(model.constraints_names, first, m, batch.names, "c");
    return first;
}

int mathoptsolverscmake::add_defined_variables(
        MathOptModel& model,
        const DefinedVariablesBatch& batch)
{
    std::size_t count = batch.number;
    const ExpressionTreesBatch& trees = batch.trees;
    int n = model.number_of_variables();
    int first = model.number_of_defined_variables();

    check_capacity(first, count, "defined variables");
    check_capacity(model.defined_variables_elements_operators.size(), trees.number_of_nodes,
            "defined-variable nodes");
    check_trees_arrays(trees, "");
    if (count == 0 && trees.number_of_nodes > 0)
        invalid("the batch has no defined variables, but operators isn't empty.");
    check_given(batch.starts, count, "starts");
    check_starts(batch.starts, count, trees.number_of_nodes, "starts", "operators");
    for (std::size_t i = 0; i < count; ++i) {
        // A defined variable may only reference defined variables with a lower index.
        check_tree(trees,
                batch.starts[i],
                row_end(batch.starts, count, trees.number_of_nodes, i),
                n, first + (int)i, false,
                "defined variable " + std::to_string(i) + " of the batch (index "
                + std::to_string(first + i) + ", which may only reference lower indices)");
    }

    int offset = (int)model.defined_variables_elements_operators.size();
    for (std::size_t i = 0; i < count; ++i)
        model.defined_variables_starts.push_back(offset + batch.starts[i]);
    append_trees(trees,
            model.defined_variables_elements_operators,
            model.defined_variables_elements_values,
            model.defined_variables_elements_variables,
            model.defined_variables_elements_number_of_children);
    return first;
}

void mathoptsolverscmake::set_objective(
        MathOptModel& model,
        const ObjectiveBatch& batch)
{
    int n = model.number_of_variables();

    // Validate every part before modifying the model.
    if (batch.set_quadratic) {
        std::size_t size = batch.number_of_quadratic_elements;
        check_capacity(0, size, "quadratic terms");
        check_given(batch.quadratic_variables_1, size, "quadratic_variables_1");
        check_given(batch.quadratic_variables_2, size, "quadratic_variables_2");
        check_given(batch.quadratic_coefficients, size, "quadratic_coefficients");
        check_variables(batch.quadratic_variables_1, size, n, "quadratic_variables_1");
        check_variables(batch.quadratic_variables_2, size, n, "quadratic_variables_2");
    }
    if (batch.set_nonlinear) {
        check_capacity(0, batch.nonlinear.number_of_nodes, "nonlinear nodes");
        check_trees_arrays(batch.nonlinear, "nonlinear_");
        check_tree(batch.nonlinear, 0, batch.nonlinear.number_of_nodes,
                n, model.number_of_defined_variables(), true, "the objective");
    }

    if (batch.set_direction)
        model.objective_direction = batch.direction;
    if (batch.set_coefficients) {
        model.objective_coefficients.resize(n);
        for (int variable_id = 0; variable_id < n; ++variable_id) {
            model.objective_coefficients[variable_id] = (batch.coefficients != nullptr)?
                batch.coefficients[variable_id]: batch.coefficient;
        }
    }
    if (batch.set_quadratic) {
        std::size_t size = batch.number_of_quadratic_elements;
        model.objective_quadratic_elements_variables_1.clear();
        model.objective_quadratic_elements_variables_2.clear();
        model.objective_quadratic_elements_coefficients.clear();
        append(model.objective_quadratic_elements_variables_1, batch.quadratic_variables_1, size);
        append(model.objective_quadratic_elements_variables_2, batch.quadratic_variables_2, size);
        append(model.objective_quadratic_elements_coefficients, batch.quadratic_coefficients, size);
    }
    if (batch.set_nonlinear) {
        model.objective_nonlinear_elements_operators.clear();
        model.objective_nonlinear_elements_values.clear();
        model.objective_nonlinear_elements_variables.clear();
        model.objective_nonlinear_elements_number_of_children.clear();
        append_trees(batch.nonlinear,
                model.objective_nonlinear_elements_operators,
                model.objective_nonlinear_elements_values,
                model.objective_nonlinear_elements_variables,
                model.objective_nonlinear_elements_number_of_children);
    }
}

/**
 * Python bindings for the solver-independent part of MathOptSolversCMake:
 * MathOptModel and the MPS / LP / JSON / .nl readers and writers.
 *
 * MathOptModel's arrays are exposed as read-only NumPy views onto the C++
 * vectors (no copy); nanobind's reference_internal policy (the default for
 * properties) makes each view keep its model alive. Names are exposed as
 * lists of str, and expression-tree operators as bytes. Readers take the
 * file contents as bytes (see read() in __init__.py, which also handles
 * .gz/.bz2 decompression) and release the GIL while parsing.
 *
 * Models can also be built from Python, in batches (see "Model building"
 * below). Growing a model may reallocate its vectors, so views obtained
 * before a modification must not be used after it.
 */
#include "mathoptsolverscmake/mathopt.hpp"
#include "mathoptsolverscmake/mathopt_mps.hpp"
#include "mathoptsolverscmake/mathopt_lp.hpp"
#include "mathoptsolverscmake/mathopt_json.hpp"
#include "mathoptsolverscmake/mathopt_nl.hpp"
#include "mathoptsolverscmake/dual.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>

#include <cstdint>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <variant>

namespace nb = nanobind;
using namespace nb::literals;
using namespace mathoptsolverscmake;

namespace
{

template <typename T>
using ReadOnlyArray = nb::ndarray<nb::numpy, const T, nb::ndim<1>>;

/** A read-only NumPy view onto one of the model's vectors. */
template <typename T, std::vector<T> MathOptModel::*Member>
ReadOnlyArray<T> view_member(const MathOptModel& model)
{
    const std::vector<T>& vector = model.*Member;
    return ReadOnlyArray<T>(vector.data(), {vector.size()});
}

/** Variable types as a NumPy array of their integer values (see VariableType). */
nb::ndarray<nb::numpy, int8_t, nb::ndim<1>> variables_types_array(const MathOptModel& model)
{
    size_t n = model.variables_types.size();
    int8_t* data = new int8_t[n];
    for (size_t i = 0; i < n; ++i)
        data[i] = static_cast<int8_t>(model.variables_types[i]);
    nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<int8_t*>(p); });
    return nb::ndarray<nb::numpy, int8_t, nb::ndim<1>>(data, {n}, owner);
}

/** Constraint classes as a NumPy array of their integer values (see ConstraintClass). */
nb::ndarray<nb::numpy, int8_t, nb::ndim<1>> constraints_classes_array(const MathOptModel& model)
{
    size_t n = model.number_of_constraints();
    int8_t* data = new int8_t[n];
    for (size_t i = 0; i < n; ++i)
        data[i] = static_cast<int8_t>(model.constraint_class(i));
    nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<int8_t*>(p); });
    return nb::ndarray<nb::numpy, int8_t, nb::ndim<1>>(data, {n}, owner);
}

nb::bytes to_bytes(const std::vector<char>& operators)
{
    return nb::bytes(operators.data(), operators.size());
}

/** Parse `data` with `reader`, without holding the GIL. */
MathOptModel read_from_bytes(
        MathOptModel (*reader)(std::istream&),
        nb::bytes data)
{
    std::string contents(data.c_str(), data.size());
    nb::gil_scoped_release release;
    std::istringstream in(contents, std::ios::in | std::ios::binary);
    return reader(in);
}

std::string write_to_string(
        void (*writer)(const MathOptModel&, std::ostream&),
        const MathOptModel& model)
{
    std::ostringstream out;
    writer(model, out);
    return out.str();
}

/*
 * Model building.
 *
 * The add_* / set_objective methods take whole batches as NumPy arrays laid
 * out like MathOptModel's own fields (CSR rows, pre-order expression trees),
 * so that building a model costs a few array copies rather than one Python
 * call per variable, constraint or term. Arrays of another dtype (or Python
 * lists) are converted once by nanobind. Each batch is fully validated
 * before the model is modified, so a rejected batch leaves it unchanged; the
 * validation and the copies run without the GIL.
 */

template <typename T>
using InputArray = nb::ndarray<const T, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

/** A borrowed, possibly absent, 1-D input array. */
template <typename T>
struct Input
{
    const T* data = nullptr;
    size_t size = 0;
    bool given = false;

    const T& operator[](size_t i) const { return data[i]; }
};

template <typename T>
Input<T> input(const std::optional<InputArray<T>>& array)
{
    Input<T> res;
    if (array) {
        res.data = array->data();
        res.size = array->shape(0);
        res.given = true;
    }
    return res;
}

/** An input that is either a scalar (applied to every entry) or an array. */
template <typename T>
struct ScalarOrInput
{
    Input<T> array;
    T scalar;

    T operator[](size_t i) const { return array.given? array.data[i]: scalar; }
};

template <typename T, typename S>
ScalarOrInput<T> scalar_or_input(
        const std::optional<std::variant<S, InputArray<T>>>& value,
        T default_value)
{
    ScalarOrInput<T> res;
    res.scalar = default_value;
    if (!value)
        return res;
    if (const S* scalar = std::get_if<S>(&*value)) {
        res.scalar = static_cast<T>(*scalar);
    } else {
        const InputArray<T>& array = std::get<InputArray<T>>(*value);
        res.array.data = array.data();
        res.array.size = array.shape(0);
        res.array.given = true;
    }
    return res;
}

/** Operator char codes, given as bytes or as an int8/uint8 array. */
using OperatorsInput = std::variant<nb::bytes, InputArray<uint8_t>, InputArray<int8_t>>;

Input<char> operators_input(const std::optional<OperatorsInput>& operators)
{
    Input<char> res;
    if (!operators)
        return res;
    res.given = true;
    if (const nb::bytes* bytes = std::get_if<nb::bytes>(&*operators)) {
        res.data = bytes->c_str();
        res.size = bytes->size();
    } else if (const auto* array = std::get_if<InputArray<uint8_t>>(&*operators)) {
        res.data = reinterpret_cast<const char*>(array->data());
        res.size = array->shape(0);
    } else {
        const auto& array_i8 = std::get<InputArray<int8_t>>(*operators);
        res.data = reinterpret_cast<const char*>(array_i8.data());
        res.size = array_i8.shape(0);
    }
    return res;
}

[[noreturn]] void invalid(const std::string& message)
{
    throw std::invalid_argument(message);
}

/** Merge the batch size implied by one argument into `count`. */
void merge_count(
        std::optional<size_t>& count,
        size_t size,
        const char* name)
{
    if (!count) {
        count = size;
    } else if (*count != size) {
        invalid(std::string(name) + " has " + std::to_string(size)
                + " entries, but " + std::to_string(*count) + " were expected.");
    }
}

template <typename T>
void merge_count(std::optional<size_t>& count, const Input<T>& in, const char* name)
{
    if (in.given)
        merge_count(count, in.size, name);
}

template <typename T>
void merge_count(std::optional<size_t>& count, const ScalarOrInput<T>& in, const char* name)
{
    merge_count(count, in.array, name);
}

template <typename T>
void check_size(const Input<T>& in, size_t size, const std::string& name, const std::string& reference)
{
    if (in.given && in.size != size) {
        invalid(name + " has " + std::to_string(in.size)
                + " entries, but " + reference + " has " + std::to_string(size) + ".");
    }
}

/** Check that the model's arrays of the given kind can grow by `added`. */
void check_capacity(size_t current, size_t added, const char* what)
{
    if (current + added > (size_t)std::numeric_limits<int>::max())
        invalid(std::string("too many ") + what + ": MathOptModel indexes them with int.");
}

/** Check the CSR `starts` of `count` rows over `number_of_elements` elements. */
void check_starts(
        const Input<int32_t>& starts,
        size_t count,
        size_t number_of_elements,
        const char* name,
        const char* elements_name)
{
    if (!starts.given) {
        if (number_of_elements > 0)
            invalid(std::string(name) + " is required with " + elements_name + ".");
        return;
    }
    if (starts.size != count) {
        invalid(std::string(name) + " has " + std::to_string(starts.size)
                + " entries, but there are " + std::to_string(count) + " rows.");
    }
    for (size_t i = 0; i < count; ++i) {
        int32_t start = starts[i];
        if ((i == 0 && start != 0) || (i > 0 && start < starts[i - 1])
                || start < 0 || (size_t)start > number_of_elements) {
            invalid(std::string(name) + "[" + std::to_string(i) + "] = " + std::to_string(start)
                    + ": starts must begin at 0, be non-decreasing and not exceed the size of "
                    + elements_name + " (" + std::to_string(number_of_elements) + ").");
        }
    }
}

void check_variables(
        const Input<int32_t>& variables,
        int number_of_variables,
        const char* name)
{
    for (size_t i = 0; i < variables.size; ++i) {
        if (variables[i] < 0 || variables[i] >= number_of_variables) {
            invalid(std::string(name) + "[" + std::to_string(i) + "] = " + std::to_string(variables[i])
                    + " is not a variable index (the model has " + std::to_string(number_of_variables)
                    + " variables).");
        }
    }
}

/** The number of children an expression-tree operator accepts. */
bool operator_arity(char op, int& min, int& max)
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

/** An expression-tree batch: operators and the arrays parallel to them. */
struct TreesInput
{
    Input<char> operators;
    Input<double> values;
    Input<int32_t> variables;
    Input<int32_t> number_of_children;

    size_t size() const { return operators.size; }
    bool given() const { return operators.given; }
};

TreesInput trees_input(
        const std::optional<OperatorsInput>& operators,
        const std::optional<InputArray<double>>& values,
        const std::optional<InputArray<int32_t>>& variables,
        const std::optional<InputArray<int32_t>>& number_of_children,
        const std::string& prefix)
{
    TreesInput res;
    res.operators = operators_input(operators);
    res.values = input(values);
    res.variables = input(variables);
    res.number_of_children = input(number_of_children);
    std::string operators_name = prefix + "operators";
    if (!res.operators.given) {
        if (res.values.given || res.variables.given || res.number_of_children.given)
            invalid(operators_name + " is required with the other " + prefix + "* arrays.");
        return res;
    }
    if (!res.number_of_children.given && res.operators.size > 0)
        invalid(prefix + "number_of_children is required with " + operators_name + ".");
    check_size(res.number_of_children, res.size(), prefix + "number_of_children", operators_name);
    check_size(res.values, res.size(), prefix + "values", operators_name);
    check_size(res.variables, res.size(), prefix + "variables", operators_name);
    return res;
}

/**
 * Check that the nodes [begin, end) of `trees` form exactly one well-formed
 * pre-order tree (or nothing, if `allow_empty`), whose variables are below
 * `number_of_variables` and defined variables below `number_of_defined_variables`.
 */
void check_tree(
        const TreesInput& trees,
        size_t begin,
        size_t end,
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
    for (size_t i = begin; i < end; ++i) {
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
            int index = trees.variables.given? trees.variables[i]: -1;
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

/** Append all the nodes of `trees` to the given model arrays. */
void append_trees(
        const TreesInput& trees,
        std::vector<char>& operators,
        std::vector<double>& values,
        std::vector<int>& variables,
        std::vector<int>& number_of_children)
{
    size_t n = trees.size();
    operators.insert(operators.end(), trees.operators.data, trees.operators.data + n);
    number_of_children.insert(number_of_children.end(),
            trees.number_of_children.data, trees.number_of_children.data + n);
    if (trees.values.given)
        values.insert(values.end(), trees.values.data, trees.values.data + n);
    else
        values.resize(values.size() + n, 0.0);
    if (trees.variables.given)
        variables.insert(variables.end(), trees.variables.data, trees.variables.data + n);
    else
        variables.resize(variables.size() + n, -1);
}

/** Append `count` names, giving the existing entries default names if needed. */
void append_names(
        std::vector<std::string>& names,
        int current_count,
        size_t count,
        const std::optional<std::vector<std::string>>& new_names,
        const char* prefix)
{
    if (!new_names && names.empty())
        return;
    for (int i = (int)names.size(); i < current_count; ++i)
        names.push_back(prefix + std::to_string(i));
    for (size_t i = 0; i < count; ++i) {
        if (new_names)
            names.push_back((*new_names)[i]);
        else
            names.push_back(prefix + std::to_string(current_count + i));
    }
}

using DoubleOrArray = std::variant<double, InputArray<double>>;
using TypeOrArray = std::variant<VariableType, InputArray<int8_t>>;
using OptionalDoubles = std::optional<InputArray<double>>;
using OptionalInts = std::optional<InputArray<int32_t>>;
using OptionalNames = std::optional<std::vector<std::string>>;

constexpr double INF = std::numeric_limits<double>::infinity();

int add_variables(
        MathOptModel& model,
        std::optional<size_t> number,
        const std::optional<DoubleOrArray>& lower_bounds_arg,
        const std::optional<DoubleOrArray>& upper_bounds_arg,
        const std::optional<TypeOrArray>& types_arg,
        const std::optional<DoubleOrArray>& objective_coefficients_arg,
        const OptionalDoubles& initial_values_arg,
        const OptionalNames& names)
{
    ScalarOrInput<double> lower_bounds = scalar_or_input(lower_bounds_arg, -INF);
    ScalarOrInput<double> upper_bounds = scalar_or_input(upper_bounds_arg, INF);
    ScalarOrInput<int8_t> types = scalar_or_input(types_arg, (int8_t)VariableType::Continuous);
    ScalarOrInput<double> objective_coefficients = scalar_or_input(objective_coefficients_arg, 0.0);
    Input<double> initial_values = input(initial_values_arg);

    std::optional<size_t> count = number;
    merge_count(count, lower_bounds, "lower_bounds");
    merge_count(count, upper_bounds, "upper_bounds");
    merge_count(count, types, "types");
    merge_count(count, objective_coefficients, "objective_coefficients");
    merge_count(count, initial_values, "initial_values");
    if (names)
        merge_count(count, names->size(), "names");
    if (!count)
        invalid("the number of variables is unknown: give 'number' or an array.");

    nb::gil_scoped_release release;
    size_t n = *count;
    int first = model.number_of_variables();
    check_capacity(first, n, "variables");
    for (size_t i = 0; i < n; ++i) {
        int8_t type = types[i];
        if (type < 0 || type > (int8_t)VariableType::Integer)
            invalid("types[" + std::to_string(i) + "] = " + std::to_string(type)
                    + " is not a VariableType value.");
    }

    model.variables_lower_bounds.reserve(first + n);
    model.variables_upper_bounds.reserve(first + n);
    model.variables_types.reserve(first + n);
    model.objective_coefficients.reserve(first + n);
    for (size_t i = 0; i < n; ++i) {
        model.variables_lower_bounds.push_back(lower_bounds[i]);
        model.variables_upper_bounds.push_back(upper_bounds[i]);
        model.variables_types.push_back(static_cast<VariableType>(types[i]));
        model.objective_coefficients.push_back(objective_coefficients[i]);
    }
    if (initial_values.given || !model.variables_initial_values.empty()) {
        model.variables_initial_values.resize(first, 0.0);
        if (initial_values.given)
            model.variables_initial_values.insert(model.variables_initial_values.end(),
                    initial_values.data, initial_values.data + n);
        else
            model.variables_initial_values.resize(first + n, 0.0);
    }
    append_names(model.variables_names, first, n, names, "x");
    return first;
}

int add_constraints(
        MathOptModel& model,
        std::optional<size_t> number,
        const std::optional<DoubleOrArray>& lower_bounds_arg,
        const std::optional<DoubleOrArray>& upper_bounds_arg,
        const OptionalInts& starts_arg,
        const OptionalInts& variables_arg,
        const OptionalDoubles& coefficients_arg,
        const OptionalInts& quadratic_starts_arg,
        const OptionalInts& quadratic_variables_1_arg,
        const OptionalInts& quadratic_variables_2_arg,
        const OptionalDoubles& quadratic_coefficients_arg,
        const OptionalInts& nonlinear_starts_arg,
        const std::optional<OperatorsInput>& nonlinear_operators_arg,
        const OptionalDoubles& nonlinear_values_arg,
        const OptionalInts& nonlinear_variables_arg,
        const OptionalInts& nonlinear_number_of_children_arg,
        const OptionalNames& names)
{
    ScalarOrInput<double> lower_bounds = scalar_or_input(lower_bounds_arg, -INF);
    ScalarOrInput<double> upper_bounds = scalar_or_input(upper_bounds_arg, INF);
    Input<int32_t> starts = input(starts_arg);
    Input<int32_t> variables = input(variables_arg);
    Input<double> coefficients = input(coefficients_arg);
    Input<int32_t> quadratic_starts = input(quadratic_starts_arg);
    Input<int32_t> quadratic_variables_1 = input(quadratic_variables_1_arg);
    Input<int32_t> quadratic_variables_2 = input(quadratic_variables_2_arg);
    Input<double> quadratic_coefficients = input(quadratic_coefficients_arg);
    Input<int32_t> nonlinear_starts = input(nonlinear_starts_arg);
    TreesInput trees = trees_input(
            nonlinear_operators_arg, nonlinear_values_arg, nonlinear_variables_arg,
            nonlinear_number_of_children_arg, "nonlinear_");

    std::optional<size_t> count = number;
    merge_count(count, lower_bounds, "lower_bounds");
    merge_count(count, upper_bounds, "upper_bounds");
    merge_count(count, starts, "starts");
    merge_count(count, quadratic_starts, "quadratic_starts");
    merge_count(count, nonlinear_starts, "nonlinear_starts");
    if (names)
        merge_count(count, names->size(), "names");
    if (!count)
        invalid("the number of constraints is unknown: give 'number' or an array.");

    // Linear terms.
    if (variables.given != coefficients.given)
        invalid("variables and coefficients must be given together.");
    check_size(coefficients, variables.size, "coefficients", "variables");

    // Quadratic terms.
    if (quadratic_variables_1.given != quadratic_variables_2.given
            || quadratic_variables_1.given != quadratic_coefficients.given)
        invalid("quadratic_variables_1, quadratic_variables_2 and quadratic_coefficients "
                "must be given together.");
    check_size(quadratic_variables_2, quadratic_variables_1.size,
            "quadratic_variables_2", "quadratic_variables_1");
    check_size(quadratic_coefficients, quadratic_variables_1.size,
            "quadratic_coefficients", "quadratic_variables_1");

    nb::gil_scoped_release release;
    size_t m = *count;
    int n = model.number_of_variables();
    int number_of_defined_variables = model.number_of_defined_variables();
    int first = model.number_of_constraints();
    check_capacity(first, m, "constraints");
    check_capacity(model.number_of_elements(), variables.size, "linear terms");
    check_capacity(model.quadratic_elements_variables_1.size(), quadratic_variables_1.size,
            "quadratic terms");
    check_capacity(model.nonlinear_elements_operators.size(), trees.size(),
            "nonlinear nodes");

    check_starts(starts, m, variables.size, "starts", "variables");
    check_variables(variables, n, "variables");
    check_starts(quadratic_starts, m, quadratic_variables_1.size,
            "quadratic_starts", "quadratic_variables_1");
    check_variables(quadratic_variables_1, n, "quadratic_variables_1");
    check_variables(quadratic_variables_2, n, "quadratic_variables_2");
    check_starts(nonlinear_starts, m, trees.size(),
            "nonlinear_starts", "nonlinear_operators");
    if (nonlinear_starts.given && trees.size() > 0) {
        for (size_t i = 0; i < m; ++i) {
            size_t begin = nonlinear_starts[i];
            size_t end = (i + 1 < m)? (size_t)nonlinear_starts[i + 1]: trees.size();
            check_tree(trees, begin, end, n, number_of_defined_variables, true,
                    "constraint " + std::to_string(i) + " of the batch");
        }
    }

    // Everything is valid: append.
    int linear_offset = model.number_of_elements();
    int quadratic_offset = (int)model.quadratic_elements_variables_1.size();
    int nonlinear_offset = (int)model.nonlinear_elements_operators.size();
    model.constraints_lower_bounds.reserve(first + m);
    model.constraints_upper_bounds.reserve(first + m);
    model.constraints_starts.reserve(first + m);
    for (size_t i = 0; i < m; ++i) {
        model.constraints_lower_bounds.push_back(lower_bounds[i]);
        model.constraints_upper_bounds.push_back(upper_bounds[i]);
        model.constraints_starts.push_back(linear_offset + (starts.given? starts[i]: 0));
    }
    model.elements_variables.insert(model.elements_variables.end(),
            variables.data, variables.data + variables.size);
    model.elements_coefficients.insert(model.elements_coefficients.end(),
            coefficients.data, coefficients.data + coefficients.size);

    model.quadratic_elements_constraints_starts.resize(first, quadratic_offset);
    model.quadratic_elements_constraints_starts.reserve(first + m);
    for (size_t i = 0; i < m; ++i)
        model.quadratic_elements_constraints_starts.push_back(
                quadratic_offset + (quadratic_starts.given? quadratic_starts[i]: 0));
    model.quadratic_elements_variables_1.insert(model.quadratic_elements_variables_1.end(),
            quadratic_variables_1.data, quadratic_variables_1.data + quadratic_variables_1.size);
    model.quadratic_elements_variables_2.insert(model.quadratic_elements_variables_2.end(),
            quadratic_variables_2.data, quadratic_variables_2.data + quadratic_variables_2.size);
    model.quadratic_elements_coefficients.insert(model.quadratic_elements_coefficients.end(),
            quadratic_coefficients.data, quadratic_coefficients.data + quadratic_coefficients.size);

    if (trees.size() > 0 || !model.nonlinear_elements_constraints_starts.empty()) {
        model.nonlinear_elements_constraints_starts.resize(first, nonlinear_offset);
        model.nonlinear_elements_constraints_starts.reserve(first + m);
        for (size_t i = 0; i < m; ++i)
            model.nonlinear_elements_constraints_starts.push_back(
                    nonlinear_offset + (nonlinear_starts.given? nonlinear_starts[i]: 0));
        append_trees(trees,
                model.nonlinear_elements_operators,
                model.nonlinear_elements_values,
                model.nonlinear_elements_variables,
                model.nonlinear_elements_number_of_children);
    }

    if (!model.constraints_functions.empty())
        model.constraints_functions.resize(first + m);
    append_names(model.constraints_names, first, m, names, "c");
    return first;
}

int add_defined_variables(
        MathOptModel& model,
        const InputArray<int32_t>& starts_arg,
        const OperatorsInput& operators_arg,
        const OptionalDoubles& values_arg,
        const OptionalInts& variables_arg,
        const OptionalInts& number_of_children_arg)
{
    Input<int32_t> starts = input(std::optional<InputArray<int32_t>>(starts_arg));
    TreesInput trees = trees_input(
            std::optional<OperatorsInput>(operators_arg), values_arg, variables_arg,
            number_of_children_arg, "");

    nb::gil_scoped_release release;
    size_t count = starts.size;
    int n = model.number_of_variables();
    int first = model.number_of_defined_variables();
    check_capacity(first, count, "defined variables");
    check_capacity(model.defined_variables_elements_operators.size(), trees.size(),
            "defined-variable nodes");
    if (count == 0 && trees.size() > 0)
        invalid("starts is empty, but operators isn't.");
    check_starts(starts, count, trees.size(), "starts", "operators");
    for (size_t i = 0; i < count; ++i) {
        size_t begin = starts[i];
        size_t end = (i + 1 < count)? (size_t)starts[i + 1]: trees.size();
        // A defined variable may only reference defined variables with a lower index.
        check_tree(trees, begin, end, n, first + (int)i, false,
                "defined variable " + std::to_string(i) + " of the batch (index "
                + std::to_string(first + i) + ", which may only reference lower indices)");
    }

    int offset = (int)model.defined_variables_elements_operators.size();
    for (size_t i = 0; i < count; ++i)
        model.defined_variables_starts.push_back(offset + starts[i]);
    append_trees(trees,
            model.defined_variables_elements_operators,
            model.defined_variables_elements_values,
            model.defined_variables_elements_variables,
            model.defined_variables_elements_number_of_children);
    return first;
}

void set_objective(
        MathOptModel& model,
        std::optional<ObjectiveDirection> direction,
        const std::optional<DoubleOrArray>& coefficients_arg,
        const OptionalInts& quadratic_variables_1_arg,
        const OptionalInts& quadratic_variables_2_arg,
        const OptionalDoubles& quadratic_coefficients_arg,
        const std::optional<OperatorsInput>& nonlinear_operators_arg,
        const OptionalDoubles& nonlinear_values_arg,
        const OptionalInts& nonlinear_variables_arg,
        const OptionalInts& nonlinear_number_of_children_arg)
{
    ScalarOrInput<double> coefficients = scalar_or_input(coefficients_arg, 0.0);
    Input<int32_t> quadratic_variables_1 = input(quadratic_variables_1_arg);
    Input<int32_t> quadratic_variables_2 = input(quadratic_variables_2_arg);
    Input<double> quadratic_coefficients = input(quadratic_coefficients_arg);
    TreesInput trees = trees_input(
            nonlinear_operators_arg, nonlinear_values_arg, nonlinear_variables_arg,
            nonlinear_number_of_children_arg, "nonlinear_");

    if (quadratic_variables_1.given != quadratic_variables_2.given
            || quadratic_variables_1.given != quadratic_coefficients.given)
        invalid("quadratic_variables_1, quadratic_variables_2 and quadratic_coefficients "
                "must be given together.");
    check_size(quadratic_variables_2, quadratic_variables_1.size,
            "quadratic_variables_2", "quadratic_variables_1");
    check_size(quadratic_coefficients, quadratic_variables_1.size,
            "quadratic_coefficients", "quadratic_variables_1");

    nb::gil_scoped_release release;
    int n = model.number_of_variables();
    if (coefficients.array.given && coefficients.array.size != (size_t)n)
        invalid("coefficients has " + std::to_string(coefficients.array.size)
                + " entries, but the model has " + std::to_string(n) + " variables.");
    check_variables(quadratic_variables_1, n, "quadratic_variables_1");
    check_variables(quadratic_variables_2, n, "quadratic_variables_2");
    if (trees.given())
        check_tree(trees, 0, trees.size(), n, model.number_of_defined_variables(), true,
                "the objective");

    if (direction)
        model.objective_direction = *direction;
    if (coefficients_arg) {
        for (int j = 0; j < n; ++j)
            model.objective_coefficients[j] = coefficients[j];
    }
    if (quadratic_variables_1.given) {
        model.objective_quadratic_elements_variables_1.assign(
                quadratic_variables_1.data, quadratic_variables_1.data + quadratic_variables_1.size);
        model.objective_quadratic_elements_variables_2.assign(
                quadratic_variables_2.data, quadratic_variables_2.data + quadratic_variables_2.size);
        model.objective_quadratic_elements_coefficients.assign(
                quadratic_coefficients.data, quadratic_coefficients.data + quadratic_coefficients.size);
    }
    if (trees.given()) {
        model.objective_nonlinear_elements_operators.clear();
        model.objective_nonlinear_elements_values.clear();
        model.objective_nonlinear_elements_variables.clear();
        model.objective_nonlinear_elements_number_of_children.clear();
        append_trees(trees,
                model.objective_nonlinear_elements_operators,
                model.objective_nonlinear_elements_values,
                model.objective_nonlinear_elements_variables,
                model.objective_nonlinear_elements_number_of_children);
    }
}

}

NB_MODULE(_core, m)
{
    m.doc() = "Python bindings for MathOptSolversCMake's model and file formats.";

    nb::enum_<ObjectiveDirection>(m, "ObjectiveDirection", nb::is_arithmetic())
        .value("Minimize", ObjectiveDirection::Minimize)
        .value("Maximize", ObjectiveDirection::Maximize);

    nb::enum_<VariableType>(m, "VariableType", nb::is_arithmetic())
        .value("Continuous", VariableType::Continuous)
        .value("Binary", VariableType::Binary)
        .value("Integer", VariableType::Integer);

    nb::enum_<ConstraintSense>(m, "ConstraintSense", nb::is_arithmetic())
        .value("LessThanOrEqualTo", ConstraintSense::LessThanOrEqualTo)
        .value("GreaterThanOrEqualTo", ConstraintSense::GreaterThanOrEqualTo)
        .value("Equality", ConstraintSense::Equality)
        .value("Range", ConstraintSense::Range)
        .value("Free", ConstraintSense::Free);

    nb::enum_<ConstraintClass>(m, "ConstraintClass", nb::is_arithmetic())
        .value("Empty", ConstraintClass::Empty)
        .value("Free", ConstraintClass::Free)
        .value("Singleton", ConstraintClass::Singleton)
        .value("Aggregation", ConstraintClass::Aggregation)
        .value("Precedence", ConstraintClass::Precedence)
        .value("VariableBound", ConstraintClass::VariableBound)
        .value("SetPartitioning", ConstraintClass::SetPartitioning)
        .value("SetPacking", ConstraintClass::SetPacking)
        .value("SetCovering", ConstraintClass::SetCovering)
        .value("Cardinality", ConstraintClass::Cardinality)
        .value("InvariantKnapsack", ConstraintClass::InvariantKnapsack)
        .value("EquationKnapsack", ConstraintClass::EquationKnapsack)
        .value("Binpacking", ConstraintClass::Binpacking)
        .value("Knapsack", ConstraintClass::Knapsack)
        .value("IntegerKnapsack", ConstraintClass::IntegerKnapsack)
        .value("MixedBinary", ConstraintClass::MixedBinary)
        .value("GeneralLinear", ConstraintClass::GeneralLinear)
        .value("GeneralQuadratic", ConstraintClass::GeneralQuadratic)
        .value("GeneralBlackBox", ConstraintClass::GeneralBlackBox)
        .value("GeneralNonlinear", ConstraintClass::GeneralNonlinear);

    nb::class_<MathOptModel>(m, "Model",
            "An optimization model. Its arrays are read-only views; see mathopt.hpp "
            "for the meaning of each field. Model() is an empty model, built with "
            "add_variables(), add_constraints(), add_defined_variables() and "
            "set_objective(); growing a model invalidates the views obtained before.")
        .def(nb::init<>())
        .def_prop_rw("objective_direction",
                [](const MathOptModel& model) { return model.objective_direction; },
                [](MathOptModel& model, ObjectiveDirection direction) {
                    model.objective_direction = direction;
                })

        // Building.
        .def("add_variables", &add_variables,
                "number"_a = nb::none(),
                "lower_bounds"_a = nb::none(),
                "upper_bounds"_a = nb::none(),
                "types"_a = nb::none(),
                "objective_coefficients"_a = nb::none(),
                "initial_values"_a = nb::none(),
                "names"_a = nb::none(),
                "Append a batch of variables and return the index of the first one.\n\n"
                "The batch size is 'number' or the length of the arrays given. "
                "lower_bounds, upper_bounds, types and objective_coefficients are each "
                "a scalar (for every variable) or an array; they default to -inf, +inf, "
                "VariableType.Continuous and 0. types arrays hold VariableType values.")
        .def("add_constraints", &add_constraints,
                "number"_a = nb::none(),
                "lower_bounds"_a = nb::none(),
                "upper_bounds"_a = nb::none(),
                "starts"_a = nb::none(),
                "variables"_a = nb::none(),
                "coefficients"_a = nb::none(),
                "quadratic_starts"_a = nb::none(),
                "quadratic_variables_1"_a = nb::none(),
                "quadratic_variables_2"_a = nb::none(),
                "quadratic_coefficients"_a = nb::none(),
                "nonlinear_starts"_a = nb::none(),
                "nonlinear_operators"_a = nb::none(),
                "nonlinear_values"_a = nb::none(),
                "nonlinear_variables"_a = nb::none(),
                "nonlinear_number_of_children"_a = nb::none(),
                "names"_a = nb::none(),
                "Append a batch of constraints, lower_bounds <= body <= upper_bounds, "
                "and return the index of the first one.\n\n"
                "Bounds are scalars or arrays (default -inf / +inf). The linear terms "
                "(variables, coefficients), the quadratic terms (quadratic_variables_1, "
                "quadratic_variables_2, quadratic_coefficients) and the expression trees "
                "(nonlinear_operators, nonlinear_values, nonlinear_variables, "
                "nonlinear_number_of_children, laid out as in mathopt.hpp) of all the "
                "batch's constraints are concatenated; starts, quadratic_starts and "
                "nonlinear_starts give the offset of each constraint's first entry in "
                "them (CSR, one entry per constraint, the first one 0). Operators are "
                "bytes or an int8/uint8 array of char codes; nonlinear_values and "
                "nonlinear_variables default to 0 and -1.")
        .def("add_defined_variables", &add_defined_variables,
                "starts"_a,
                "operators"_a,
                "values"_a = nb::none(),
                "variables"_a = nb::none(),
                "number_of_children"_a = nb::none(),
                "Append a batch of defined variables (common subexpressions), one "
                "expression tree each, and return the index of the first one. starts "
                "gives the offset of each tree in the concatenated arrays. A defined "
                "variable may only reference defined variables with a lower index.")
        .def("set_objective", &set_objective,
                "direction"_a = nb::none(),
                "coefficients"_a = nb::none(),
                "quadratic_variables_1"_a = nb::none(),
                "quadratic_variables_2"_a = nb::none(),
                "quadratic_coefficients"_a = nb::none(),
                "nonlinear_operators"_a = nb::none(),
                "nonlinear_values"_a = nb::none(),
                "nonlinear_variables"_a = nb::none(),
                "nonlinear_number_of_children"_a = nb::none(),
                "Set the objective: each part given (direction, the dense linear "
                "coefficients, the quadratic terms, the expression tree) replaces the "
                "current one; the others are kept. Pass empty arrays to remove a part.")
        .def_prop_ro("number_of_variables",
                [](const MathOptModel& model) { return model.number_of_variables(); })
        .def_prop_ro("number_of_constraints",
                [](const MathOptModel& model) { return model.number_of_constraints(); })
        .def_prop_ro("number_of_elements",
                [](const MathOptModel& model) { return model.number_of_elements(); })

        // Variables.
        .def_prop_ro("variables_lower_bounds",
                view_member<double, &MathOptModel::variables_lower_bounds>)
        .def_prop_ro("variables_upper_bounds",
                view_member<double, &MathOptModel::variables_upper_bounds>)
        .def_prop_ro("variables_types", variables_types_array, nb::rv_policy::move,
                "Variable types as an int8 array of VariableType values.")
        .def_prop_ro("variables_names",
                [](const MathOptModel& model) { return model.variables_names; },
                "Variable names (empty if the file format doesn't store them).")
        .def_prop_ro("variables_initial_values",
                view_member<double, &MathOptModel::variables_initial_values>)
        .def("variable_name", &MathOptModel::variable_name, "variable_id"_a)

        // Constraints.
        .def_prop_ro("constraints_lower_bounds",
                view_member<double, &MathOptModel::constraints_lower_bounds>)
        .def_prop_ro("constraints_upper_bounds",
                view_member<double, &MathOptModel::constraints_upper_bounds>)
        .def_prop_ro("constraints_names",
                [](const MathOptModel& model) { return model.constraints_names; },
                "Constraint names (empty if the file format doesn't store them).")
        .def("constraint_name", &MathOptModel::constraint_name, "constraint_id"_a)
        .def("constraint_sense", &MathOptModel::constraint_sense, "constraint_id"_a)
        .def("constraint_class", &MathOptModel::constraint_class, "constraint_id"_a)
        .def_prop_ro("constraints_classes", constraints_classes_array, nb::rv_policy::move,
                "Constraint classes as an int8 array of ConstraintClass values.")
        .def("constraint_end", &MathOptModel::constraint_end, "constraint_id"_a)
        .def("quadratic_constraint_end", &MathOptModel::quadratic_constraint_end, "constraint_id"_a)
        .def("nonlinear_constraint_end", &MathOptModel::nonlinear_constraint_end, "constraint_id"_a)
        .def("format_constraint",
                [](const MathOptModel& model, int constraint_id) {
                    std::ostringstream os;
                    model.format_constraint(os, constraint_id);
                    return os.str();
                },
                "constraint_id"_a)

        // Linear structures.
        .def_prop_ro("objective_coefficients",
                view_member<double, &MathOptModel::objective_coefficients>)
        .def_prop_ro("constraints_starts",
                view_member<int, &MathOptModel::constraints_starts>)
        .def_prop_ro("elements_variables",
                view_member<int, &MathOptModel::elements_variables>)
        .def_prop_ro("elements_coefficients",
                view_member<double, &MathOptModel::elements_coefficients>)

        // Quadratic structures.
        .def_prop_ro("objective_quadratic_elements_variables_1",
                view_member<int, &MathOptModel::objective_quadratic_elements_variables_1>)
        .def_prop_ro("objective_quadratic_elements_variables_2",
                view_member<int, &MathOptModel::objective_quadratic_elements_variables_2>)
        .def_prop_ro("objective_quadratic_elements_coefficients",
                view_member<double, &MathOptModel::objective_quadratic_elements_coefficients>)
        .def_prop_ro("quadratic_elements_constraints_starts",
                view_member<int, &MathOptModel::quadratic_elements_constraints_starts>)
        .def_prop_ro("quadratic_elements_variables_1",
                view_member<int, &MathOptModel::quadratic_elements_variables_1>)
        .def_prop_ro("quadratic_elements_variables_2",
                view_member<int, &MathOptModel::quadratic_elements_variables_2>)
        .def_prop_ro("quadratic_elements_coefficients",
                view_member<double, &MathOptModel::quadratic_elements_coefficients>)

        // Nonlinear structures.
        .def_prop_ro("objective_nonlinear_elements_operators",
                [](const MathOptModel& model) { return to_bytes(model.objective_nonlinear_elements_operators); },
                "Operator char codes of the objective's expression tree (see mathopt.hpp), as bytes.")
        .def_prop_ro("objective_nonlinear_elements_values",
                view_member<double, &MathOptModel::objective_nonlinear_elements_values>)
        .def_prop_ro("objective_nonlinear_elements_variables",
                view_member<int, &MathOptModel::objective_nonlinear_elements_variables>)
        .def_prop_ro("objective_nonlinear_elements_number_of_children",
                view_member<int, &MathOptModel::objective_nonlinear_elements_number_of_children>)
        .def_prop_ro("nonlinear_elements_constraints_starts",
                view_member<int, &MathOptModel::nonlinear_elements_constraints_starts>)
        .def_prop_ro("nonlinear_elements_operators",
                [](const MathOptModel& model) { return to_bytes(model.nonlinear_elements_operators); },
                "Operator char codes of the constraints' expression trees (see mathopt.hpp), as bytes.")
        .def_prop_ro("nonlinear_elements_values",
                view_member<double, &MathOptModel::nonlinear_elements_values>)
        .def_prop_ro("nonlinear_elements_variables",
                view_member<int, &MathOptModel::nonlinear_elements_variables>)
        .def_prop_ro("nonlinear_elements_number_of_children",
                view_member<int, &MathOptModel::nonlinear_elements_number_of_children>)

        // Defined variables (common subexpressions).
        .def_prop_ro("number_of_defined_variables",
                [](const MathOptModel& model) { return model.number_of_defined_variables(); })
        .def("defined_variable_end", &MathOptModel::defined_variable_end, "defined_variable_id"_a)
        .def_prop_ro("defined_variables_starts",
                view_member<int, &MathOptModel::defined_variables_starts>)
        .def_prop_ro("defined_variables_elements_operators",
                [](const MathOptModel& model) { return to_bytes(model.defined_variables_elements_operators); },
                "Operator char codes of the defined variables' expression trees (see mathopt.hpp), as bytes.")
        .def_prop_ro("defined_variables_elements_values",
                view_member<double, &MathOptModel::defined_variables_elements_values>)
        .def_prop_ro("defined_variables_elements_variables",
                view_member<int, &MathOptModel::defined_variables_elements_variables>)
        .def_prop_ro("defined_variables_elements_number_of_children",
                view_member<int, &MathOptModel::defined_variables_elements_number_of_children>)

        // Model kind.
        .def("has_non_continuous_variables", &MathOptModel::has_non_continuous_variables)
        .def("has_defined_variables", &MathOptModel::has_defined_variables)
        .def("has_quadratic", &MathOptModel::has_quadratic)
        .def("has_nonlinear", &MathOptModel::has_nonlinear)
        .def("is_lp", &MathOptModel::is_lp)
        .def("is_milp", &MathOptModel::is_milp)

        // Evaluations and checks.
        .def("evaluate_objective", &MathOptModel::evaluate_objective, "solution"_a)
        .def("evaluate_constraint", &MathOptModel::evaluate_constraint, "solution"_a, "constraint_id"_a)
        .def("check", &MathOptModel::check, "verbosity_level"_a = 0,
                "Check that the model's arrays are consistent.")
        .def("check_solution", &MathOptModel::check_solution, "solution"_a, "verbosity_level"_a = 0)
        .def("format",
                [](const MathOptModel& model, int verbosity_level) {
                    std::ostringstream os;
                    model.format(os, verbosity_level);
                    return os.str();
                },
                "verbosity_level"_a = 1)
        .def("__repr__",
                [](const MathOptModel& model) {
                    return "<mathoptmodel.Model: "
                        + std::to_string(model.number_of_variables()) + " variables, "
                        + std::to_string(model.number_of_constraints()) + " constraints>";
                });

    // Readers, from the raw file contents.
    m.def("read_mps_bytes",
            [](nb::bytes data) { return read_from_bytes(read_mps, data); },
            "data"_a, "Read a model from the contents of a free-format MPS file.");
    m.def("read_lp_bytes",
            [](nb::bytes data) { return read_from_bytes(read_lp, data); },
            "data"_a, "Read a model from the contents of a CPLEX-style .lp file.");
    m.def("read_json_bytes",
            [](nb::bytes data) { return read_from_bytes(read_json, data); },
            "data"_a, "Read a model from the contents of a MathOptSolversCMake JSON file.");
    m.def("read_nl_bytes",
            [](nb::bytes data) { return read_from_bytes(read_nl, data); },
            "data"_a, "Read a model from the contents of an AMPL .nl file.");

    // Writers, to a str.
    m.def("write_mps_string",
            [](const MathOptModel& model) { return write_to_string(write_mps, model); },
            "model"_a, "Write a model as a free-format MPS file, returned as a str.");
    m.def("write_lp_string",
            [](const MathOptModel& model) { return write_to_string(write_lp, model); },
            "model"_a, "Write a model as a CPLEX-style .lp file, returned as a str.");
    m.def("write_json_string",
            [](const MathOptModel& model) { return write_to_string(write_json, model); },
            "model"_a, "Write a model as a MathOptSolversCMake JSON file, returned as a str.");
    m.def("write_nl_string",
            [](const MathOptModel& model) { return write_to_string(write_nl, model); },
            "model"_a, "Write a model as an AMPL .nl file (text format), returned as a str.");

    m.def("dual", &dual, "model"_a, "Build the dual of a linear program.");
    m.def("expand_defined_variables", &expand_defined_variables, "model"_a,
            "Return an equivalent model whose expression trees reference no defined variable.");
}

#include "mathoptsolverscmake/mathoptmodel_c.h"

#include "mathoptsolverscmake/mathopt.hpp"
#include "mathoptsolverscmake/mathopt_build.hpp"
#include "mathoptsolverscmake/mathopt_mps.hpp"
#include "mathoptsolverscmake/mathopt_lp.hpp"
#include "mathoptsolverscmake/mathopt_json.hpp"
#include "mathoptsolverscmake/mathopt_nl.hpp"
#include "mathoptsolverscmake/dual.hpp"

#include <cstdlib>
#include <cstring>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace mathoptsolverscmake;

struct mom_model
{
    MathOptModel model;
};

namespace
{

static_assert(sizeof(int) == sizeof(int32_t), "MathOptModel indexes are int32.");
static_assert(sizeof(char) == sizeof(int8_t), "Operators are int8.");

void set_error(
        char* error_message,
        size_t error_message_size,
        const char* message)
{
    if (error_message == nullptr || error_message_size == 0)
        return;
    std::strncpy(error_message, message, error_message_size - 1);
    error_message[error_message_size - 1] = '\0';
}

/** Run 'function', translating exceptions into a status and an error message. */
template <typename Function>
int guarded(
        char* error_message,
        size_t error_message_size,
        Function function)
{
    try {
        function();
        return MOM_OK;
    } catch (const std::invalid_argument& e) {
        set_error(error_message, error_message_size, e.what());
        return MOM_INVALID_ARGUMENT;
    } catch (const std::bad_alloc&) {
        set_error(error_message, error_message_size, "out of memory");
        return MOM_ERROR;
    } catch (const std::exception& e) {
        set_error(error_message, error_message_size, e.what());
        return MOM_ERROR;
    } catch (...) {
        set_error(error_message, error_message_size, "unknown error");
        return MOM_ERROR;
    }
}

void check_model(const mom_model* model)
{
    if (model == nullptr)
        throw std::invalid_argument("model is null.");
}

void check_count(int64_t count, const char* name)
{
    if (count < 0)
        throw std::invalid_argument(std::string(name) + " is negative.");
}

/** Copy 'contents' into a buffer allocated with malloc. */
void output_string(
        const std::string& contents,
        char** data,
        int64_t* size)
{
    if (data == nullptr || size == nullptr)
        throw std::invalid_argument("data and size must not be null.");
    char* buffer = static_cast<char*>(std::malloc(contents.size() + 1));
    if (buffer == nullptr)
        throw std::bad_alloc();
    std::memcpy(buffer, contents.data(), contents.size());
    buffer[contents.size()] = '\0';
    *data = buffer;
    *size = (int64_t)contents.size();
}

/** Build the names of a batch from concatenated bytes and offsets. */
bool input_names(
        int64_t number,
        const char* names_data,
        const int64_t* names_offsets,
        std::vector<std::string>& names)
{
    if (names_offsets == nullptr)
        return false;
    if (names_data == nullptr && names_offsets[number] > 0)
        throw std::invalid_argument("names_data is null.");
    names.reserve(number);
    for (int64_t i = 0; i < number; ++i) {
        int64_t begin = names_offsets[i];
        int64_t end = names_offsets[i + 1];
        if (begin < 0 || end < begin)
            throw std::invalid_argument("names_offsets must be non-negative and non-decreasing.");
        names.emplace_back(names_data + begin, names_data + end);
    }
    return true;
}

ExpressionTreesBatch trees_batch(
        int64_t number_of_nodes,
        const int8_t* operators,
        const double* values,
        const int32_t* variables,
        const int32_t* number_of_children)
{
    check_count(number_of_nodes, "number_of_nodes");
    ExpressionTreesBatch res;
    res.number_of_nodes = (size_t)number_of_nodes;
    res.operators = reinterpret_cast<const char*>(operators);
    res.values = values;
    res.variables = variables;
    res.number_of_children = number_of_children;
    return res;
}

template <typename T>
void copy_vector(const std::vector<T>& vector, void* out)
{
    if (!vector.empty())
        std::memcpy(out, vector.data(), vector.size() * sizeof(T));
}

}

/*
 * Life cycle.
 */

mom_model* mom_model_create(void)
{
    return new (std::nothrow) mom_model();
}

void mom_model_free(mom_model* model)
{
    delete model;
}

mom_model* mom_model_copy(const mom_model* model)
{
    if (model == nullptr)
        return nullptr;
    try {
        return new mom_model(*model);
    } catch (...) {
        return nullptr;
    }
}

void mom_free(void* data)
{
    std::free(data);
}

/*
 * Reading and writing.
 */

int mom_read(
        const char* data,
        int64_t size,
        int format,
        mom_model** model,
        char* error_message,
        size_t error_message_size)
{
    return guarded(error_message, error_message_size, [&]() {
        if (model == nullptr)
            throw std::invalid_argument("model must not be null.");
        check_count(size, "size");
        if (data == nullptr && size > 0)
            throw std::invalid_argument("data is null.");
        std::string contents(data, data + size);
        std::istringstream in(contents, std::ios::in | std::ios::binary);
        mom_model* res = new mom_model();
        try {
            switch (format) {
            case MOM_FORMAT_MPS: res->model = read_mps(in); break;
            case MOM_FORMAT_LP: res->model = read_lp(in); break;
            case MOM_FORMAT_JSON: res->model = read_json(in); break;
            case MOM_FORMAT_NL: res->model = read_nl(in); break;
            default: throw std::invalid_argument("unknown format " + std::to_string(format) + ".");
            }
        } catch (...) {
            delete res;
            throw;
        }
        *model = res;
    });
}

int mom_write(
        const mom_model* model,
        int format,
        char** data,
        int64_t* size,
        char* error_message,
        size_t error_message_size)
{
    return guarded(error_message, error_message_size, [&]() {
        check_model(model);
        std::ostringstream out;
        switch (format) {
        case MOM_FORMAT_MPS: write_mps(model->model, out); break;
        case MOM_FORMAT_LP: write_lp(model->model, out); break;
        case MOM_FORMAT_JSON: write_json(model->model, out); break;
        case MOM_FORMAT_NL: write_nl(model->model, out); break;
        default: throw std::invalid_argument("unknown format " + std::to_string(format) + ".");
        }
        output_string(out.str(), data, size);
    });
}

int mom_read_nl_statistics(
        const char* data,
        int64_t size,
        int64_t* values,
        char* error_message,
        size_t error_message_size)
{
    return guarded(error_message, error_message_size, [&]() {
        check_count(size, "size");
        if (data == nullptr && size > 0)
            throw std::invalid_argument("data is null.");
        if (values == nullptr)
            throw std::invalid_argument("values must not be null.");
        std::istringstream in(std::string(data, data + size), std::ios::in | std::ios::binary);
        NlStatistics statistics = read_nl_statistics(in);
        values[MOM_NL_BINARY] = statistics.binary;
        values[MOM_NL_NUMBER_OF_VARIABLES] = statistics.number_of_variables;
        values[MOM_NL_NUMBER_OF_CONSTRAINTS] = statistics.number_of_constraints;
        values[MOM_NL_NUMBER_OF_OBJECTIVES] = statistics.number_of_objectives;
        values[MOM_NL_NUMBER_OF_NONLINEAR_CONSTRAINTS] = statistics.number_of_nonlinear_constraints;
        values[MOM_NL_NUMBER_OF_NONLINEAR_OBJECTIVES] = statistics.number_of_nonlinear_objectives;
        values[MOM_NL_NUMBER_OF_NONLINEAR_VARIABLES] = statistics.number_of_nonlinear_variables;
        values[MOM_NL_NUMBER_OF_BINARY_VARIABLES] = statistics.number_of_binary_variables;
        values[MOM_NL_NUMBER_OF_INTEGER_VARIABLES] = statistics.number_of_integer_variables;
        values[MOM_NL_NUMBER_OF_ELEMENTS] = statistics.number_of_elements;
        values[MOM_NL_NUMBER_OF_OBJECTIVE_ELEMENTS] = statistics.number_of_objective_elements;
        values[MOM_NL_NUMBER_OF_DEFINED_VARIABLES] = statistics.number_of_defined_variables;
        values[MOM_NL_NUMBER_OF_LOGICAL_CONSTRAINTS] = statistics.number_of_logical_constraints;
        values[MOM_NL_NUMBER_OF_COMPLEMENTARITY_CONSTRAINTS] = statistics.number_of_complementarity_constraints;
    });
}

int mom_write_file(
        const mom_model* model,
        int format,
        const char* path,
        char* error_message,
        size_t error_message_size)
{
    return guarded(error_message, error_message_size, [&]() {
        check_model(model);
        if (path == nullptr)
            throw std::invalid_argument("path must not be null.");
        std::string file_path(path);
        switch (format) {
        case MOM_FORMAT_MPS: write_mps(model->model, file_path); break;
        case MOM_FORMAT_LP: write_lp(model->model, file_path); break;
        case MOM_FORMAT_JSON: write_json(model->model, file_path); break;
        case MOM_FORMAT_NL: write_nl(model->model, file_path); break;
        default: throw std::invalid_argument("unknown format " + std::to_string(format) + ".");
        }
    });
}

int mom_format(
        const mom_model* model,
        int verbosity_level,
        char** data,
        int64_t* size,
        char* error_message,
        size_t error_message_size)
{
    return guarded(error_message, error_message_size, [&]() {
        check_model(model);
        std::ostringstream out;
        model->model.format(out, verbosity_level);
        output_string(out.str(), data, size);
    });
}

int mom_format_constraint(
        const mom_model* model,
        int32_t constraint_id,
        char** data,
        int64_t* size,
        char* error_message,
        size_t error_message_size)
{
    return guarded(error_message, error_message_size, [&]() {
        check_model(model);
        if (constraint_id < 0 || constraint_id >= model->model.number_of_constraints())
            throw std::invalid_argument("invalid constraint index " + std::to_string(constraint_id) + ".");
        std::ostringstream out;
        model->model.format_constraint(out, constraint_id);
        output_string(out.str(), data, size);
    });
}

/*
 * Inspection.
 */

int32_t mom_number_of_variables(const mom_model* model)
{
    return (model != nullptr)? model->model.number_of_variables(): 0;
}

int32_t mom_number_of_constraints(const mom_model* model)
{
    return (model != nullptr)? model->model.number_of_constraints(): 0;
}

int32_t mom_number_of_elements(const mom_model* model)
{
    return (model != nullptr)? model->model.number_of_elements(): 0;
}

int32_t mom_number_of_defined_variables(const mom_model* model)
{
    return (model != nullptr)? model->model.number_of_defined_variables(): 0;
}

int32_t mom_objective_direction(const mom_model* model)
{
    return (model != nullptr)? (int32_t)model->model.objective_direction: 0;
}

int64_t mom_array_size(
        const mom_model* model,
        int array)
{
    if (model == nullptr)
        return -1;
    const MathOptModel& m = model->model;
    switch (array) {
    case MOM_VARIABLES_LOWER_BOUNDS: return m.variables_lower_bounds.size();
    case MOM_VARIABLES_UPPER_BOUNDS: return m.variables_upper_bounds.size();
    case MOM_VARIABLES_TYPES: return m.variables_types.size();
    case MOM_VARIABLES_INITIAL_VALUES: return m.variables_initial_values.size();
    case MOM_CONSTRAINTS_LOWER_BOUNDS: return m.constraints_lower_bounds.size();
    case MOM_CONSTRAINTS_UPPER_BOUNDS: return m.constraints_upper_bounds.size();
    case MOM_CONSTRAINTS_CLASSES: return m.number_of_constraints();
    case MOM_OBJECTIVE_COEFFICIENTS: return m.objective_coefficients.size();
    case MOM_CONSTRAINTS_STARTS: return m.constraints_starts.size();
    case MOM_ELEMENTS_VARIABLES: return m.elements_variables.size();
    case MOM_ELEMENTS_COEFFICIENTS: return m.elements_coefficients.size();
    case MOM_OBJECTIVE_QUADRATIC_ELEMENTS_VARIABLES_1: return m.objective_quadratic_elements_variables_1.size();
    case MOM_OBJECTIVE_QUADRATIC_ELEMENTS_VARIABLES_2: return m.objective_quadratic_elements_variables_2.size();
    case MOM_OBJECTIVE_QUADRATIC_ELEMENTS_COEFFICIENTS: return m.objective_quadratic_elements_coefficients.size();
    case MOM_QUADRATIC_ELEMENTS_CONSTRAINTS_STARTS: return m.quadratic_elements_constraints_starts.size();
    case MOM_QUADRATIC_ELEMENTS_VARIABLES_1: return m.quadratic_elements_variables_1.size();
    case MOM_QUADRATIC_ELEMENTS_VARIABLES_2: return m.quadratic_elements_variables_2.size();
    case MOM_QUADRATIC_ELEMENTS_COEFFICIENTS: return m.quadratic_elements_coefficients.size();
    case MOM_OBJECTIVE_NONLINEAR_ELEMENTS_OPERATORS: return m.objective_nonlinear_elements_operators.size();
    case MOM_OBJECTIVE_NONLINEAR_ELEMENTS_VALUES: return m.objective_nonlinear_elements_values.size();
    case MOM_OBJECTIVE_NONLINEAR_ELEMENTS_VARIABLES: return m.objective_nonlinear_elements_variables.size();
    case MOM_OBJECTIVE_NONLINEAR_ELEMENTS_NUMBER_OF_CHILDREN: return m.objective_nonlinear_elements_number_of_children.size();
    case MOM_NONLINEAR_ELEMENTS_CONSTRAINTS_STARTS: return m.nonlinear_elements_constraints_starts.size();
    case MOM_NONLINEAR_ELEMENTS_OPERATORS: return m.nonlinear_elements_operators.size();
    case MOM_NONLINEAR_ELEMENTS_VALUES: return m.nonlinear_elements_values.size();
    case MOM_NONLINEAR_ELEMENTS_VARIABLES: return m.nonlinear_elements_variables.size();
    case MOM_NONLINEAR_ELEMENTS_NUMBER_OF_CHILDREN: return m.nonlinear_elements_number_of_children.size();
    case MOM_DEFINED_VARIABLES_STARTS: return m.defined_variables_starts.size();
    case MOM_DEFINED_VARIABLES_ELEMENTS_OPERATORS: return m.defined_variables_elements_operators.size();
    case MOM_DEFINED_VARIABLES_ELEMENTS_VALUES: return m.defined_variables_elements_values.size();
    case MOM_DEFINED_VARIABLES_ELEMENTS_VARIABLES: return m.defined_variables_elements_variables.size();
    case MOM_DEFINED_VARIABLES_ELEMENTS_NUMBER_OF_CHILDREN: return m.defined_variables_elements_number_of_children.size();
    default: return -1;
    }
}

int mom_copy_array(
        const mom_model* model,
        int array,
        void* out,
        char* error_message,
        size_t error_message_size)
{
    return guarded(error_message, error_message_size, [&]() {
        check_model(model);
        if (out == nullptr && mom_array_size(model, array) > 0)
            throw std::invalid_argument("out is null.");
        const MathOptModel& m = model->model;
        switch (array) {
        case MOM_VARIABLES_LOWER_BOUNDS: copy_vector(m.variables_lower_bounds, out); break;
        case MOM_VARIABLES_UPPER_BOUNDS: copy_vector(m.variables_upper_bounds, out); break;
        case MOM_VARIABLES_TYPES: {
            int8_t* types = static_cast<int8_t*>(out);
            for (size_t i = 0; i < m.variables_types.size(); ++i)
                types[i] = (int8_t)m.variables_types[i];
            break;
        }
        case MOM_VARIABLES_INITIAL_VALUES: copy_vector(m.variables_initial_values, out); break;
        case MOM_CONSTRAINTS_LOWER_BOUNDS: copy_vector(m.constraints_lower_bounds, out); break;
        case MOM_CONSTRAINTS_UPPER_BOUNDS: copy_vector(m.constraints_upper_bounds, out); break;
        case MOM_CONSTRAINTS_CLASSES: {
            int8_t* classes = static_cast<int8_t*>(out);
            for (int i = 0; i < m.number_of_constraints(); ++i)
                classes[i] = (int8_t)m.constraint_class(i);
            break;
        }
        case MOM_OBJECTIVE_COEFFICIENTS: copy_vector(m.objective_coefficients, out); break;
        case MOM_CONSTRAINTS_STARTS: copy_vector(m.constraints_starts, out); break;
        case MOM_ELEMENTS_VARIABLES: copy_vector(m.elements_variables, out); break;
        case MOM_ELEMENTS_COEFFICIENTS: copy_vector(m.elements_coefficients, out); break;
        case MOM_OBJECTIVE_QUADRATIC_ELEMENTS_VARIABLES_1: copy_vector(m.objective_quadratic_elements_variables_1, out); break;
        case MOM_OBJECTIVE_QUADRATIC_ELEMENTS_VARIABLES_2: copy_vector(m.objective_quadratic_elements_variables_2, out); break;
        case MOM_OBJECTIVE_QUADRATIC_ELEMENTS_COEFFICIENTS: copy_vector(m.objective_quadratic_elements_coefficients, out); break;
        case MOM_QUADRATIC_ELEMENTS_CONSTRAINTS_STARTS: copy_vector(m.quadratic_elements_constraints_starts, out); break;
        case MOM_QUADRATIC_ELEMENTS_VARIABLES_1: copy_vector(m.quadratic_elements_variables_1, out); break;
        case MOM_QUADRATIC_ELEMENTS_VARIABLES_2: copy_vector(m.quadratic_elements_variables_2, out); break;
        case MOM_QUADRATIC_ELEMENTS_COEFFICIENTS: copy_vector(m.quadratic_elements_coefficients, out); break;
        case MOM_OBJECTIVE_NONLINEAR_ELEMENTS_OPERATORS: copy_vector(m.objective_nonlinear_elements_operators, out); break;
        case MOM_OBJECTIVE_NONLINEAR_ELEMENTS_VALUES: copy_vector(m.objective_nonlinear_elements_values, out); break;
        case MOM_OBJECTIVE_NONLINEAR_ELEMENTS_VARIABLES: copy_vector(m.objective_nonlinear_elements_variables, out); break;
        case MOM_OBJECTIVE_NONLINEAR_ELEMENTS_NUMBER_OF_CHILDREN: copy_vector(m.objective_nonlinear_elements_number_of_children, out); break;
        case MOM_NONLINEAR_ELEMENTS_CONSTRAINTS_STARTS: copy_vector(m.nonlinear_elements_constraints_starts, out); break;
        case MOM_NONLINEAR_ELEMENTS_OPERATORS: copy_vector(m.nonlinear_elements_operators, out); break;
        case MOM_NONLINEAR_ELEMENTS_VALUES: copy_vector(m.nonlinear_elements_values, out); break;
        case MOM_NONLINEAR_ELEMENTS_VARIABLES: copy_vector(m.nonlinear_elements_variables, out); break;
        case MOM_NONLINEAR_ELEMENTS_NUMBER_OF_CHILDREN: copy_vector(m.nonlinear_elements_number_of_children, out); break;
        case MOM_DEFINED_VARIABLES_STARTS: copy_vector(m.defined_variables_starts, out); break;
        case MOM_DEFINED_VARIABLES_ELEMENTS_OPERATORS: copy_vector(m.defined_variables_elements_operators, out); break;
        case MOM_DEFINED_VARIABLES_ELEMENTS_VALUES: copy_vector(m.defined_variables_elements_values, out); break;
        case MOM_DEFINED_VARIABLES_ELEMENTS_VARIABLES: copy_vector(m.defined_variables_elements_variables, out); break;
        case MOM_DEFINED_VARIABLES_ELEMENTS_NUMBER_OF_CHILDREN: copy_vector(m.defined_variables_elements_number_of_children, out); break;
        default: throw std::invalid_argument("unknown array " + std::to_string(array) + ".");
        }
    });
}

namespace
{

const std::vector<std::string>* names_list(const mom_model* model, int names)
{
    switch (names) {
    case MOM_VARIABLES_NAMES: return &model->model.variables_names;
    case MOM_CONSTRAINTS_NAMES: return &model->model.constraints_names;
    default: return nullptr;
    }
}

}

int64_t mom_names_size(
        const mom_model* model,
        int names,
        int64_t* total_size)
{
    int64_t total = 0;
    int64_t count = 0;
    if (model != nullptr) {
        const std::vector<std::string>* list = names_list(model, names);
        if (list != nullptr) {
            count = list->size();
            for (const std::string& name: *list)
                total += name.size();
        }
    }
    if (total_size != nullptr)
        *total_size = total;
    return count;
}

int mom_copy_names(
        const mom_model* model,
        int names,
        char* data,
        int64_t* offsets,
        char* error_message,
        size_t error_message_size)
{
    return guarded(error_message, error_message_size, [&]() {
        check_model(model);
        const std::vector<std::string>* list = names_list(model, names);
        if (list == nullptr)
            throw std::invalid_argument("unknown names list " + std::to_string(names) + ".");
        if (offsets == nullptr)
            throw std::invalid_argument("offsets is null.");
        int64_t offset = 0;
        for (size_t i = 0; i < list->size(); ++i) {
            offsets[i] = offset;
            const std::string& name = (*list)[i];
            if (!name.empty()) {
                if (data == nullptr)
                    throw std::invalid_argument("data is null.");
                std::memcpy(data + offset, name.data(), name.size());
            }
            offset += name.size();
        }
        offsets[list->size()] = offset;
    });
}

int32_t mom_constraint_sense(const mom_model* model, int32_t constraint_id)
{
    if (model == nullptr || constraint_id < 0 || constraint_id >= model->model.number_of_constraints())
        return -1;
    return (int32_t)model->model.constraint_sense(constraint_id);
}

int32_t mom_constraint_class(const mom_model* model, int32_t constraint_id)
{
    if (model == nullptr || constraint_id < 0 || constraint_id >= model->model.number_of_constraints())
        return -1;
    try {
        return (int32_t)model->model.constraint_class(constraint_id);
    } catch (...) {
        return -1;
    }
}

int32_t mom_check(const mom_model* model, int verbosity_level)
{
    if (model == nullptr)
        return 0;
    try {
        return model->model.check(verbosity_level)? 1: 0;
    } catch (...) {
        return 0;
    }
}

int mom_evaluate_objective(
        const mom_model* model,
        const double* solution,
        double* value,
        char* error_message,
        size_t error_message_size)
{
    return guarded(error_message, error_message_size, [&]() {
        check_model(model);
        if (solution == nullptr || value == nullptr)
            throw std::invalid_argument("solution and value must not be null.");
        std::vector<double> x(solution, solution + model->model.number_of_variables());
        *value = model->model.evaluate_objective(x);
    });
}

int mom_evaluate_constraint(
        const mom_model* model,
        const double* solution,
        int32_t constraint_id,
        double* value,
        char* error_message,
        size_t error_message_size)
{
    return guarded(error_message, error_message_size, [&]() {
        check_model(model);
        if (solution == nullptr || value == nullptr)
            throw std::invalid_argument("solution and value must not be null.");
        if (constraint_id < 0 || constraint_id >= model->model.number_of_constraints())
            throw std::invalid_argument("invalid constraint index " + std::to_string(constraint_id) + ".");
        std::vector<double> x(solution, solution + model->model.number_of_variables());
        *value = model->model.evaluate_constraint(x, constraint_id);
    });
}

int mom_check_solution(
        const mom_model* model,
        const double* solution,
        int verbosity_level,
        int32_t* feasible,
        char* error_message,
        size_t error_message_size)
{
    return guarded(error_message, error_message_size, [&]() {
        check_model(model);
        if (solution == nullptr || feasible == nullptr)
            throw std::invalid_argument("solution and feasible must not be null.");
        std::vector<double> x(solution, solution + model->model.number_of_variables());
        *feasible = model->model.check_solution(x, verbosity_level)? 1: 0;
    });
}

/*
 * Building.
 */

int mom_add_variables(
        mom_model* model,
        int64_t number,
        const double* lower_bounds, double lower_bound,
        const double* upper_bounds, double upper_bound,
        const int8_t* types, int8_t type,
        const double* objective_coefficients, double objective_coefficient,
        const double* initial_values,
        const char* names_data, const int64_t* names_offsets,
        int32_t* first,
        char* error_message,
        size_t error_message_size)
{
    return guarded(error_message, error_message_size, [&]() {
        check_model(model);
        check_count(number, "number");
        if (types == nullptr && (type < 0 || type > (int8_t)VariableType::SemiInteger))
            throw std::invalid_argument("type = " + std::to_string(type) + " is not a VariableType value.");
        std::vector<std::string> names;
        bool has_names = input_names(number, names_data, names_offsets, names);

        VariablesBatch batch;
        batch.number = (size_t)number;
        batch.lower_bounds = lower_bounds;
        batch.lower_bound = lower_bound;
        batch.upper_bounds = upper_bounds;
        batch.upper_bound = upper_bound;
        batch.types = reinterpret_cast<const signed char*>(types);
        batch.type = static_cast<VariableType>(type);
        batch.objective_coefficients = objective_coefficients;
        batch.objective_coefficient = objective_coefficient;
        batch.initial_values = initial_values;
        batch.names = has_names? &names: nullptr;
        int32_t res = add_variables(model->model, batch);
        if (first != nullptr)
            *first = res;
    });
}

int mom_add_constraints(
        mom_model* model,
        int64_t number,
        const double* lower_bounds, double lower_bound,
        const double* upper_bounds, double upper_bound,
        const int32_t* starts,
        int64_t number_of_elements,
        const int32_t* variables,
        const double* coefficients,
        const int32_t* quadratic_starts,
        int64_t number_of_quadratic_elements,
        const int32_t* quadratic_variables_1,
        const int32_t* quadratic_variables_2,
        const double* quadratic_coefficients,
        const int32_t* nonlinear_starts,
        int64_t number_of_nonlinear_nodes,
        const int8_t* nonlinear_operators,
        const double* nonlinear_values,
        const int32_t* nonlinear_variables,
        const int32_t* nonlinear_number_of_children,
        const char* names_data, const int64_t* names_offsets,
        int32_t* first,
        char* error_message,
        size_t error_message_size)
{
    return guarded(error_message, error_message_size, [&]() {
        check_model(model);
        check_count(number, "number");
        check_count(number_of_elements, "number_of_elements");
        check_count(number_of_quadratic_elements, "number_of_quadratic_elements");
        std::vector<std::string> names;
        bool has_names = input_names(number, names_data, names_offsets, names);

        ConstraintsBatch batch;
        batch.number = (size_t)number;
        batch.lower_bounds = lower_bounds;
        batch.lower_bound = lower_bound;
        batch.upper_bounds = upper_bounds;
        batch.upper_bound = upper_bound;
        batch.starts = starts;
        batch.number_of_elements = (size_t)number_of_elements;
        batch.variables = variables;
        batch.coefficients = coefficients;
        batch.quadratic_starts = quadratic_starts;
        batch.number_of_quadratic_elements = (size_t)number_of_quadratic_elements;
        batch.quadratic_variables_1 = quadratic_variables_1;
        batch.quadratic_variables_2 = quadratic_variables_2;
        batch.quadratic_coefficients = quadratic_coefficients;
        batch.nonlinear_starts = nonlinear_starts;
        batch.nonlinear = trees_batch(number_of_nonlinear_nodes, nonlinear_operators,
                nonlinear_values, nonlinear_variables, nonlinear_number_of_children);
        batch.names = has_names? &names: nullptr;
        int32_t res = add_constraints(model->model, batch);
        if (first != nullptr)
            *first = res;
    });
}

int mom_add_defined_variables(
        mom_model* model,
        int64_t number,
        const int32_t* starts,
        int64_t number_of_nodes,
        const int8_t* operators,
        const double* values,
        const int32_t* variables,
        const int32_t* number_of_children,
        int32_t* first,
        char* error_message,
        size_t error_message_size)
{
    return guarded(error_message, error_message_size, [&]() {
        check_model(model);
        check_count(number, "number");
        DefinedVariablesBatch batch;
        batch.number = (size_t)number;
        batch.starts = starts;
        batch.trees = trees_batch(number_of_nodes, operators, values, variables, number_of_children);
        int32_t res = add_defined_variables(model->model, batch);
        if (first != nullptr)
            *first = res;
    });
}

int mom_set_objective(
        mom_model* model,
        int32_t set_direction, int32_t direction,
        int32_t set_coefficients, const double* coefficients, double coefficient,
        int32_t set_quadratic,
        int64_t number_of_quadratic_elements,
        const int32_t* quadratic_variables_1,
        const int32_t* quadratic_variables_2,
        const double* quadratic_coefficients,
        int32_t set_nonlinear,
        int64_t number_of_nonlinear_nodes,
        const int8_t* nonlinear_operators,
        const double* nonlinear_values,
        const int32_t* nonlinear_variables,
        const int32_t* nonlinear_number_of_children,
        char* error_message,
        size_t error_message_size)
{
    return guarded(error_message, error_message_size, [&]() {
        check_model(model);
        if (set_direction && direction != (int32_t)ObjectiveDirection::Minimize
                && direction != (int32_t)ObjectiveDirection::Maximize)
            throw std::invalid_argument("direction = " + std::to_string(direction)
                    + " is not an ObjectiveDirection value.");
        check_count(number_of_quadratic_elements, "number_of_quadratic_elements");

        ObjectiveBatch batch;
        batch.set_direction = (set_direction != 0);
        batch.direction = static_cast<ObjectiveDirection>(direction);
        batch.set_coefficients = (set_coefficients != 0);
        batch.coefficients = coefficients;
        batch.coefficient = coefficient;
        batch.set_quadratic = (set_quadratic != 0);
        batch.number_of_quadratic_elements = (size_t)number_of_quadratic_elements;
        batch.quadratic_variables_1 = quadratic_variables_1;
        batch.quadratic_variables_2 = quadratic_variables_2;
        batch.quadratic_coefficients = quadratic_coefficients;
        batch.set_nonlinear = (set_nonlinear != 0);
        batch.nonlinear = trees_batch(number_of_nonlinear_nodes, nonlinear_operators,
                nonlinear_values, nonlinear_variables, nonlinear_number_of_children);
        set_objective(model->model, batch);
    });
}

/*
 * Transformations.
 */

int mom_dual(
        const mom_model* model,
        mom_model** res,
        char* error_message,
        size_t error_message_size)
{
    return guarded(error_message, error_message_size, [&]() {
        check_model(model);
        if (res == nullptr)
            throw std::invalid_argument("dual must not be null.");
        mom_model* output = new mom_model();
        try {
            output->model = dual(model->model);
        } catch (...) {
            delete output;
            throw;
        }
        *res = output;
    });
}

int mom_expand_defined_variables(
        const mom_model* model,
        mom_model** res,
        char* error_message,
        size_t error_message_size)
{
    return guarded(error_message, error_message_size, [&]() {
        check_model(model);
        if (res == nullptr)
            throw std::invalid_argument("expanded must not be null.");
        mom_model* output = new mom_model();
        try {
            output->model = expand_defined_variables(model->model);
        } catch (...) {
            delete output;
            throw;
        }
        *res = output;
    });
}

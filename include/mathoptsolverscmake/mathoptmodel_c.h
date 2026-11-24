/**
 * C API of MathOptModel: building, reading, writing and inspecting models.
 *
 * It wraps the solver-independent part of the library (mathopt.hpp,
 * mathopt_build.hpp and the MPS / LP / JSON / .nl readers and writers) for
 * other languages, such as Java through the Foreign Function & Memory API.
 *
 * Conventions:
 * - Functions that can fail return a status (MOM_OK on success) and take a
 *   'char* error_message, size_t error_message_size' buffer (which may be
 *   null) that receives a null-terminated error message on failure.
 * - Arrays are passed as pointers plus a count given separately; a null
 *   array means "not given" (see mathopt_build.hpp).
 * - Indices are int32_t (as in MathOptModel); sizes are int64_t.
 * - Strings returned through 'char** data' are allocated by the library and
 *   must be released with mom_free().
 * - A model must not be used concurrently from several threads if one of
 *   them modifies it.
 */
#ifndef MATHOPTMODEL_C_H
#define MATHOPTMODEL_C_H

#include <stddef.h>
#include <stdint.h>

#if defined(_WIN32)
#  if defined(MATHOPTMODEL_C_BUILD)
#    define MOM_API __declspec(dllexport)
#  else
#    define MOM_API __declspec(dllimport)
#  endif
#else
#  define MOM_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct mom_model mom_model;

/** Status codes. */
enum
{
    MOM_OK = 0,
    /** Invalid argument (for example an invalid batch or a malformed file). */
    MOM_INVALID_ARGUMENT = 1,
    /** Any other error. */
    MOM_ERROR = 2,
};

/** File formats. */
enum
{
    MOM_FORMAT_MPS = 0,
    MOM_FORMAT_LP = 1,
    MOM_FORMAT_JSON = 2,
    MOM_FORMAT_NL = 3,
};

/**
 * Model arrays, for mom_array_size() and mom_copy_array(). The element type
 * of each array is given in brackets: [double], [int32] or [int8] (operator
 * char codes, VariableType and ConstraintClass values).
 */
enum
{
    MOM_VARIABLES_LOWER_BOUNDS = 0,                               /* [double] */
    MOM_VARIABLES_UPPER_BOUNDS = 1,                               /* [double] */
    MOM_VARIABLES_TYPES = 2,                                      /* [int8]   */
    MOM_VARIABLES_INITIAL_VALUES = 3,                             /* [double] */
    MOM_CONSTRAINTS_LOWER_BOUNDS = 4,                             /* [double] */
    MOM_CONSTRAINTS_UPPER_BOUNDS = 5,                             /* [double] */
    MOM_CONSTRAINTS_CLASSES = 6,                                  /* [int8]   */
    MOM_OBJECTIVE_COEFFICIENTS = 7,                               /* [double] */
    MOM_CONSTRAINTS_STARTS = 8,                                   /* [int32]  */
    MOM_ELEMENTS_VARIABLES = 9,                                   /* [int32]  */
    MOM_ELEMENTS_COEFFICIENTS = 10,                               /* [double] */
    MOM_OBJECTIVE_QUADRATIC_ELEMENTS_VARIABLES_1 = 11,            /* [int32]  */
    MOM_OBJECTIVE_QUADRATIC_ELEMENTS_VARIABLES_2 = 12,            /* [int32]  */
    MOM_OBJECTIVE_QUADRATIC_ELEMENTS_COEFFICIENTS = 13,           /* [double] */
    MOM_QUADRATIC_ELEMENTS_CONSTRAINTS_STARTS = 14,               /* [int32]  */
    MOM_QUADRATIC_ELEMENTS_VARIABLES_1 = 15,                      /* [int32]  */
    MOM_QUADRATIC_ELEMENTS_VARIABLES_2 = 16,                      /* [int32]  */
    MOM_QUADRATIC_ELEMENTS_COEFFICIENTS = 17,                     /* [double] */
    MOM_OBJECTIVE_NONLINEAR_ELEMENTS_OPERATORS = 18,              /* [int8]   */
    MOM_OBJECTIVE_NONLINEAR_ELEMENTS_VALUES = 19,                 /* [double] */
    MOM_OBJECTIVE_NONLINEAR_ELEMENTS_VARIABLES = 20,              /* [int32]  */
    MOM_OBJECTIVE_NONLINEAR_ELEMENTS_NUMBER_OF_CHILDREN = 21,     /* [int32]  */
    MOM_NONLINEAR_ELEMENTS_CONSTRAINTS_STARTS = 22,               /* [int32]  */
    MOM_NONLINEAR_ELEMENTS_OPERATORS = 23,                        /* [int8]   */
    MOM_NONLINEAR_ELEMENTS_VALUES = 24,                           /* [double] */
    MOM_NONLINEAR_ELEMENTS_VARIABLES = 25,                        /* [int32]  */
    MOM_NONLINEAR_ELEMENTS_NUMBER_OF_CHILDREN = 26,               /* [int32]  */
    MOM_DEFINED_VARIABLES_STARTS = 27,                            /* [int32]  */
    MOM_DEFINED_VARIABLES_ELEMENTS_OPERATORS = 28,                /* [int8]   */
    MOM_DEFINED_VARIABLES_ELEMENTS_VALUES = 29,                   /* [double] */
    MOM_DEFINED_VARIABLES_ELEMENTS_VARIABLES = 30,                /* [int32]  */
    MOM_DEFINED_VARIABLES_ELEMENTS_NUMBER_OF_CHILDREN = 31,       /* [int32]  */
};

/** Name lists, for mom_names_size() and mom_copy_names(). */
enum
{
    MOM_VARIABLES_NAMES = 0,
    MOM_CONSTRAINTS_NAMES = 1,
};

/*
 * Life cycle.
 */

/** Create an empty model (null on allocation failure). */
MOM_API mom_model* mom_model_create(void);

/** Free a model (no-op on null). */
MOM_API void mom_model_free(mom_model* model);

/** Copy a model (null on allocation failure). */
MOM_API mom_model* mom_model_copy(const mom_model* model);

/** Free a buffer returned by the library. */
MOM_API void mom_free(void* data);

/*
 * Reading and writing.
 */

/** Read a model from the contents of a file of the given format. */
MOM_API int mom_read(
        const char* data,
        int64_t size,
        int format,
        mom_model** model,
        char* error_message,
        size_t error_message_size);

/** Write a model in the given format; release '*data' with mom_free(). */
MOM_API int mom_write(
        const mom_model* model,
        int format,
        char** data,
        int64_t* size,
        char* error_message,
        size_t error_message_size);

/** Format the model (MathOptModel::format); release '*data' with mom_free(). */
MOM_API int mom_format(
        const mom_model* model,
        int verbosity_level,
        char** data,
        int64_t* size,
        char* error_message,
        size_t error_message_size);

/** Format a constraint; release '*data' with mom_free(). */
MOM_API int mom_format_constraint(
        const mom_model* model,
        int32_t constraint_id,
        char** data,
        int64_t* size,
        char* error_message,
        size_t error_message_size);

/*
 * Inspection.
 */

MOM_API int32_t mom_number_of_variables(const mom_model* model);
MOM_API int32_t mom_number_of_constraints(const mom_model* model);
MOM_API int32_t mom_number_of_elements(const mom_model* model);
MOM_API int32_t mom_number_of_defined_variables(const mom_model* model);

/** ObjectiveDirection value (0: minimize, 1: maximize). */
MOM_API int32_t mom_objective_direction(const mom_model* model);

/** Number of entries of a model array (MOM_VARIABLES_LOWER_BOUNDS...), -1 if unknown. */
MOM_API int64_t mom_array_size(
        const mom_model* model,
        int array);

/** Copy a model array into 'out', which must hold mom_array_size() entries of its type. */
MOM_API int mom_copy_array(
        const mom_model* model,
        int array,
        void* out,
        char* error_message,
        size_t error_message_size);

/**
 * Number of names of a list (MOM_VARIABLES_NAMES / MOM_CONSTRAINTS_NAMES),
 * 0 if the model has none, and their total size in bytes (UTF-8, not
 * null-terminated) in '*total_size'.
 */
MOM_API int64_t mom_names_size(
        const mom_model* model,
        int names,
        int64_t* total_size);

/**
 * Copy the names of a list, concatenated, into 'data' (total_size bytes),
 * and the offset of each name in 'offsets' (mom_names_size() + 1 entries,
 * the last one total_size).
 */
MOM_API int mom_copy_names(
        const mom_model* model,
        int names,
        char* data,
        int64_t* offsets,
        char* error_message,
        size_t error_message_size);

/** ConstraintSense value of a constraint, -1 if the index is invalid. */
MOM_API int32_t mom_constraint_sense(const mom_model* model, int32_t constraint_id);

/** ConstraintClass value of a constraint, -1 if the index is invalid. */
MOM_API int32_t mom_constraint_class(const mom_model* model, int32_t constraint_id);

/** 1 if the model's arrays are consistent (MathOptModel::check), 0 otherwise. */
MOM_API int32_t mom_check(const mom_model* model, int verbosity_level);

/** Objective value of a solution (number_of_variables entries). */
MOM_API int mom_evaluate_objective(
        const mom_model* model,
        const double* solution,
        double* value,
        char* error_message,
        size_t error_message_size);

/** Value of a constraint's body for a solution (number_of_variables entries). */
MOM_API int mom_evaluate_constraint(
        const mom_model* model,
        const double* solution,
        int32_t constraint_id,
        double* value,
        char* error_message,
        size_t error_message_size);

/** 1 if a solution (number_of_variables entries) is feasible, 0 otherwise. */
MOM_API int mom_check_solution(
        const mom_model* model,
        const double* solution,
        int verbosity_level,
        int32_t* feasible,
        char* error_message,
        size_t error_message_size);

/*
 * Building, in batches (see mathopt_build.hpp for the layout and the
 * validation). '*first' receives the index of the batch's first entry. Names
 * are given as concatenated UTF-8 bytes 'names_data' and 'names_offsets'
 * (number + 1 entries), or null.
 */

MOM_API int mom_add_variables(
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
        size_t error_message_size);

MOM_API int mom_add_constraints(
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
        size_t error_message_size);

MOM_API int mom_add_defined_variables(
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
        size_t error_message_size);

/**
 * Update the objective: each part whose 'set_*' flag is non-zero replaces
 * the current one; the others are kept.
 */
MOM_API int mom_set_objective(
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
        size_t error_message_size);

/*
 * Transformations.
 */

/** Build the dual of a linear program. */
MOM_API int mom_dual(
        const mom_model* model,
        mom_model** dual,
        char* error_message,
        size_t error_message_size);

/** Build an equivalent model without defined variables. */
MOM_API int mom_expand_defined_variables(
        const mom_model* model,
        mom_model** expanded,
        char* error_message,
        size_t error_message_size);

#ifdef __cplusplus
}
#endif

#endif

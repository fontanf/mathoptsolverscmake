#pragma once

/**
 * Building a MathOptModel in batches.
 *
 * Each function appends (or replaces) a whole batch of entries given as raw
 * arrays laid out like MathOptModel's own fields: CSR rows for the linear and
 * quadratic terms of constraints, pre-order element arrays for expression
 * trees (see mathopt.hpp). They are the common implementation of the Python
 * and C APIs' model building.
 *
 * Every batch is fully validated before the model is modified: on an invalid
 * batch, std::invalid_argument is thrown and the model is left unchanged.
 * Argument names in the error messages are the batch fields' names.
 *
 * A null pointer means "not given". For the per-entry arrays of
 * VariablesBatch and ConstraintsBatch (bounds, types, objective
 * coefficients), the corresponding scalar is then used for every entry.
 */

#include "mathoptsolverscmake/mathopt.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace mathoptsolverscmake
{

struct VariablesBatch
{
    /** Number of variables of the batch. */
    std::size_t number = 0;

    const double* lower_bounds = nullptr;
    double lower_bound = -std::numeric_limits<double>::infinity();

    const double* upper_bounds = nullptr;
    double upper_bound = std::numeric_limits<double>::infinity();

    /** VariableType values. */
    const signed char* types = nullptr;
    VariableType type = VariableType::Continuous;

    const double* objective_coefficients = nullptr;
    double objective_coefficient = 0.0;

    /** Initial values (0 for every variable if not given). */
    const double* initial_values = nullptr;

    /** Names (default names "x<index>" if not given). */
    const std::vector<std::string>* names = nullptr;
};

/** Expression trees, concatenated: element arrays of number_of_nodes entries. */
struct ExpressionTreesBatch
{
    std::size_t number_of_nodes = 0;

    /** Operator char codes (see mathopt.hpp). */
    const char* operators = nullptr;

    /** Constant values ('k' nodes); 0 if not given. */
    const double* values = nullptr;

    /** Variable ('v') and defined-variable ('d') indices; -1 if not given. */
    const int* variables = nullptr;

    const int* number_of_children = nullptr;
};

struct ConstraintsBatch
{
    /** Number of constraints of the batch. */
    std::size_t number = 0;

    const double* lower_bounds = nullptr;
    double lower_bound = -std::numeric_limits<double>::infinity();

    const double* upper_bounds = nullptr;
    double upper_bound = std::numeric_limits<double>::infinity();

    /*
     * Linear terms: 'starts' gives the offset of each constraint's first term
     * in 'variables' / 'coefficients' (number entries, the first one 0).
     */

    const int* starts = nullptr;
    std::size_t number_of_elements = 0;
    const int* variables = nullptr;
    const double* coefficients = nullptr;

    /* Quadratic terms, laid out the same way. */

    const int* quadratic_starts = nullptr;
    std::size_t number_of_quadratic_elements = 0;
    const int* quadratic_variables_1 = nullptr;
    const int* quadratic_variables_2 = nullptr;
    const double* quadratic_coefficients = nullptr;

    /*
     * Expression trees: 'nonlinear_starts' gives the offset of each
     * constraint's tree (empty for a constraint without one).
     */

    const int* nonlinear_starts = nullptr;
    ExpressionTreesBatch nonlinear;

    /** Names (default names "c<index>" if not given). */
    const std::vector<std::string>* names = nullptr;
};

struct DefinedVariablesBatch
{
    /** Number of defined variables of the batch. */
    std::size_t number = 0;

    /** Offset of each defined variable's tree in 'trees'. */
    const int* starts = nullptr;

    ExpressionTreesBatch trees;
};

/** Append a batch of variables; return the index of the first one. */
int add_variables(
        MathOptModel& model,
        const VariablesBatch& batch);

/** Append a batch of constraints; return the index of the first one. */
int add_constraints(
        MathOptModel& model,
        const ConstraintsBatch& batch);

/**
 * Append a batch of defined variables; return the index of the first one.
 *
 * A defined variable may only reference defined variables with a lower index.
 */
int add_defined_variables(
        MathOptModel& model,
        const DefinedVariablesBatch& batch);

/**
 * Objective update: each part whose 'set_*' flag is true replaces the current
 * one; the others are kept. An empty quadratic part or expression tree
 * removes it.
 */
struct ObjectiveBatch
{
    bool set_direction = false;
    ObjectiveDirection direction = ObjectiveDirection::Minimize;

    /**
     * Linear coefficients: 'coefficients' has one entry per variable; if
     * null, every coefficient is set to 'coefficient'.
     */
    bool set_coefficients = false;
    const double* coefficients = nullptr;
    double coefficient = 0.0;

    bool set_quadratic = false;
    std::size_t number_of_quadratic_elements = 0;
    const int* quadratic_variables_1 = nullptr;
    const int* quadratic_variables_2 = nullptr;
    const double* quadratic_coefficients = nullptr;

    bool set_nonlinear = false;
    ExpressionTreesBatch nonlinear;
};

/** Update the objective; all parts are validated before any is applied. */
void set_objective(
        MathOptModel& model,
        const ObjectiveBatch& batch);

}

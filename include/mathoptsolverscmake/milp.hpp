#pragma once

#include "mathoptsolverscmake/common.hpp"

#include <vector>
#include <limits>
#include <istream>

#ifdef HIGHS_FOUND
#include "Highs.h"
#endif

#ifdef CBC_FOUND
#include "coin/OsiCbcSolverInterface.hpp"
#endif

#ifdef XPRESS_FOUND
#include "xprs.h"
#endif

namespace mathoptsolverscmake
{

enum class VariableType
{
    Continuous,
    Binary,
    Integer,
};

std::istream& operator>>(
        std::istream& in,
        VariableType& variable_type);

std::ostream& operator<<(
        std::ostream& os,
        VariableType variable_type);

enum class ConstraintSense
{
    LessThanOrEqualTo,
    GreaterThanOrEqualTo,
    Equality,
    Range,
    Free,
};

std::ostream& operator<<(
        std::ostream& os,
        ConstraintSense constraint_sense);

/**
 * Constraint classes.
 *
 * Based on https://miplib.zib.de/statistics.html
 */
enum class ConstraintClass
{
    /** Linear constraint with no variables. */
    Empty,
    /** Linear constraint with no finite side. */
    Free,
    /** Linear constraint with a single variable. */
    Singleton,
    /** Linear constraint of the type ax+by=c. */
    Aggregation,
    /** Linear constraint of the type ax−ay≤b, where x and y must have the same type. */
    Precedence,
    /** Linear constraint of the form ax+by≤c,x∈{0,1}. */
    VariableBound,
    /** Linear constraint of the form ∑xi=1,xi∈{0,1}∀i. */
    SetPartitioning,
    /** Linear constraint of the form ∑xi≤1,xi∈{0,1}∀i. */
    SetPacking,
    /** Linear constraint of the form ∑xi≥1,xi∈{0,1}∀i. */
    SetCovering,
    /** Linear constraint of the form ∑xi=k,xi∈{0,1}∀i,k≥2. */
    Cardinality,
    /** Linear constraint of the form ∑xi≤b,xi∈{0,1}∀i,b∈N≥2. */
    InvariantKnapsack,
    /** Linear constraint of the form ∑aixi=b,xi∈{0,1}∀i,b∈N≥2. */
    EquationKnapsack,
    /** Linear constraint of the form ∑aixi+ax≤a,x,xi∈{0,1}∀i,a∈N≥2. */
    Binpacking,
    /** Linear constraint of the form ∑akxk≤b,xi∈{0,1}∀i,b∈N≥2. */
    Knapsack,
    /** Linear constraint of the form ∑akxk≤b,xi∈Z∀i,b∈N. */
    IntegerKnapsack,
    /** Linear constraint of the form ∑akxk+∑pjsj{≤,=}b,xi∈{0,1}∀i,sj cont. ∀j. */
    MixedBinary,
    /** Linear constraint with no special structure. */
    GeneralLinear,
};

std::ostream& operator<<(
        std::ostream& os,
        ConstraintClass constraint_class);

struct MilpModel
{
    MilpModel(
            int number_of_variables = 0,
            int number_of_constraints = 0,
            int number_of_elements = 0):
        variables_lower_bounds(number_of_variables, -std::numeric_limits<double>::infinity()),
        variables_upper_bounds(number_of_variables, +std::numeric_limits<double>::infinity()),
        variables_types(number_of_variables),
        objective_coefficients(number_of_variables),
        constraints_lower_bounds(number_of_constraints, -std::numeric_limits<double>::infinity()),
        constraints_upper_bounds(number_of_constraints, +std::numeric_limits<double>::infinity()),
        constraints_starts(number_of_constraints),
        elements_variables(number_of_elements),
        elements_coefficients(number_of_elements)
    { }

    /** Get the number of variables. */
    int number_of_variables() const { return variables_lower_bounds.size(); }

    /** Get the number of constraints. */
    int number_of_constraints() const { return constraints_lower_bounds.size(); }

    /** Get the number of elements. */
    int number_of_elements() const { return elements_variables.size(); }

    /** Get the name of a variable. */
    std::string variable_name(int variable_id) const;

    /** Get the name of a constraint. */
    std::string constraint_name(int constraint_id) const;

    /** Get the end element of a constraint. */
    int constraint_end(int constraint_id) const;

    /** Get the number of variables involved in a given constraint. */
    int number_of_variables(int constraint_id) const;

    /** Get the number of variables of a given type in a given constraint. */
    int number_of_variables(
            int constraint_id,
            VariableType variable_type) const;

    /** Get the sense of a given constraint. */
    ConstraintSense constraint_sense(int constraint_id) const;

    /** Find the type of a constraint. */
    ConstraintClass constraint_class(int constraint_id) const;

    /** Format a single constraint. */
    std::ostream& format_constraint(
            std::ostream& os,
            int constraint_id) const;

    /** Format the model. */
    std::ostream& format(
            std::ostream& os,
            int verbosity_level = 1) const;

    /** Format a solution. */
    std::ostream& format_solution(
            std::ostream& os,
            const std::vector<double>& solution,
            int verbosity_level = 1) const;

    /** Get the objective value of a given solution. */
    double evaluate_objective(
            const std::vector<double>& solution) const;

    /** Get the value of a constraint in a given solution. */
    double evaluate_constraint(
            const std::vector<double>& solution,
            int constraint_id) const;

    /** Check if a variable is feasible. */
    bool check_solution_variable(
            const std::vector<double>& solution,
            int variable_id,
            int verbosity_level) const;

    /** Check if a constraint is feasible. */
    bool check_solution_constraint(
            const std::vector<double>& solution,
            int constraint_id,
            int verbosity_level) const;

    /** Check if a solution is feasible. */
    bool check_solution(
            const std::vector<double>& solution,
            int verbosity_level = 0) const;


    ObjectiveDirection objective_direction;

    std::vector<double> variables_lower_bounds;
    std::vector<double> variables_upper_bounds;
    std::vector<VariableType> variables_types;
    std::vector<std::string> variables_names;

    std::vector<double> objective_coefficients;

    std::vector<double> constraints_lower_bounds;
    std::vector<double> constraints_upper_bounds;
    std::vector<int> constraints_starts;
    std::vector<std::string> constraints_names;
    std::vector<int> elements_variables;
    std::vector<double> elements_coefficients;

    double feasiblity_tolerance = 0.0;
    double integrality_tolerance = 0.0;
};

#ifdef CBC_FOUND

void load(
        CbcModel& cbc_model,
        const MilpModel& model);

void reduce_printout(
        CbcModel& cbc_model);

void set_time_limit(
        CbcModel& cbc_model,
        double time_limit);

void solve(
        CbcModel& cbc_model);

double get_solution_value(
        const CbcModel& cbc_model);

std::vector<double> get_solution(
        const CbcModel& cbc_model);

double get_bound(
        const CbcModel& cbc_model);

int get_number_of_nodes(
        const CbcModel& cbc_model);

#endif

#ifdef HIGHS_FOUND

void load(
        Highs& highs_model,
        const MilpModel& model);

void reduce_printout(
        Highs& highs_model);

void set_time_limit(
        Highs& highs_model,
        double time_limit);

void set_log_file(
        Highs& highs_model,
        const std::string& log_file);

void solve(
        Highs& highs_model);

std::vector<double> get_solution(
        const Highs& highs_model);

double get_bound(
        const Highs& highs_model);

#endif

#ifdef XPRESS_FOUND

void load(
        XPRSprob& xpress_model,
        const MilpModel& model);

void set_time_limit(
        XPRSprob& xpress_model,
        double time_limit);

void set_log_file(
        XPRSprob& xpress_model,
        const std::string& log_file);

void write_mps(
        const XPRSprob& xpress_model,
        const std::string& mps_file);

void solve(
        XPRSprob& xpress_model);

double get_solution_value(
        const XPRSprob& xpress_model);

std::vector<double> get_solution(
        const XPRSprob& xpress_model);

double get_bound(
        const XPRSprob& xpress_model);

int get_number_of_nodes(
        const XPRSprob& xpress_model);

#endif

}

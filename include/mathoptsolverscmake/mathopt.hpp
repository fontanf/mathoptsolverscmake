#pragma once

#if defined(__GNUC__) || defined(__clang__)
#  define FUNC_SIGNATURE std::string(__PRETTY_FUNCTION__)
#elif defined(_MSC_VER)
#  define FUNC_SIGNATURE std::string(__FUNCSIG__)
#else
#  define FUNC_SIGNATURE std::string(__func__)
#endif

#include <vector>
#include <limits>
#include <istream>
#include <functional>

namespace mathoptsolverscmake
{

enum class SolverName
{
    Cbc,
    Highs,
    Xpress,
    Knitro,
    Dlib,
    ConicBundle,
};

std::istream& operator>>(
        std::istream& in,
        SolverName& solver_name);

std::ostream& operator<<(
        std::ostream& os,
        SolverName solver_name);

enum class ObjectiveDirection
{
    Minimize,
    Maximize,
};

std::istream& operator>>(
        std::istream& in,
        ObjectiveDirection& objective_direction);

std::ostream& operator<<(
        std::ostream& os,
        ObjectiveDirection objective_direction);

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
    /** Linear quadratic with no special structure. */
    GeneralQuadratic,
    /** Black-box constraint with no special structure. */
    GeneralBlackBox,
    /** Nonlinear constraint expressed as an expression tree. */
    GeneralNonlinear,
};

std::ostream& operator<<(
        std::ostream& os,
        ConstraintClass constraint_class);

struct BlackBoxFunctionOutput
{
    double objective_value = 0;

    std::vector<double> gradient;
};

using BlackBoxFunction = std::function<BlackBoxFunctionOutput(const std::vector<double>&)>;

struct MathOptModel
{
    MathOptModel(
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
        elements_coefficients(number_of_elements),
        quadratic_elements_constraints_starts(number_of_constraints),
        constraints_functions(number_of_constraints)
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

    /** Get the end quadratic element of a constraint. */
    int quadratic_constraint_end(int constraint_id) const;

    /** Get the end nonlinear element of a constraint. */
    int nonlinear_constraint_end(int constraint_id) const;

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

    /*
     * Evaluations
     */

    /** Get the objective value of a given solution. */
    double evaluate_objective(
            const std::vector<double>& solution) const;

    /** Get the value of a constraint in a given solution. */
    double evaluate_constraint(
            const std::vector<double>& solution,
            int constraint_id) const;

    /*
     * Checks
     */

    /** Check if the model is consistent. */
    bool check(int verbosity_level = 0) const;

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

    /*
     * Methods to check what components the model has
     */

    /** Return true if all variables are continuous. */
    bool has_non_continuous_variables() const
    {
        for (const VariableType type: variables_types)
            if (type != VariableType::Continuous)
                return true;
        return false;
    }

    /** Return true if the model has any quadratic terms (objective or constraints). */
    bool has_quadratic() const
    {
        return !this->objective_quadratic_elements_variables_1.empty()
            || !this->quadratic_elements_variables_1.empty();
    }

    /** Return true if the model has any nonlinear expression tree terms (objective or constraints). */
    bool has_nonlinear() const
    {
        return !this->objective_nonlinear_elements_operators.empty()
            || !this->nonlinear_elements_operators.empty();
    }

    /** Return true if the model has any black-box objective or constraint functions. */
    bool has_black_box() const
    {
        if (bool(this->objective_function))
            return true;
        for (const auto& constraint_function: this->constraints_functions)
            if ((bool)constraint_function)
                return true;
        return false;
    }

    /*
     * Method to check the type of model
     */

    /** Return true if the model is a linear program (MILP with only continuous variables). */
    bool is_lp() const
    {
        return !this->has_black_box()
            && !this->has_quadratic()
            && !this->has_nonlinear()
            && !this->has_non_continuous_variables();
    }

    /** Return true if the model is a mixed-integer linear program (no black-box, no quadratic, no nonlinear). An LP is a special case. */
    bool is_milp() const { return !this->has_black_box() && !this->has_quadratic() && !this->has_nonlinear(); }

    /** Return true if the model is only box-constrained. */
    bool is_box_constrained() const { return this->constraints_lower_bounds.empty(); }


    /** Write a solution to a file in HiGHS raw solution format. */
    void write_solution(
            const std::vector<double>& solution,
            const std::string& solution_file) const;


    ObjectiveDirection objective_direction;

    /*
     * Variables
     */

    std::vector<double> variables_lower_bounds;
    std::vector<double> variables_upper_bounds;
    std::vector<VariableType> variables_types;
    std::vector<std::string> variables_names;
    std::vector<double> variables_initial_values;

    /*
     * Constraints
     */

    std::vector<double> constraints_lower_bounds;
    std::vector<double> constraints_upper_bounds;
    std::vector<std::string> constraints_names;

    /*
     * Linear structures
     */

    std::vector<double> objective_coefficients;
    std::vector<int> constraints_starts;
    std::vector<int> elements_variables;
    std::vector<double> elements_coefficients;

    /*
     * Quadratic structures
     */

    std::vector<int> objective_quadratic_elements_variables_1;
    std::vector<int> objective_quadratic_elements_variables_2;
    std::vector<double> objective_quadratic_elements_coefficients;
    std::vector<int> quadratic_elements_constraints_starts;
    std::vector<int> quadratic_elements_variables_1;
    std::vector<int> quadratic_elements_variables_2;
    std::vector<double> quadratic_elements_coefficients;

    /*
     * Nonlinear structures (expression trees stored as parallel element arrays).
     *
     * Element operators use the following char codes:
     *   '+', '-', '*', '/' — binary arithmetic
     *   'n'                — unary negation
     *   'e'                — exp, 'l' — log (natural), 'q' — sqrt
     *   's'                — sin, 'c' — cos, 't' — tan
     *   'p'                — pow (binary)
     *   'v'                — variable (index stored in *_elements_variables)
     *   'k'                — constant  (value  stored in *_elements_values)
     *
     * *_elements_left / *_elements_right hold absolute indices into the element
     * arrays; -1 means no child. By convention the root of each tree is the last
     * element in its range, so no separate root index is needed.
     */

    std::vector<char> objective_nonlinear_elements_operators;
    std::vector<double> objective_nonlinear_elements_values;
    std::vector<int> objective_nonlinear_elements_variables;
    std::vector<int> objective_nonlinear_elements_left;
    std::vector<int> objective_nonlinear_elements_right;

    std::vector<int> nonlinear_elements_constraints_starts;
    std::vector<char> nonlinear_elements_operators;
    std::vector<double> nonlinear_elements_values;
    std::vector<int> nonlinear_elements_variables;
    std::vector<int> nonlinear_elements_left;
    std::vector<int> nonlinear_elements_right;

    /*
     * Black-box structures
     */

    BlackBoxFunction objective_function;
    std::vector<BlackBoxFunction> constraints_functions;

    /*
     * Other parameters
     */

    double feasibility_tolerance = 0.0;
    double integrality_tolerance = 0.0;
};

}

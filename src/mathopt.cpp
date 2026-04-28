#include "mathoptsolverscmake/mathopt.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace mathoptsolverscmake;

std::istream& mathoptsolverscmake::operator>>(
        std::istream& in,
        SolverName& solver_name)
{
    std::string token;
    in >> token;
    if (token == "cbc"
            || token == "Cbc"
            || token == "CBC") {
        solver_name = SolverName::Cbc;
    } else if (token == "highs"
            || token == "Highs"
            || token == "HiGHS"
            || token == "HIGHS") {
        solver_name = SolverName::Highs;
    } else if (token == "xpress"
            || token == "Xpress"
            || token == "XPRESS") {
        solver_name = SolverName::Xpress;
    } else if (token == "Knitro"
            || token == "knitro") {
        solver_name = SolverName::Knitro;
    } else if (token == "dlib") {
        solver_name = SolverName::Dlib;
    } else if (token == "ConicBundle"
            || token == "conicbundle") {
        solver_name = SolverName::ConicBundle;
    } else  {
        in.setstate(std::ios_base::failbit);
    }
    return in;
}

std::ostream& mathoptsolverscmake::operator<<(
        std::ostream& os,
        SolverName solver_name)
{
    switch (solver_name) {
    case SolverName::Cbc: {
        os << "Cbc";
        break;
    } case SolverName::Highs: {
        os << "HiGHS";
        break;
    } case SolverName::Xpress: {
        os << "XPRESS";
        break;
    } case SolverName::Knitro: {
        os << "Knitro";
        break;
    } case SolverName::Dlib: {
        os << "dlib";
        break;
    } case SolverName::ConicBundle: {
        os << "ConicBundle";
        break;
    }
    }
    return os;
}

std::istream& mathoptsolverscmake::operator>>(
        std::istream& in,
        ObjectiveDirection& objective_direction)
{
    std::string token;
    std::getline(in, token);
    if (token == "min"
            || token == "Min"
            || token == "minimize"
            || token == "Minimize") {
        objective_direction = ObjectiveDirection::Minimize;
    } else if (token == "max"
            || token == "Max"
            || token == "maximize"
            || token == "Maximize") {
        objective_direction = ObjectiveDirection::Maximize;
    } else  {
        //throw std::invalid_argument(
        //        FUNC_SIGNATURE + ": "
        //        "invalid input; "
        //        "in: " + token + ".");
        in.setstate(std::ios_base::failbit);
    }
    return in;
}

std::ostream& mathoptsolverscmake::operator<<(
        std::ostream& os,
        ObjectiveDirection objective_direction)
{
    switch (objective_direction) {
    case ObjectiveDirection::Minimize: {
        os << "Minimize";
        break;
    } case ObjectiveDirection::Maximize: {
        os << "Maximize";
        break;
    }
    }
    return os;
}

std::istream& mathoptsolverscmake::operator>>(
        std::istream& in,
        VariableType& variable_type)
{
    std::string token;
    std::getline(in, token);
    if (token == "continuous"
            || token == "Continous"
            || token == "c"
            || token == "C") {
        variable_type = VariableType::Continuous;
    } else if (token == "binary"
            || token == "Binary"
            || token == "b"
            || token == "B") {
        variable_type = VariableType::Binary;
    } else if (token == "integer"
            || token == "Integer"
            || token == "i"
            || token == "I") {
        variable_type = VariableType::Integer;
    } else  {
        //throw std::invalid_argument(
        //        FUNC_SIGNATURE + ": "
        //        "invalid input; "
        //        "in: " + token + ".");
        in.setstate(std::ios_base::failbit);
    }
    return in;
}

std::ostream& mathoptsolverscmake::operator<<(
        std::ostream& os,
        VariableType variable_type)
{
    switch (variable_type) {
    case VariableType::Continuous: {
        os << "Continuous";
        break;
    } case VariableType::Binary: {
        os << "Binary";
        break;
    } case VariableType::Integer: {
        os << "Integer";
        break;
    }
    }
    return os;
}

std::ostream& mathoptsolverscmake::operator<<(
        std::ostream& os,
        ConstraintSense constraint_sense)
{
    switch (constraint_sense) {
    case ConstraintSense::LessThanOrEqualTo: {
        os << "L";
        break;
    } case ConstraintSense::GreaterThanOrEqualTo: {
        os << "G";
        break;
    } case ConstraintSense::Equality: {
        os << "E";
        break;
    } case ConstraintSense::Range: {
        os << "R";
        break;
    } case ConstraintSense::Free: {
        os << "F";
        break;
    }
    }
    return os;
}

std::ostream& mathoptsolverscmake::operator<<(
        std::ostream& os,
        ConstraintClass constraint_class)
{
    switch (constraint_class) {
    case ConstraintClass::Empty: {
        os << "Empty";
        break;
    } case ConstraintClass::Free: {
        os << "Free";
        break;
    } case ConstraintClass::Singleton: {
        os << "Singleton";
        break;
    } case ConstraintClass::Aggregation: {
        os << "Aggregation";
        break;
    } case ConstraintClass::Precedence: {
        os << "Precedence";
        break;
    } case ConstraintClass::VariableBound: {
        os << "Variable bound";
        break;
    } case ConstraintClass::SetPartitioning: {
        os << "Set partitioning";
        break;
    } case ConstraintClass::SetPacking: {
        os << "Set packing";
        break;
    } case ConstraintClass::SetCovering: {
        os << "Set covering";
        break;
    } case ConstraintClass::Cardinality: {
        os << "Cardinality";
        break;
    } case ConstraintClass::InvariantKnapsack: {
        os << "Invariant knapsack";
        break;
    } case ConstraintClass::EquationKnapsack: {
        os << "Equation knapsack";
        break;
    } case ConstraintClass::Binpacking: {
        os << "Bin packing";
        break;
    } case ConstraintClass::Knapsack: {
        os << "Knapsack";
        break;
    } case ConstraintClass::IntegerKnapsack: {
        os << "Integer knapsack";
        break;
    } case ConstraintClass::MixedBinary: {
        os << "Mixed binary";
        break;
    } case ConstraintClass::GeneralLinear: {
        os << "General linear";
        break;
    } case ConstraintClass::GeneralQuadratic: {
        os << "General quadratic";
        break;
    } case ConstraintClass::GeneralBlackBox: {
        os << "General black-box";
        break;
    } case ConstraintClass::GeneralNonlinear: {
        os << "General nonlinear";
        break;
    }
    }
    return os;
}

std::string MathOptModel::variable_name(int variable_id) const
{
    if (this->variables_names.empty())
        return "x" + std::to_string(variable_id);
    return this->variables_names[variable_id];
}

std::string MathOptModel::constraint_name(int constraint_id) const
{
    if (this->constraints_names.empty())
        return "c" + std::to_string(constraint_id);
    return this->constraints_names[constraint_id];
}

int MathOptModel::constraint_end(int constraint_id) const
{
    return (constraint_id == this->number_of_constraints() - 1)?
        this->number_of_elements():
        this->constraints_starts[constraint_id + 1];
}

int MathOptModel::quadratic_constraint_end(int constraint_id) const
{
    return (constraint_id == this->number_of_constraints() - 1)?
        (int)this->quadratic_elements_variables_1.size():
        this->quadratic_elements_constraints_starts[constraint_id + 1];
}

int MathOptModel::nonlinear_constraint_end(int constraint_id) const
{
    return (constraint_id == this->number_of_constraints() - 1)?
        (int)this->nonlinear_elements_operators.size():
        this->nonlinear_elements_constraints_starts[constraint_id + 1];
}

int MathOptModel::number_of_variables(int constraint_id) const
{
    return this->constraint_end(constraint_id) - this->constraints_starts[constraint_id];
}

int MathOptModel::number_of_variables(
        int constraint_id,
        VariableType variable_type) const
{
    int res = 0;
    for (int element_id = this->constraints_starts[constraint_id];
            element_id < this->constraint_end(constraint_id);
            ++element_id) {
        int variable_id = this->elements_variables[element_id];
        if (variable_type == this->variables_types[variable_id])
            res++;
    }
    return res;
}

ConstraintSense MathOptModel::constraint_sense(int constraint_id) const
{
    if (this->constraints_lower_bounds[constraint_id]
            == -std::numeric_limits<double>::infinity()
            && this->constraints_upper_bounds[constraint_id]
            == std::numeric_limits<double>::infinity()) {
        return ConstraintSense::Free;
    }

    if (this->constraints_lower_bounds[constraint_id]
            == -std::numeric_limits<double>::infinity()) {
        return ConstraintSense::LessThanOrEqualTo;
    }

    if (this->constraints_upper_bounds[constraint_id]
            == std::numeric_limits<double>::infinity()) {
        return ConstraintSense::GreaterThanOrEqualTo;
    }

    if (this->constraints_lower_bounds[constraint_id]
            == this->constraints_upper_bounds[constraint_id]) {
        return ConstraintSense::Equality;
    }

    return ConstraintSense::Range;
}

ConstraintClass MathOptModel::constraint_class(int constraint_id) const
{
    if (!this->constraints_functions.empty()
            && this->constraints_functions[constraint_id]) {
        return ConstraintClass::GeneralBlackBox;
    }

    if (!this->nonlinear_elements_operators.empty()
            && this->nonlinear_constraint_end(constraint_id)
            > this->nonlinear_elements_constraints_starts[constraint_id]) {
        return ConstraintClass::GeneralNonlinear;
    }

    if (!this->quadratic_elements_variables_1.empty()
            && this->quadratic_constraint_end(constraint_id)
            > this->quadratic_elements_constraints_starts[constraint_id]) {
        return ConstraintClass::GeneralQuadratic;
    }

    int number_of_variables = this->number_of_variables(constraint_id);
    ConstraintSense constraint_sense = this->constraint_sense(constraint_id);
    double lower_bound = this->constraints_lower_bounds[constraint_id];
    double upper_bound = this->constraints_upper_bounds[constraint_id];

    if (number_of_variables == 0)
        return ConstraintClass::Empty;

    if (constraint_sense == ConstraintSense::Free)
        return ConstraintClass::Free;

    if (number_of_variables == 1)
        return ConstraintClass::Singleton;

    if (constraint_sense == ConstraintSense::Equality
            && number_of_variables == 2) {
        return ConstraintClass::Aggregation;
    }

    if (number_of_variables == 2
            && constraint_sense != ConstraintSense::Range) {
        int constraint_start = this->constraints_starts[constraint_id];
        int variable_1_id = this->elements_variables[constraint_start];
        int variable_2_id = this->elements_variables[constraint_start + 1];
        if (this->variables_types.empty()
                || (this->variables_types[variable_1_id] == this->variables_types[variable_2_id])) {
            return ConstraintClass::Precedence;
        }
    }

    if (number_of_variables == 2
            && constraint_sense != ConstraintSense::Range
            && !this->variables_types.empty()) {
        int constraint_start = this->constraints_starts[constraint_id];
        int variable_1_id = this->elements_variables[constraint_start];
        int variable_2_id = this->elements_variables[constraint_start + 1];
        if (this->variables_types[variable_1_id] == VariableType::Binary
                || this->variables_types[variable_2_id] == VariableType::Binary) {
            return ConstraintClass::VariableBound;
        }
    }

    int number_of_binary_variables = this->number_of_variables(constraint_id, VariableType::Binary);
    bool all_coefficient_one = true;
    for (int element_id = this->constraints_starts[constraint_id];
            element_id < this->constraint_end(constraint_id);
            ++element_id) {
        double coefficient = this->elements_coefficients[element_id];
        if (coefficient != 1.0)
            all_coefficient_one = false;
    }

    if (lower_bound == 1
            && upper_bound == 1
            && all_coefficient_one
            && number_of_binary_variables == number_of_variables) {
        return ConstraintClass::SetPartitioning;
    }

    if (constraint_sense == ConstraintSense::LessThanOrEqualTo
            && upper_bound == 1
            && all_coefficient_one
            && number_of_binary_variables == number_of_variables) {
        return ConstraintClass::SetPacking;
    }

    if (constraint_sense == ConstraintSense::GreaterThanOrEqualTo
            && lower_bound == 1
            && all_coefficient_one
            && number_of_binary_variables == number_of_variables) {
        return ConstraintClass::SetCovering;
    }

    if (constraint_sense == ConstraintSense::Equality
            && all_coefficient_one
            && number_of_binary_variables == number_of_variables) {
        return ConstraintClass::Cardinality;
    }

    if (constraint_sense == ConstraintSense::LessThanOrEqualTo
            && all_coefficient_one
            && number_of_binary_variables == number_of_variables) {
        return ConstraintClass::InvariantKnapsack;
    }

    if (constraint_sense == ConstraintSense::Equality
            && number_of_binary_variables == number_of_variables) {
        return ConstraintClass::EquationKnapsack;
    }

    if (constraint_sense == ConstraintSense::LessThanOrEqualTo
            && number_of_binary_variables == number_of_variables) {
        bool ok = false;
        for (int element_id = this->constraints_starts[constraint_id];
                element_id < this->constraint_end(constraint_id);
                ++element_id) {
            double coefficient = this->elements_coefficients[element_id];
            if (coefficient == upper_bound)
                ok = true;
        }
        if (ok)
            return ConstraintClass::Binpacking;
    }

    if (constraint_sense == ConstraintSense::LessThanOrEqualTo
            && number_of_binary_variables == number_of_variables) {
        return ConstraintClass::Knapsack;
    }

    int number_of_integer_variables = this->number_of_variables(constraint_id, VariableType::Integer);
    if (constraint_sense == ConstraintSense::LessThanOrEqualTo
            && number_of_binary_variables + number_of_integer_variables == number_of_variables) {
        return ConstraintClass::IntegerKnapsack;
    }

    if (number_of_integer_variables == 0)
        return ConstraintClass::MixedBinary;

    return ConstraintClass::GeneralLinear;
}

std::ostream& MathOptModel::format_constraint(
        std::ostream& os,
        int constraint_id) const
{
    os << this->constraint_name(constraint_id) << ": "
        << this->constraints_lower_bounds[constraint_id] << " <=";

    bool first = true;

    if (!this->constraints_functions.empty()
            && bool(this->constraints_functions[constraint_id])) {
        os << " g_" << constraint_id << "(x)";
        first = false;
    }

    if (!this->quadratic_elements_variables_1.empty()) {
        for (int element_id = this->quadratic_elements_constraints_starts[constraint_id];
                element_id < this->quadratic_constraint_end(constraint_id);
                ++element_id) {

            int variable_id_1 = this->quadratic_elements_variables_1[element_id];
            int variable_id_2 = this->quadratic_elements_variables_2[element_id];
            double coefficient = this->quadratic_elements_coefficients[element_id];

            if (first) {
                if (coefficient == 1) {
                    os << " " << this->variable_name(variable_id_1) << " * " << this->variable_name(variable_id_2);
                } else if (coefficient == -1) {
                    os << " -" << this->variable_name(variable_id_1) << " * " << this->variable_name(variable_id_2);
                } else if (coefficient > 0) {
                    os << " " << coefficient << " " << this->variable_name(variable_id_1) << " * " << this->variable_name(variable_id_2);
                } else {
                    os << " -" << -coefficient << " " << this->variable_name(variable_id_1) << " * " << this->variable_name(variable_id_2);
                }
                first = false;
            } else {
                if (coefficient == 1) {
                    os << " + " << this->variable_name(variable_id_1) << " * " << this->variable_name(variable_id_2);
                } else if (coefficient == -1) {
                    os << " - " << this->variable_name(variable_id_1) << " * " << this->variable_name(variable_id_2);
                } else if (coefficient > 0) {
                    os << " + " << coefficient << " " << this->variable_name(variable_id_1) << " * " << this->variable_name(variable_id_2);
                } else {
                    os << " - " << -coefficient << " " << this->variable_name(variable_id_1) << " * " << this->variable_name(variable_id_2);
                }
            }
        }
    }

    if (!this->elements_variables.empty()) {
        for (int element_id = this->constraints_starts[constraint_id];
                element_id < this->constraint_end(constraint_id);
                ++element_id) {

            double variable_id = this->elements_variables[element_id];
            double coefficient = this->elements_coefficients[element_id];

            if (first) {
                // Print constraint lower bound.
                if (coefficient == 1) {
                    os << " " << this->variable_name(variable_id);
                } else if (coefficient == -1) {
                    os << " -" << this->variable_name(variable_id);
                } else {
                    os << " " << coefficient << " " << this->variable_name(variable_id);
                }
                first = false;
            } else {
                if (coefficient == 1) {
                    os << " +" << " " << this->variable_name(variable_id);
                } else if (coefficient == -1) {
                    os << " -" << " " << this->variable_name(variable_id);
                } else if (coefficient > 0) {
                    os << " + " << coefficient << " " << this->variable_name(variable_id);
                } else {
                    os << " - " << -coefficient << " " << this->variable_name(variable_id);
                }
            }
        }
    }

    os << " <= " << this->constraints_upper_bounds[constraint_id] << std::endl;
    return os;
}

std::ostream& MathOptModel::format(
        std::ostream& os,
        int verbosity_level) const
{
    if (verbosity_level == 0)
        return os;

    if (verbosity_level >= 1) {
        os
            << "Number of variables:    " << this->number_of_variables() << std::endl
            << "Number of constraints:  " << this->number_of_constraints() << std::endl
            << "Number of elements:     " << this->number_of_elements() << std::endl
            << "Objective:              " << this->objective_direction << std::endl
            << "Has non-continuous:     " << this->has_non_continuous_variables() << std::endl
            << "Has quadratic:          " << this->has_quadratic() << std::endl
            << "Has nonlinear:          " << this->has_nonlinear() << std::endl
            << "Has black-box:          " << this->has_black_box() << std::endl
            << "Is LP:                  " << this->is_lp() << std::endl
            << "Is MILP:                " << this->is_milp() << std::endl
            << "Is box-constrained:     " << this->is_box_constrained() << std::endl
            ;
    }

    if (verbosity_level >= 2) {
        // Print variables.
        os << std::right << std::endl
            << std::setw(12) << "Variable"
            << std::setw(24) << "Name"
            << std::setw(12) << "Lower"
            << std::setw(12) << "Upper"
            << std::setw(12) << "Type"
            << std::endl
            << std::setw(12) << "--------"
            << std::setw(24) << "----"
            << std::setw(12) << "-----"
            << std::setw(12) << "-----"
            << std::setw(12) << "----"
            << std::endl;
        for (int variable_id = 0;
                variable_id < this->number_of_variables();
                ++variable_id) {
            os
                << std::setw(12) << variable_id
                << std::setw(24) << this->variable_name(variable_id)
                << std::setw(12) << this->variables_lower_bounds[variable_id]
                << std::setw(12) << this->variables_upper_bounds[variable_id]
                << std::setw(12) << this->variables_types[variable_id]
                << std::endl;
        }
        os << std::endl;

        // Print constraints.
        os << std::right << std::endl
            << std::setw(8) << "Constr."
            << std::setw(24) << "Name"
            << std::setw(20) << "Class"
            << std::setw(6) << "Sense"
            << std::setw(12) << "Lower"
            << std::setw(12) << "Upper"
            << std::setw(8) << "# var."
            << std::setw(8) << "# cont."
            << std::setw(8) << "# bin."
            << std::setw(8) << "# int."
            << std::endl
            << std::setw(8) << "-------"
            << std::setw(24) << "----"
            << std::setw(20) << "-----"
            << std::setw(6) << "-----"
            << std::setw(12) << "-----"
            << std::setw(12) << "-----"
            << std::setw(8) << "------"
            << std::setw(8) << "-------"
            << std::setw(8) << "------"
            << std::setw(8) << "------"
            << std::endl;
        for (int constraint_id = 0;
                constraint_id < this->number_of_constraints();
                ++constraint_id) {
            os
                << std::setw(8) << constraint_id
                << std::setw(24) << this->constraint_name(constraint_id)
                << std::setw(20) << this->constraint_class(constraint_id)
                << std::setw(6) << this->constraint_sense(constraint_id)
                << std::setw(12) << this->constraints_lower_bounds[constraint_id]
                << std::setw(12) << this->constraints_upper_bounds[constraint_id]
                << std::setw(8) << this->number_of_variables(constraint_id)
                << std::setw(8) << this->number_of_variables(constraint_id, VariableType::Continuous)
                << std::setw(8) << this->number_of_variables(constraint_id, VariableType::Binary)
                << std::setw(8) << this->number_of_variables(constraint_id, VariableType::Integer)
                << std::endl;
        }
        os << std::endl;
    }

    if (verbosity_level >= 3) {
        // Print objective.
        os << "obj:";
        bool first = true;
        if (bool(this->objective_function)) {
            os << " f(x)";
            first = false;
        }
        for (int element_id = 0;
                element_id < (int)this->objective_quadratic_elements_variables_1.size();
                ++element_id) {
            int variable_id_1 = this->objective_quadratic_elements_variables_1[element_id];
            int variable_id_2 = this->objective_quadratic_elements_variables_2[element_id];
            double coefficient = this->objective_quadratic_elements_coefficients[element_id];
            if (first) {
                if (coefficient == 1) {
                    os << " " << this->variable_name(variable_id_1) << " * " << this->variable_name(variable_id_2);
                } else if (coefficient == -1) {
                    os << " -" << this->variable_name(variable_id_1) << " * " << this->variable_name(variable_id_2);
                } else if (coefficient > 0) {
                    os << " " << coefficient << " " << this->variable_name(variable_id_1) << " * " << this->variable_name(variable_id_2);
                } else {
                    os << " -" << -coefficient << " " << this->variable_name(variable_id_1) << " * " << this->variable_name(variable_id_2);
                }
                first = false;
            } else {
                if (coefficient == 1) {
                    os << " + " << this->variable_name(variable_id_1) << " * " << this->variable_name(variable_id_2);
                } else if (coefficient == -1) {
                    os << " - " << this->variable_name(variable_id_1) << " * " << this->variable_name(variable_id_2);
                } else if (coefficient > 0) {
                    os << " + " << coefficient << " " << this->variable_name(variable_id_1) << " * " << this->variable_name(variable_id_2);
                } else {
                    os << " - " << -coefficient << " " << this->variable_name(variable_id_1) << " * " << this->variable_name(variable_id_2);
                }
            }
        }
        for (int variable_id = 0;
                variable_id < this->number_of_variables();
                ++variable_id) {
            double coefficient = this->objective_coefficients[variable_id];
            if (coefficient == 0)
                continue;
            if (first) {
                if (coefficient == 1) {
                    os << " " << this->variable_name(variable_id);
                } else if (coefficient == -1) {
                    os << " -" << this->variable_name(variable_id);
                } else if (coefficient > 0) {
                    os << " " << coefficient << " " << this->variable_name(variable_id);
                } else {
                    os << " -" << -coefficient << " " << this->variable_name(variable_id);
                }
                first = false;
            } else {
                if (coefficient == 0) {
                } else if (coefficient == 1) {
                    os << " + " << this->variable_name(variable_id);
                } else if (coefficient == -1) {
                    os << " - " << this->variable_name(variable_id);
                } else if (coefficient > 0) {
                    os << " + " << coefficient << " " << this->variable_name(variable_id);
                } else {
                    os << " - " << -coefficient << " " << this->variable_name(variable_id);
                }
            }
        }
        os << std::endl;

        // Print matrix.
        for (int constraint_id = 0;
                constraint_id < this->number_of_constraints();
                ++constraint_id) {
            this->format_constraint(os, constraint_id);
        }
    }
    return os;
}

std::ostream& MathOptModel::format_solution(
        std::ostream& os,
        const std::vector<double>& solution,
        int verbosity_level) const
{
    if (verbosity_level == 0)
        return os;

    if (verbosity_level >= 1) {
        os
            << "Objective:  " << this->evaluate_objective(solution) << std::endl
            << "Feasible:   " << this->check_solution(solution) << std::endl
            ;
    }

    if (verbosity_level >= 2) {
        // Print variables.
        os << std::right << std::endl
            << std::setw(8) << "Var."
            << std::setw(24) << "Name"
            << std::setw(12) << "Type"
            << std::setw(12) << "Lower"
            << std::setw(12) << "Value"
            << std::setw(12) << "Upper"
            << std::setw(8) << "Feas."
            << std::endl
            << std::setw(8) << "----"
            << std::setw(24) << "----"
            << std::setw(12) << "-----"
            << std::setw(12) << "-----"
            << std::setw(12) << "----"
            << std::setw(8) << "-----"
            << std::endl;
        for (int variable_id = 0;
                variable_id < this->number_of_variables();
                ++variable_id) {
            os
                << std::setw(8) << variable_id
                << std::setw(24) << this->variable_name(variable_id)
                << std::setw(12) << this->variables_types[variable_id]
                << std::setw(12) << this->variables_lower_bounds[variable_id]
                << std::setw(12) << solution[variable_id]
                << std::setw(12) << this->variables_upper_bounds[variable_id]
                << std::setw(8) << check_solution_variable(solution, variable_id, 0)
                << std::endl;
        }

        // Print constraints.
        os << std::right << std::endl
            << std::setw(8) << "Constr."
            << std::setw(24) << "Name"
            << std::setw(12) << "Lower"
            << std::setw(12) << "Value"
            << std::setw(12) << "Upper"
            << std::setw(8) << "Feas."
            << std::endl
            << std::setw(8) << "-------"
            << std::setw(24) << "----"
            << std::setw(12) << "-----"
            << std::setw(12) << "-----"
            << std::setw(12) << "-----"
            << std::setw(8) << "-----"
            << std::endl;
        for (int constraint_id = 0;
                constraint_id < this->number_of_constraints();
                ++constraint_id) {
            os
                << std::setw(8) << constraint_id
                << std::setw(24) << this->constraint_name(constraint_id)
                << std::setw(12) << this->constraints_lower_bounds[constraint_id]
                << std::setw(12) << this->evaluate_constraint(solution, constraint_id)
                << std::setw(12) << this->constraints_upper_bounds[constraint_id]
                << std::setw(8) << this->check_solution_constraint(solution, constraint_id, 0)
                << std::endl;
        }
        os << std::endl;
    }

    return os;
}

static double evaluate_nonlinear_element(
        const std::vector<char>& operators,
        const std::vector<double>& values,
        const std::vector<int>& variables,
        const std::vector<int>& parent,
        const std::vector<double>& solution,
        int start,
        int end)
{
    const int n = end - start;

    // Arity is needed for n-ary + and *; derive it from the parent array.
    std::vector<int> arity(n, 0);
    for (int i = start + 1; i < end; ++i)
        arity[parent[i] - start]++;

    // Evaluate by scanning from the last leaf back to the root (index start).
    // Pre-order storage means this is a post-order traversal: every child is
    // processed before its parent.  The left subtree ends up on top of the
    // stack (lower indices, therefore processed last), so for binary operators
    // first-pop == left operand and second-pop == right operand.
    std::vector<double> stk;
    stk.reserve(n);
    for (int element_id = end - 1; element_id >= start; --element_id) {
        const int i = element_id - start;
        switch (operators[element_id]) {
        case 'k': stk.push_back(values[element_id]); break;
        case 'v': stk.push_back(solution[variables[element_id]]); break;
        case 'n': { double a = stk.back(); stk.pop_back(); stk.push_back(-a); break; }
        case 'e': { double a = stk.back(); stk.pop_back(); stk.push_back(std::exp(a)); break; }
        case 'l': { double a = stk.back(); stk.pop_back(); stk.push_back(std::log(a)); break; }
        case 'q': { double a = stk.back(); stk.pop_back(); stk.push_back(std::sqrt(a)); break; }
        case 's': { double a = stk.back(); stk.pop_back(); stk.push_back(std::sin(a)); break; }
        case 'c': { double a = stk.back(); stk.pop_back(); stk.push_back(std::cos(a)); break; }
        case 't': { double a = stk.back(); stk.pop_back(); stk.push_back(std::tan(a)); break; }
        case '+': {
            double sum = 0;
            for (int c = 0; c < arity[i]; ++c) { sum += stk.back(); stk.pop_back(); }
            stk.push_back(sum);
            break;
        }
        case '*': {
            double prod = 1;
            for (int c = 0; c < arity[i]; ++c) { prod *= stk.back(); stk.pop_back(); }
            stk.push_back(prod);
            break;
        }
        case '-': { double l = stk.back(); stk.pop_back(); double r = stk.back(); stk.pop_back(); stk.push_back(l - r); break; }
        case '/': { double l = stk.back(); stk.pop_back(); double r = stk.back(); stk.pop_back(); stk.push_back(l / r); break; }
        case 'p': { double l = stk.back(); stk.pop_back(); double r = stk.back(); stk.pop_back(); stk.push_back(std::pow(l, r)); break; }
        default:  stk.push_back(0.0); break;
        }
    }
    return stk.back();
}

double MathOptModel::evaluate_objective(
        const std::vector<double>& solution) const
{
    double value = 0.0;
    if (!this->objective_coefficients.empty()) {
        for (int variable_id = 0;
                variable_id < this->number_of_variables();
                ++variable_id) {
            value += solution[variable_id] * this->objective_coefficients[variable_id];
        }
    }
    if (!this->objective_quadratic_elements_variables_1.empty()) {
        for (int element_id = 0;
                element_id < (int)this->objective_quadratic_elements_variables_1.size();
                ++element_id) {
            int variable_1_id = this->objective_quadratic_elements_variables_1[element_id];
            int variable_2_id = this->objective_quadratic_elements_variables_2[element_id];
            double coefficient = this->objective_quadratic_elements_coefficients[element_id];
            value += coefficient * solution[variable_1_id] * solution[variable_2_id];
        }
    }
    if (!this->objective_nonlinear_elements_operators.empty()) {
        value += evaluate_nonlinear_element(
                this->objective_nonlinear_elements_operators,
                this->objective_nonlinear_elements_values,
                this->objective_nonlinear_elements_variables,
                this->objective_nonlinear_elements_parent,
                solution,
                0,
                (int)this->objective_nonlinear_elements_operators.size());
    }
    if (this->objective_function)
        value += this->objective_function(solution).objective_value;
    return value;
}

double MathOptModel::evaluate_constraint(
        const std::vector<double>& solution,
        int constraint_id) const
{
    double value = 0.0;
    double compensation = 0.0;
    if (!this->elements_variables.empty()) {
        for (int element_id = this->constraints_starts[constraint_id];
                element_id < this->constraint_end(constraint_id);
                ++element_id) {
            double variable_id = this->elements_variables[element_id];
            double coefficient = this->elements_coefficients[element_id];
            double term = solution[variable_id] * coefficient - compensation;
            double new_value = value + term;
            compensation = (new_value - value) - term;
            value = new_value;
        }
    }
    if (!this->quadratic_elements_variables_1.empty()) {
        for (int element_id = this->quadratic_elements_constraints_starts[constraint_id];
                element_id < this->quadratic_constraint_end(constraint_id);
                ++element_id) {
            int variable_1_id = this->quadratic_elements_variables_1[element_id];
            int variable_2_id = this->quadratic_elements_variables_2[element_id];
            double coefficient = this->quadratic_elements_coefficients[element_id];
            value += coefficient * solution[variable_1_id] * solution[variable_2_id];
        }
    }
    if (!this->nonlinear_elements_operators.empty()) {
        int start = this->nonlinear_elements_constraints_starts[constraint_id];
        int end = this->nonlinear_constraint_end(constraint_id);
        if (end > start) {
            value += evaluate_nonlinear_element(
                    this->nonlinear_elements_operators,
                    this->nonlinear_elements_values,
                    this->nonlinear_elements_variables,
                    this->nonlinear_elements_parent,
                    solution,
                    start,
                    end);
        }
    }
    if (!this->constraints_functions.empty()
            && this->constraints_functions[constraint_id]) {
        value += this->constraints_functions[constraint_id](solution).objective_value;
    }
    return value;
}

bool MathOptModel::check_solution_variable(
        const std::vector<double>& solution,
        int variable_id,
        int verbosity_level) const
{
    bool feasible = true;
    double value = solution[variable_id];
    if (value < this->variables_lower_bounds[variable_id] - this->feasibility_tolerance) {
        if (verbosity_level > 0) {
            std::stringstream ss;
            ss << "violated variable lower bound; "
                << "variable_id: " << variable_id << "; "
                << "name: " << this->variable_name(variable_id) << "; "
                << "lower_bound: " << this->variables_lower_bounds[variable_id] << "; "
                << "value: " << value << ".";
            std::cout << ss.str() << std::endl;
        }
        feasible = false;
    }
    if (value > this->variables_upper_bounds[variable_id] + this->feasibility_tolerance) {
        if (verbosity_level > 0) {
            std::stringstream ss;
            ss << "violated variable upper bound; "
                << "variable_id: " << variable_id << "; "
                << "name: " << this->variable_name(variable_id) << "; "
                << "upper_bound: " << this->variables_upper_bounds[variable_id] << "; "
                << "value: " << value << ".";
            std::cout << ss.str() << std::endl;
        }
        feasible = false;
    }
    double fractionality = std::abs(value - std::round(value));
    if (this->variables_types[variable_id] != VariableType::Continuous
            && fractionality > this->integrality_tolerance) {
        if (verbosity_level > 0) {
            std::stringstream ss;
            ss << "violated variable integrality; "
                << "variable_id: " << variable_id << "; "
                << "name: " << this->variable_name(variable_id) << "; "
                << "type: " << this->variables_types[variable_id] << "; "
                << "value: " << value << ".";
            std::cout << ss.str() << std::endl;
        }
        feasible = false;
    }
    return feasible;
}

bool MathOptModel::check_solution_constraint(
        const std::vector<double>& solution,
        int constraint_id,
        int verbosity_level) const
{
    bool feasible = true;
    double value = this->evaluate_constraint(solution, constraint_id);
    if (value < this->constraints_lower_bounds[constraint_id] - this->feasibility_tolerance) {
        if (verbosity_level > 0) {
            std::stringstream ss;
            this->format_constraint(ss, constraint_id);
            ss << "violated constraint lower bound; "
                << "constraint_id: " << constraint_id << "; "
                << "name: " << this->constraint_name(constraint_id) << "; "
                << "lower_bound: " << this->constraints_lower_bounds[constraint_id] << "; "
                << "value: " << value << ".";
            std::cout << ss.str() << std::endl;
        }
        feasible = false;
    }
    if (value > this->constraints_upper_bounds[constraint_id] + this->feasibility_tolerance) {
        if (verbosity_level > 0) {
            std::stringstream ss;
            this->format_constraint(ss, constraint_id);
            ss << "violated constraint upper bound; "
                << "constraint_id: " << constraint_id << "; "
                << "name: " << this->constraint_name(constraint_id) << "; "
                << "upper_bound: " << this->constraints_upper_bounds[constraint_id] << "; "
                << "value: " << value << ".";
            std::cout << ss.str() << std::endl;
        }
        feasible = false;
    }
    return feasible;
}

bool MathOptModel::check_solution(
        const std::vector<double>& solution,
        int verbosity_level) const
{
    bool feasible = true;

    // Check solution size.
    if (solution.size() != this->number_of_variables()) {
        if (verbosity_level > 0) {
            std::stringstream ss;
            ss << "wrong solution size; "
                << "solution.size(): " << solution.size() << "; "
                << "number_of_variables(): " << this->number_of_variables() << ".";
            std::cout << ss.str() << std::endl;
        }
        feasible = false;
        return feasible;
    }

    // Check variable bounds and integrity.
    for (int variable_id = 0;
            variable_id < this->number_of_variables();
            ++variable_id) {
        feasible &= check_solution_variable(
                solution,
                variable_id,
                verbosity_level);
    }

    for (int constraint_id = 0;
            constraint_id < this->number_of_constraints();
            ++constraint_id) {
        feasible &= check_solution_constraint(
                solution,
                constraint_id,
                verbosity_level);
    }

    return feasible;
}

bool MathOptModel::check(int verbosity_level) const
{
    bool ok = true;

    // Check variable vector sizes.
    if ((int)this->variables_upper_bounds.size() != this->number_of_variables()) {
        if (verbosity_level > 0) {
            std::cout << "inconsistent variables_upper_bounds size: "
                << this->variables_upper_bounds.size() << " != " << this->number_of_variables() << "." << std::endl;
        }
        ok = false;
    }
    if ((int)this->variables_types.size() != this->number_of_variables()) {
        if (verbosity_level > 0) {
            std::cout << "inconsistent variables_types size: "
                << this->variables_types.size() << " != " << this->number_of_variables() << "." << std::endl;
        }
        ok = false;
    }
    if ((int)this->objective_coefficients.size() != this->number_of_variables()) {
        if (verbosity_level > 0) {
            std::cout << "inconsistent objective_coefficients size: "
                << this->objective_coefficients.size() << " != " << this->number_of_variables() << "." << std::endl;
        }
        ok = false;
    }
    if (!this->variables_names.empty() && (int)this->variables_names.size() != this->number_of_variables()) {
        if (verbosity_level > 0) {
            std::cout << "inconsistent variables_names size: "
                << this->variables_names.size() << " != " << this->number_of_variables() << "." << std::endl;
        }
        ok = false;
    }

    // Check constraint vector sizes.
    if ((int)this->constraints_upper_bounds.size() != this->number_of_constraints()) {
        if (verbosity_level > 0) {
            std::cout << "inconsistent constraints_upper_bounds size: "
                << this->constraints_upper_bounds.size() << " != " << this->number_of_constraints() << "." << std::endl;
        }
        ok = false;
    }
    if ((int)this->constraints_starts.size() != this->number_of_constraints()) {
        if (verbosity_level > 0) {
            std::cout << "inconsistent constraints_starts size: "
                << this->constraints_starts.size() << " != " << this->number_of_constraints() << "." << std::endl;
        }
        ok = false;
    }
    if (!this->constraints_names.empty() && (int)this->constraints_names.size() != this->number_of_constraints()) {
        if (verbosity_level > 0) {
            std::cout << "inconsistent constraints_names size: "
                << this->constraints_names.size() << " != " << this->number_of_constraints() << "." << std::endl;
        }
        ok = false;
    }

    // Check element vector sizes.
    if ((int)this->elements_coefficients.size() != this->number_of_elements()) {
        if (verbosity_level > 0) {
            std::cout << "inconsistent elements_coefficients size: "
                << this->elements_coefficients.size() << " != " << this->number_of_elements() << "." << std::endl;
        }
        ok = false;
    }

    // Check quadratic objective vector sizes.
    if (this->objective_quadratic_elements_variables_2.size()
            != this->objective_quadratic_elements_variables_1.size()) {
        if (verbosity_level > 0) {
            std::cout << "inconsistent objective_quadratic_elements_variables_2 size: "
                << this->objective_quadratic_elements_variables_2.size()
                << " != " << this->objective_quadratic_elements_variables_1.size() << "." << std::endl;
        }
        ok = false;
    }
    if (this->objective_quadratic_elements_coefficients.size()
            != this->objective_quadratic_elements_variables_1.size()) {
        if (verbosity_level > 0) {
            std::cout << "inconsistent objective_quadratic_elements_coefficients size: "
                << this->objective_quadratic_elements_coefficients.size()
                << " != " << this->objective_quadratic_elements_variables_1.size() << "." << std::endl;
        }
        ok = false;
    }

    // Check quadratic constraint vector sizes.
    if (!this->quadratic_elements_variables_1.empty()) {
        if ((int)this->quadratic_elements_constraints_starts.size() != this->number_of_constraints()) {
            if (verbosity_level > 0) {
                std::cout << "inconsistent quadratic_elements_constraints_starts size: "
                    << this->quadratic_elements_constraints_starts.size()
                    << " != " << this->number_of_constraints() << "." << std::endl;
            }
            ok = false;
        }
        if (this->quadratic_elements_variables_2.size() != this->quadratic_elements_variables_1.size()) {
            if (verbosity_level > 0) {
                std::cout << "inconsistent quadratic_elements_variables_2 size: "
                    << this->quadratic_elements_variables_2.size()
                    << " != " << this->quadratic_elements_variables_1.size() << "." << std::endl;
            }
            ok = false;
        }
        if (this->quadratic_elements_coefficients.size() != this->quadratic_elements_variables_1.size()) {
            if (verbosity_level > 0) {
                std::cout << "inconsistent quadratic_elements_coefficients size: "
                    << this->quadratic_elements_coefficients.size()
                    << " != " << this->quadratic_elements_variables_1.size() << "." << std::endl;
            }
            ok = false;
        }
    }

    // Check black-box constraint vector size.
    if (!this->constraints_functions.empty()
            && (int)this->constraints_functions.size() != this->number_of_constraints()) {
        if (verbosity_level > 0) {
            std::cout << "inconsistent constraints_functions size: "
                << this->constraints_functions.size()
                << " != " << this->number_of_constraints() << "." << std::endl;
        }
        ok = false;
    }

    // Check variables_initial_values size.
    if (!this->variables_initial_values.empty()
            && (int)this->variables_initial_values.size() != this->number_of_variables()) {
        if (verbosity_level > 0) {
            std::cout << "inconsistent variables_initial_values size: "
                << this->variables_initial_values.size()
                << " != " << this->number_of_variables() << "." << std::endl;
        }
        ok = false;
    }

    // Check nonlinear objective element vector sizes.
    if (this->objective_nonlinear_elements_values.size()
            != this->objective_nonlinear_elements_operators.size()) {
        if (verbosity_level > 0)
            std::cout << "inconsistent objective_nonlinear_elements_values size: "
                << this->objective_nonlinear_elements_values.size()
                << " != " << this->objective_nonlinear_elements_operators.size() << "." << std::endl;
        ok = false;
    }
    if (this->objective_nonlinear_elements_variables.size()
            != this->objective_nonlinear_elements_operators.size()) {
        if (verbosity_level > 0)
            std::cout << "inconsistent objective_nonlinear_elements_variables size: "
                << this->objective_nonlinear_elements_variables.size()
                << " != " << this->objective_nonlinear_elements_operators.size() << "." << std::endl;
        ok = false;
    }
    if (this->objective_nonlinear_elements_parent.size()
            != this->objective_nonlinear_elements_operators.size()) {
        if (verbosity_level > 0)
            std::cout << "inconsistent objective_nonlinear_elements_parent size: "
                << this->objective_nonlinear_elements_parent.size()
                << " != " << this->objective_nonlinear_elements_operators.size() << "." << std::endl;
        ok = false;
    }

    // Check nonlinear constraint element vector sizes.
    if (!this->nonlinear_elements_operators.empty()) {
        if ((int)this->nonlinear_elements_constraints_starts.size() != this->number_of_constraints()) {
            if (verbosity_level > 0)
                std::cout << "inconsistent nonlinear_elements_constraints_starts size: "
                    << this->nonlinear_elements_constraints_starts.size()
                    << " != " << this->number_of_constraints() << "." << std::endl;
            ok = false;
        }
        if (this->nonlinear_elements_values.size() != this->nonlinear_elements_operators.size()) {
            if (verbosity_level > 0)
                std::cout << "inconsistent nonlinear_elements_values size: "
                    << this->nonlinear_elements_values.size()
                    << " != " << this->nonlinear_elements_operators.size() << "." << std::endl;
            ok = false;
        }
        if (this->nonlinear_elements_variables.size() != this->nonlinear_elements_operators.size()) {
            if (verbosity_level > 0)
                std::cout << "inconsistent nonlinear_elements_variables size: "
                    << this->nonlinear_elements_variables.size()
                    << " != " << this->nonlinear_elements_operators.size() << "." << std::endl;
            ok = false;
        }
        if (this->nonlinear_elements_parent.size() != this->nonlinear_elements_operators.size()) {
            if (verbosity_level > 0)
                std::cout << "inconsistent nonlinear_elements_parent size: "
                    << this->nonlinear_elements_parent.size()
                    << " != " << this->nonlinear_elements_operators.size() << "." << std::endl;
            ok = false;
        }
    }

    // Stop here if sizes are wrong to avoid out-of-bounds access below.
    if (!ok)
        return false;

    // Check constraints_starts validity.
    for (int constraint_id = 0;
            constraint_id < this->number_of_constraints();
            ++constraint_id) {
        int start = this->constraints_starts[constraint_id];
        if (start < 0 || start > this->number_of_elements()) {
            if (verbosity_level > 0) {
                std::cout << "constraints_starts[" << constraint_id << "] = " << start
                    << " is out of range [0, " << this->number_of_elements() << "]." << std::endl;
            }
            ok = false;
        }
        if (constraint_id > 0 && start < this->constraints_starts[constraint_id - 1]) {
            if (verbosity_level > 0) {
                std::cout << "constraints_starts is not non-decreasing at constraint "
                    << constraint_id << "." << std::endl;
            }
            ok = false;
        }
    }

    // Check elements_variables validity and duplicate variables per constraint.
    for (int constraint_id = 0; constraint_id < this->number_of_constraints(); ++constraint_id) {
        int start = this->constraints_starts[constraint_id];
        int end = this->constraint_end(constraint_id);
        std::vector<bool> seen(this->number_of_variables(), false);
        for (int element_id = start; element_id < end; ++element_id) {
            int variable_id = this->elements_variables[element_id];
            if (variable_id < 0 || variable_id >= this->number_of_variables()) {
                if (verbosity_level > 0) {
                    std::cout << "elements_variables[" << element_id << "] = " << variable_id
                        << " is out of range [0, " << this->number_of_variables() << ") in constraint "
                        << constraint_id << "." << std::endl;
                }
                ok = false;
                continue;
            }
            if (seen[variable_id]) {
                if (verbosity_level > 0) {
                    std::cout << "variable " << variable_id
                        << " appears more than once in constraint " << constraint_id << "." << std::endl;
                }
                ok = false;
            }
            seen[variable_id] = true;
        }
    }

    return ok;
}

void MathOptModel::write_solution(
        const std::vector<double>& solution,
        const std::string& solution_file) const
{
    std::ofstream file(solution_file);
    file << "Model status\n"
         << "Feasible\n"
         << "\n"
         << "# Primal solution values\n"
         << "Feasible\n";
    file << std::setprecision(17);
    file << "Objective " << this->evaluate_objective(solution) << "\n";
    file << "# Columns " << this->number_of_variables() << "\n";
    for (int variable_id = 0;
            variable_id < this->number_of_variables();
            ++variable_id) {
        file << this->variable_name(variable_id)
             << " " << solution[variable_id] << "\n";
    }
    file << "# Rows " << this->number_of_constraints() << "\n";
    for (int constraint_id = 0;
            constraint_id < this->number_of_constraints();
            ++constraint_id) {
        file << this->constraint_name(constraint_id)
             << " " << this->evaluate_constraint(solution, constraint_id) << "\n";
    }
}

#include "mathoptsolverscmake/milp.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <numeric>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace mathoptsolverscmake;

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
    }
    }
    return os;
}

std::string MilpModel::variable_name(int variable_id) const
{
    if (this->variables_names.empty())
        return "x" + std::to_string(variable_id);
    return this->variables_names[variable_id];
}

std::string MilpModel::constraint_name(int constraint_id) const
{
    if (this->constraints_names.empty())
        return "c" + std::to_string(constraint_id);
    return this->constraints_names[constraint_id];
}

int MilpModel::constraint_end(int constraint_id) const
{
    return (constraint_id == this->number_of_constraints() - 1)?
        this->number_of_elements():
        this->constraints_starts[constraint_id + 1];
}

int MilpModel::number_of_variables(int constraint_id) const
{
    return this->constraint_end(constraint_id) - this->constraints_starts[constraint_id] - 1;
}

int MilpModel::number_of_variables(
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

ConstraintSense MilpModel::constraint_sense(int constraint_id) const
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

ConstraintClass MilpModel::constraint_class(int constraint_id) const
{
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
            && number_of_variables) {
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

std::ostream& MilpModel::format_constraint(
        std::ostream& os,
        int constraint_id) const
{
    os << this->constraint_name(constraint_id) << ": "
        << this->constraints_lower_bounds[constraint_id] << " <=";

    bool first = true;
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
    os << " <= " << this->constraints_upper_bounds[constraint_id] << std::endl;
    return os;
}

std::ostream& MilpModel::format(
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

std::ostream& MilpModel::format_solution(
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

double MilpModel::evaluate_objective(
        const std::vector<double>& solution) const
{
    double value = 0.0;
    for (int variable_id = 0;
            variable_id < this->number_of_variables();
            ++variable_id) {
        double coefficient = this->objective_coefficients[variable_id];
        value += solution[variable_id] * coefficient;
    }
    return value;
}

double MilpModel::evaluate_constraint(
        const std::vector<double>& solution,
        int constraint_id) const
{
    double value = 0.0;
    for (int element_id = this->constraints_starts[constraint_id];
            element_id < this->constraint_end(constraint_id);
            ++element_id) {
        double variable_id = this->elements_variables[element_id];
        double coefficient = this->elements_coefficients[element_id];
        value += solution[variable_id] * coefficient;
    }
    return value;
}

bool MilpModel::check_solution_variable(
        const std::vector<double>& solution,
        int variable_id,
        int verbosity_level) const
{
    bool feasible = true;
    double value = solution[variable_id];
    if (value < this->variables_lower_bounds[variable_id] - this->feasiblity_tolerance) {
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
    if (value > this->variables_upper_bounds[variable_id] + this->feasiblity_tolerance) {
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

bool MilpModel::check_solution_constraint(
        const std::vector<double>& solution,
        int constraint_id,
        int verbosity_level) const
{
    bool feasible = true;
    double value = this->evaluate_constraint(solution, constraint_id);
    if (value < this->constraints_lower_bounds[constraint_id] - this->feasiblity_tolerance) {
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
    if (value > this->constraints_upper_bounds[constraint_id] + this->feasiblity_tolerance) {
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

bool MilpModel::check_solution(
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

#ifdef CBC_FOUND

void mathoptsolverscmake::load(
        CbcModel& cbc_model,
        const MilpModel& model)
{
    std::vector<int> constraints_lengths(model.number_of_constraints(), 0);
    for (int constraint_id = 0;
            constraint_id < model.number_of_constraints() - 1;
            ++constraint_id) {
        constraints_lengths[constraint_id] = model.constraints_starts[constraint_id + 1] - model.constraints_starts[constraint_id];
    }
    constraints_lengths[model.number_of_constraints() - 1] = model.number_of_elements() - model.constraints_starts[model.number_of_constraints() - 1];
    CoinPackedMatrix matrix(
            false,
            model.number_of_variables(),
            model.number_of_constraints(),
            model.number_of_elements(),
            model.elements_coefficients.data(),
            model.elements_variables.data(),
            model.constraints_starts.data(),
            constraints_lengths.data());
    cbc_model.solver()->loadProblem(
            matrix,
            model.variables_lower_bounds.data(),
            model.variables_upper_bounds.data(),
            model.objective_coefficients.data(),
            model.constraints_lower_bounds.data(),
            model.constraints_upper_bounds.data());

    // Set objective direction.
    if (model.objective_direction == ObjectiveDirection::Minimize) {
        cbc_model.setObjSense(1);
    } else {
        cbc_model.setObjSense(-1);
    }

    // Set variable types.
    std::vector<int> integer_variables;
    for (int variable_id = 0;
            variable_id < model.number_of_variables();
            ++variable_id) {
        if (model.variables_types[variable_id] != VariableType::Continuous)
            integer_variables.push_back(variable_id);
    }
    cbc_model.solver()->setInteger(integer_variables.data(), integer_variables.size());
}

void mathoptsolverscmake::reduce_printout(
        CbcModel& cbc_model)
{
    cbc_model.setLogLevel(0);
    cbc_model.messageHandler()->setLogLevel(0);
    cbc_model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
}

void mathoptsolverscmake::set_time_limit(
        CbcModel& cbc_model,
        double time_limit)
{
    cbc_model.setMaximumSeconds(time_limit);
}

void mathoptsolverscmake::solve(
        CbcModel& cbc_model)
{
    cbc_model.branchAndBound();
}

double mathoptsolverscmake::get_solution_value(
        const CbcModel& cbc_model)
{
    if (cbc_model.bestSolution() == nullptr)
        return cbc_model.getObjSense() * std::numeric_limits<double>::infinity();
    return cbc_model.getObjValue();
}

std::vector<double> mathoptsolverscmake::get_solution(
        const CbcModel& cbc_model)
{
    if (cbc_model.bestSolution() == nullptr)
        return {};
    return std::vector<double>(
            cbc_model.bestSolution(),
            cbc_model.bestSolution() + cbc_model.getNumCols());
}

double mathoptsolverscmake::get_bound(
        const CbcModel& cbc_model)
{
    return cbc_model.getBestPossibleObjValue();
}

int mathoptsolverscmake::get_number_of_nodes(
        const CbcModel& cbc_model)
{
    return cbc_model.getNodeCount();
}

#endif

#ifdef HIGHS_FOUND

void mathoptsolverscmake::load(
        Highs& highs_model,
        const MilpModel& model)
{
    HighsStatus highs_status;
    highs_status = highs_model.addCols(
            model.number_of_variables(),
            model.objective_coefficients.data(),
            model.variables_lower_bounds.data(),
            model.variables_upper_bounds.data(),
            0,
            NULL,
            NULL,
            NULL);
    highs_status = highs_model.addRows(
            model.number_of_constraints(),
            model.constraints_lower_bounds.data(),
            model.constraints_upper_bounds.data(),
            model.elements_coefficients.size(),
            model.constraints_starts.data(),
            model.elements_variables.data(),
            model.elements_coefficients.data());

    // Set objective sense.
    if (model.objective_direction == ObjectiveDirection::Minimize) {
        highs_model.changeObjectiveSense(ObjSense::kMinimize);
    } else {
        highs_model.changeObjectiveSense(ObjSense::kMaximize);
    }

    // Set variable types.
    std::vector<int> variables_indices(model.number_of_variables(), 0);
    std::iota(variables_indices.begin(), variables_indices.end(), 0);
    std::vector<HighsVarType> variables_types(model.number_of_variables(), HighsVarType::kContinuous);
    for (int variable_id = 0;
            variable_id < model.number_of_variables();
            ++variable_id) {
        switch (model.variables_types[variable_id]) {
        case VariableType::Continuous:
            variables_types[variable_id] = HighsVarType::kContinuous;
            break;
        case VariableType::Binary:
            variables_types[variable_id] = HighsVarType::kInteger;
            break;
        case VariableType::Integer:
            variables_types[variable_id] = HighsVarType::kInteger;
            break;
        }
    }
    highs_model.changeColsIntegrality(
            variables_indices.size(),
            variables_indices.data(),
            variables_types.data());
}

void mathoptsolverscmake::reduce_printout(
        Highs& highs_model)
{
    HighsStatus return_status = highs_model.setOptionValue(
            "log_to_console",
            false);
}

void mathoptsolverscmake::set_time_limit(
        Highs& highs_model,
        double time_limit)
{
    HighsStatus highs_status = highs_model.setOptionValue(
            "time_limit",
            time_limit);
}

void mathoptsolverscmake::set_log_file(
        Highs& highs_model,
        const std::string& log_file)
{
    HighsStatus return_status = highs_model.setOptionValue(
            "log_file",
            log_file);
}

void mathoptsolverscmake::solve(
        Highs& highs_model)
{
    HighsStatus return_status = highs_model.run();
}

std::vector<double> mathoptsolverscmake::get_solution(
        const Highs& highs_model)
{
    return highs_model.getSolution().col_value;
}

double mathoptsolverscmake::get_bound(
        const Highs& highs_model)
{
    return highs_model.getInfo().mip_dual_bound;
}

#endif

#ifdef XPRESS_FOUND

void mathoptsolverscmake::load(
        XPRSprob& xpress_model,
        const MilpModel& model)
{
    int status = 0;

    // Add variables.
    status = XPRSaddcols(
            xpress_model,
            model.number_of_variables(),
            0,  // ncoeffs
            model.objective_coefficients.data(),
            NULL,  // starts
            NULL,  // constraintind
            NULL,  // constraintcoef
            model.variables_lower_bounds.data(),
            model.variables_upper_bounds.data());
    if (status) {
        throw std::runtime_error(
                "mathoptsolverscmake::load: "
                "error loading variables; "
                "status: " + std::to_string(status) + ".");
    }

    // Add constraints.
    std::vector<char> constraints_types(model.number_of_constraints());
    std::vector<double> rhs(model.number_of_constraints());
    std::vector<double> rng(model.number_of_constraints());
    for (int constraint_id = 0;
            constraint_id < model.number_of_constraints();
            ++constraint_id) {
        if (model.constraints_lower_bounds[constraint_id] == -std::numeric_limits<double>::infinity()) {
            constraints_types[constraint_id] = 'L';
            rhs[constraint_id] = model.constraints_upper_bounds[constraint_id];
        } else if (model.constraints_upper_bounds[constraint_id] == std::numeric_limits<double>::infinity()) {
            constraints_types[constraint_id] = 'G';
            rhs[constraint_id] = model.constraints_lower_bounds[constraint_id];
        } else {
            constraints_types[constraint_id] = 'R';
            rhs[constraint_id] = model.constraints_upper_bounds[constraint_id];
            rng[constraint_id] = model.constraints_upper_bounds[constraint_id] - model.constraints_lower_bounds[constraint_id];
        }
    }
    status = XPRSaddrows(
            xpress_model,
            model.number_of_constraints(),
            model.number_of_elements(),
            constraints_types.data(),
            rhs.data(),
            rng.data(),
            model.constraints_starts.data(),
            model.elements_variables.data(),
            model.elements_coefficients.data());
    if (status) {
        throw std::runtime_error(
                "mathoptsolverscmake::load: "
                "error loading constraints; "
                "status: " + std::to_string(status) + ".");
    }

    // Set objective sense.
    if (model.objective_direction == ObjectiveDirection::Minimize) {
        XPRSchgobjsense(xpress_model, XPRS_OBJ_MINIMIZE);
    } else {
        XPRSchgobjsense(xpress_model, XPRS_OBJ_MAXIMIZE);
    }

    // Set variable types
    std::vector<int> variables_indices(model.number_of_variables(), 0);
    std::iota(variables_indices.begin(), variables_indices.end(), 0);
    std::vector<char> variables_types(model.number_of_variables(), 'A');
    for (int variable_id = 0;
            variable_id < model.number_of_variables();
            ++variable_id) {
        switch (model.variables_types[variable_id]) {
        case VariableType::Continuous:
            variables_types[variable_id] = 'B';
            break;
        case VariableType::Binary:
            variables_types[variable_id] = 'B';
            break;
        case VariableType::Integer:
            variables_types[variable_id] = 'I';
            break;
        }
    }
    status = XPRSchgcoltype(
            xpress_model,
            variables_indices.size(),
            variables_indices.data(),
            variables_types.data());
    if (status) {
        throw std::runtime_error(
                "mathoptsolverscmake::load: "
                "error setting variable types; "
                "status: " + std::to_string(status) + ".");
    }

}

void mathoptsolverscmake::set_time_limit(
        XPRSprob& xpress_model,
        double time_limit)
{
    XPRSsetdblcontrol(xpress_model, XPRS_TIMELIMIT, time_limit);
}

void mathoptsolverscmake::set_log_file(
        XPRSprob& xpress_model,
        const std::string& log_file)
{
    XPRSsetlogfile(xpress_model, log_file.c_str());
}

void mathoptsolverscmake::write_mps(
        const XPRSprob& xpress_model,
        const std::string& mps_file)
{
    XPRSwriteprob(xpress_model, mps_file.c_str(), "ot");
}

void mathoptsolverscmake::solve(
        XPRSprob& xpress_model)
{
    int status = XPRSmipoptimize(xpress_model, "");
    if (status) {
        throw std::runtime_error(
                "mathoptsolverscmake::solve: "
                "error solving model; "
                "status: " + std::to_string(status) + ".");
    }
}

double mathoptsolverscmake::get_solution_value(
        const XPRSprob& xpress_model)
{
    double objective_value = 0;
    XPRSgetdblattrib(xpress_model, XPRS_MIPOBJVAL, &objective_value);
    return objective_value;
}

std::vector<double> mathoptsolverscmake::get_solution(
        const XPRSprob& xpress_model)
{
    int number_of_variables = 0;
    XPRSgetintattrib(xpress_model, XPRS_ORIGINALCOLS, &number_of_variables);
    std::vector<double> solution(number_of_variables);
    int status = XPRSgetmipsol(
            xpress_model,
            solution.data(),
            NULL);
    if (status) {
        throw std::runtime_error(
                "mathoptsolverscmake::get_solution: "
                "error retrieving solution; "
                "status: " + std::to_string(status) + ".");
    }
    return solution;
}

double mathoptsolverscmake::get_bound(
        const XPRSprob& xpress_model)
{
    double bound = 0.0;
    XPRSgetdblattrib(xpress_model, XPRS_BESTBOUND, &bound);
    return bound;
}

int get_number_of_nodes(
        const XPRSprob& xpress_model)
{
    int number_of_nodes = 0;
    XPRSgetintattrib(xpress_model, XPRS_NODES, &number_of_nodes);
    return number_of_nodes;
}

#endif

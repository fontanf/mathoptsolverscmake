/**
 * Formatting the expressions of a model (linear, quadratic and expression
 * tree terms) in infix notation, and their polynomial degree.
 */

#include "mathoptsolverscmake/mathopt.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace mathoptsolverscmake;

namespace
{

const double inf = std::numeric_limits<double>::infinity();

/** Precedence levels: a sum, a product, a negation, a power, an atom. */
const int SUM = 1;
const int PRODUCT = 2;
const int NEGATION = 3;
const int POWER = 4;
const int ATOM = 5;

/** A number, with up to 15 significant digits. */
std::string format_number(double value)
{
    if (value == inf)
        return "inf";
    if (value == -inf)
        return "-inf";
    char text[32];
    std::snprintf(text, sizeof(text), "%.15g", value);
    std::string result(text);
    // snprintf follows the C locale's decimal separator.
    for (char& c: result)
        if (c == ',')
            c = '.';
    return result;
}

struct Formatted
{
    std::string text;
    int precedence;
};

/**
 * True if an expression starts with a minus sign that applies to the whole
 * expression: a negation, or a product whose first factor is negative ("-42
 * * x"), so that "a + e" can be written "a - e'".
 */
bool is_negative(const Formatted& formatted)
{
    return !formatted.text.empty()
        && formatted.text[0] == '-'
        && formatted.precedence >= PRODUCT;
}

std::string parenthesize(const Formatted& formatted, int minimum_precedence)
{
    if (formatted.precedence < minimum_precedence)
        return "(" + formatted.text + ")";
    return formatted.text;
}

/** One of the expression trees of a model (see mathopt.hpp). */
struct Tree
{
    const std::vector<char>& operators;
    const std::vector<double>& values;
    const std::vector<int>& variables;
    const std::vector<int>& number_of_children;
};

class ExpressionFormatter
{

public:

    ExpressionFormatter(
            const MathOptModel& model,
            bool expand_defined_variables):
        model_(model),
        expand_defined_variables_(expand_defined_variables),
        defined_variables_(model.number_of_defined_variables()),
        defined_variables_done_(model.number_of_defined_variables(), false)
    { }

    /** Formats the tree whose root is at 'node'; sets 'node' after its last node. */
    Formatted format(const Tree& tree, int& node)
    {
        char op = tree.operators[node];
        int number_of_children = tree.number_of_children[node];
        int current = node;
        ++node;
        std::vector<Formatted> children;
        for (int child = 0; child < number_of_children; ++child)
            children.push_back(format(tree, node));

        switch (op) {
        case 'k': {
            double value = tree.values[current];
            return {format_number(value), (value < 0)? NEGATION: ATOM};
        } case 'v': {
            return {model_.variable_name(tree.variables[current]), ATOM};
        } case 'd': {
            int defined_variable_id = tree.variables[current];
            if (!expand_defined_variables_)
                return {"d" + std::to_string(defined_variable_id), ATOM};
            return defined_variable(defined_variable_id);
        } case '+': {
            std::string text = children[0].text;
            for (size_t child = 1; child < children.size(); ++child) {
                if (is_negative(children[child]))
                    text += " - " + children[child].text.substr(1);
                else
                    text += " + " + children[child].text;
            }
            return {text, SUM};
        } case '-': {
            return {children[0].text + " - " + parenthesize(children[1], PRODUCT), SUM};
        } case '*': {
            std::string text = parenthesize(children[0], PRODUCT);
            for (size_t child = 1; child < children.size(); ++child)
                text += " * " + parenthesize(children[child], NEGATION);
            return {text, PRODUCT};
        } case '/': {
            return {parenthesize(children[0], PRODUCT) + " / " + parenthesize(children[1], NEGATION), PRODUCT};
        } case 'n': {
            return {"-" + parenthesize(children[0], PRODUCT), NEGATION};
        } case 'p': {
            return {parenthesize(children[0], ATOM) + " ^ " + parenthesize(children[1], ATOM), POWER};
        } case 'm':
          case 'M': {
            std::string text = (op == 'm')? "min(": "max(";
            for (size_t child = 0; child < children.size(); ++child)
                text += ((child > 0)? ", ": "") + children[child].text;
            return {text + ")", ATOM};
        } default: {
            std::string name = function_name(op);
            return {name + "(" + (children.empty()? std::string(): children[0].text) + ")", ATOM};
        }
        }
    }

    /** Formats a whole tree (an empty one is the constant 0). */
    Formatted format_tree(const Tree& tree, int start, int end)
    {
        if (start >= end)
            return {"0", ATOM};
        int node = start;
        return format(tree, node);
    }

    /** The expression of a defined variable. */
    Formatted defined_variable(int defined_variable_id)
    {
        if (!defined_variables_done_[defined_variable_id]) {
            Tree tree{
                model_.defined_variables_elements_operators,
                model_.defined_variables_elements_values,
                model_.defined_variables_elements_variables,
                model_.defined_variables_elements_number_of_children};
            defined_variables_[defined_variable_id] = format_tree(
                    tree,
                    model_.defined_variables_starts[defined_variable_id],
                    model_.defined_variable_end(defined_variable_id));
            defined_variables_done_[defined_variable_id] = true;
        }
        return defined_variables_[defined_variable_id];
    }

private:

    static std::string function_name(char op)
    {
        switch (op) {
        case 'e': return "exp";
        case 'l': return "log";
        case 'g': return "log10";
        case 'q': return "sqrt";
        case 's': return "sin";
        case 'c': return "cos";
        case 't': return "tan";
        case 'h': return "tanh";
        case 'a': return "abs";
        default: return std::string("op_") + op;
        }
    }

    const MathOptModel& model_;
    bool expand_defined_variables_;
    std::vector<Formatted> defined_variables_;
    std::vector<bool> defined_variables_done_;

};

/** Appends a term "coefficient * factors" to a sum, with its sign. */
void append_term(
        std::string& text,
        double coefficient,
        const std::string& factors)
{
    if (coefficient == 0.0)
        return;
    double magnitude = std::abs(coefficient);
    std::string term = (magnitude == 1.0)? factors: format_number(magnitude) + " * " + factors;
    if (text.empty())
        text = (coefficient < 0)? "-" + term: term;
    else
        text += ((coefficient < 0)? " - ": " + ") + term;
}

/** Appends an expression (a tree, a black-box function...) to a sum. */
void append_expression(
        std::string& text,
        const Formatted& expression)
{
    if (text.empty()) {
        text = expression.text;
    } else if (is_negative(expression)) {
        text += " - " + expression.text.substr(1);
    } else {
        text += " + " + expression.text;
    }
}

std::string quadratic_factors(const MathOptModel& model, int variable_1_id, int variable_2_id)
{
    if (variable_1_id == variable_2_id)
        return model.variable_name(variable_1_id) + " ^ 2";
    return model.variable_name(variable_1_id) + " * " + model.variable_name(variable_2_id);
}

/*
 * Degrees.
 */

class DegreeComputer
{

public:

    explicit DegreeComputer(const MathOptModel& model):
        model_(model),
        defined_variables_(model.number_of_defined_variables(), -1)
    { }

    /** The degree of the tree whose root is at 'node'; sets 'node' after its last node. */
    double degree(const Tree& tree, int& node)
    {
        char op = tree.operators[node];
        int number_of_children = tree.number_of_children[node];
        int current = node;
        ++node;
        std::vector<double> children;
        std::vector<int> children_roots;
        for (int child = 0; child < number_of_children; ++child) {
            children_roots.push_back(node);
            children.push_back(degree(tree, node));
        }

        switch (op) {
        case 'k': return 0;
        case 'v': return 1;
        case 'd': return defined_variable(tree.variables[current]);
        case '+':
        case '-': return *std::max_element(children.begin(), children.end());
        case 'n': return children[0];
        case '*': {
            double sum = 0;
            for (double child: children)
                sum += child;
            return sum;
        } case '/': {
            return (children[1] == 0)? children[0]: inf;
        } case 'p': {
            if (children[0] == 0 && children[1] == 0)
                return 0;
            // A polynomial raised to a constant non-negative integer power.
            int exponent_root = children_roots[1];
            if (tree.operators[exponent_root] == 'k') {
                double exponent = tree.values[exponent_root];
                if (exponent >= 0 && exponent == std::floor(exponent))
                    return children[0] * exponent;
            }
            return inf;
        } default: {
            // Functions, min and max: constant only if their arguments are.
            for (double child: children)
                if (child != 0)
                    return inf;
            return 0;
        }
        }
    }

    double degree_tree(const Tree& tree, int start, int end)
    {
        if (start >= end)
            return 0;
        int node = start;
        return degree(tree, node);
    }

    double defined_variable(int defined_variable_id)
    {
        if (defined_variables_[defined_variable_id] < 0) {
            Tree tree{
                model_.defined_variables_elements_operators,
                model_.defined_variables_elements_values,
                model_.defined_variables_elements_variables,
                model_.defined_variables_elements_number_of_children};
            defined_variables_[defined_variable_id] = degree_tree(
                    tree,
                    model_.defined_variables_starts[defined_variable_id],
                    model_.defined_variable_end(defined_variable_id));
        }
        return defined_variables_[defined_variable_id];
    }

private:

    const MathOptModel& model_;
    std::vector<double> defined_variables_;

};

void check_constraint_id(const MathOptModel& model, int constraint_id)
{
    if (constraint_id < 0 || constraint_id >= model.number_of_constraints()) {
        throw std::invalid_argument(
                "constraint " + std::to_string(constraint_id) + " doesn't exist.");
    }
}

}

std::string MathOptModel::format_constraint_expression(
        int constraint_id,
        bool expand_defined_variables) const
{
    check_constraint_id(*this, constraint_id);
    std::string text;
    for (int element_id = this->constraints_starts[constraint_id];
            element_id < this->constraint_end(constraint_id);
            ++element_id) {
        append_term(
                text,
                this->elements_coefficients[element_id],
                this->variable_name(this->elements_variables[element_id]));
    }
    if (!this->quadratic_elements_constraints_starts.empty()) {
        for (int element_id = this->quadratic_elements_constraints_starts[constraint_id];
                element_id < this->quadratic_constraint_end(constraint_id);
                ++element_id) {
            append_term(
                    text,
                    this->quadratic_elements_coefficients[element_id],
                    quadratic_factors(
                        *this,
                        this->quadratic_elements_variables_1[element_id],
                        this->quadratic_elements_variables_2[element_id]));
        }
    }
    if (!this->nonlinear_elements_constraints_starts.empty()) {
        int start = this->nonlinear_elements_constraints_starts[constraint_id];
        int end = this->nonlinear_constraint_end(constraint_id);
        if (start < end) {
            ExpressionFormatter formatter(*this, expand_defined_variables);
            Tree tree{
                this->nonlinear_elements_operators,
                this->nonlinear_elements_values,
                this->nonlinear_elements_variables,
                this->nonlinear_elements_number_of_children};
            append_expression(text, formatter.format_tree(tree, start, end));
        }
    }
    if (!this->constraints_functions.empty() && bool(this->constraints_functions[constraint_id]))
        append_expression(text, {"g" + std::to_string(constraint_id) + "(x)", ATOM});
    return text.empty()? "0": text;
}

std::string MathOptModel::format_objective_expression(
        bool expand_defined_variables) const
{
    std::string text;
    for (int variable_id = 0; variable_id < (int)this->objective_coefficients.size(); ++variable_id)
        append_term(text, this->objective_coefficients[variable_id], this->variable_name(variable_id));
    for (int element_id = 0;
            element_id < (int)this->objective_quadratic_elements_coefficients.size();
            ++element_id) {
        append_term(
                text,
                this->objective_quadratic_elements_coefficients[element_id],
                quadratic_factors(
                    *this,
                    this->objective_quadratic_elements_variables_1[element_id],
                    this->objective_quadratic_elements_variables_2[element_id]));
    }
    if (!this->objective_nonlinear_elements_operators.empty()) {
        ExpressionFormatter formatter(*this, expand_defined_variables);
        Tree tree{
            this->objective_nonlinear_elements_operators,
            this->objective_nonlinear_elements_values,
            this->objective_nonlinear_elements_variables,
            this->objective_nonlinear_elements_number_of_children};
        append_expression(
                text,
                formatter.format_tree(tree, 0, (int)this->objective_nonlinear_elements_operators.size()));
    }
    if (bool(this->objective_function))
        append_expression(text, {"f(x)", ATOM});
    return text.empty()? "0": text;
}

std::string MathOptModel::format_defined_variable_expression(
        int defined_variable_id,
        bool expand_defined_variables) const
{
    if (defined_variable_id < 0 || defined_variable_id >= this->number_of_defined_variables()) {
        throw std::invalid_argument(
                "defined variable " + std::to_string(defined_variable_id) + " doesn't exist.");
    }
    ExpressionFormatter formatter(*this, expand_defined_variables);
    Tree tree{
        this->defined_variables_elements_operators,
        this->defined_variables_elements_values,
        this->defined_variables_elements_variables,
        this->defined_variables_elements_number_of_children};
    return formatter.format_tree(
            tree,
            this->defined_variables_starts[defined_variable_id],
            this->defined_variable_end(defined_variable_id)).text;
}

double MathOptModel::constraint_degree(int constraint_id) const
{
    check_constraint_id(*this, constraint_id);
    if (!this->constraints_functions.empty() && bool(this->constraints_functions[constraint_id]))
        return inf;
    double degree = 0;
    for (int element_id = this->constraints_starts[constraint_id];
            element_id < this->constraint_end(constraint_id);
            ++element_id) {
        if (this->elements_coefficients[element_id] != 0.0)
            degree = 1;
    }
    if (!this->quadratic_elements_constraints_starts.empty()) {
        for (int element_id = this->quadratic_elements_constraints_starts[constraint_id];
                element_id < this->quadratic_constraint_end(constraint_id);
                ++element_id) {
            if (this->quadratic_elements_coefficients[element_id] != 0.0)
                degree = 2;
        }
    }
    if (!this->nonlinear_elements_constraints_starts.empty()) {
        DegreeComputer computer(*this);
        Tree tree{
            this->nonlinear_elements_operators,
            this->nonlinear_elements_values,
            this->nonlinear_elements_variables,
            this->nonlinear_elements_number_of_children};
        degree = (std::max)(degree, computer.degree_tree(
                    tree,
                    this->nonlinear_elements_constraints_starts[constraint_id],
                    this->nonlinear_constraint_end(constraint_id)));
    }
    return degree;
}

double MathOptModel::objective_degree() const
{
    if (bool(this->objective_function))
        return inf;
    double degree = 0;
    for (double coefficient: this->objective_coefficients)
        if (coefficient != 0.0)
            degree = 1;
    for (double coefficient: this->objective_quadratic_elements_coefficients)
        if (coefficient != 0.0)
            degree = 2;
    if (!this->objective_nonlinear_elements_operators.empty()) {
        DegreeComputer computer(*this);
        Tree tree{
            this->objective_nonlinear_elements_operators,
            this->objective_nonlinear_elements_values,
            this->objective_nonlinear_elements_variables,
            this->objective_nonlinear_elements_number_of_children};
        degree = (std::max)(degree, computer.degree_tree(
                    tree, 0, (int)this->objective_nonlinear_elements_operators.size()));
    }
    return degree;
}

std::ostream& MathOptModel::format_constraint(
        std::ostream& os,
        int constraint_id) const
{
    std::string expression = this->format_constraint_expression(constraint_id);
    double lower = this->constraints_lower_bounds[constraint_id];
    double upper = this->constraints_upper_bounds[constraint_id];
    os << this->constraint_name(constraint_id) << ": ";
    if (lower == upper) {
        os << expression << " = " << format_number(lower);
    } else if (lower != -inf && upper != inf) {
        os << format_number(lower) << " <= " << expression << " <= " << format_number(upper);
    } else if (upper != inf) {
        os << expression << " <= " << format_number(upper);
    } else if (lower != -inf) {
        os << expression << " >= " << format_number(lower);
    } else {
        os << expression;
    }
    os << std::endl;
    return os;
}

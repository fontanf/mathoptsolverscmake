#include "mathoptsolverscmake/mathopt_nonlinear.hpp"

#include <utility>

using namespace mathoptsolverscmake;

namespace
{

NonlinearExpression leaf(char op, double value, int variable_id)
{
    return [op, value, variable_id](
            std::vector<char>& operators,
            std::vector<double>& values,
            std::vector<int>& variables,
            std::vector<int>& parents,
            int parent_id)
    {
        operators.push_back(op);
        values.push_back(value);
        variables.push_back(variable_id);
        parents.push_back(parent_id);
    };
}

NonlinearExpression unary(char op, NonlinearExpression operand)
{
    return [op, operand](
            std::vector<char>& operators,
            std::vector<double>& values,
            std::vector<int>& variables,
            std::vector<int>& parents,
            int parent_id)
    {
        int node_id = (int)operators.size();
        operators.push_back(op);
        values.push_back(0.0);
        variables.push_back(-1);
        parents.push_back(parent_id);
        operand(operators, values, variables, parents, node_id);
    };
}

NonlinearExpression nary(char op, std::vector<NonlinearExpression> operands)
{
    return [op, operands](
            std::vector<char>& operators,
            std::vector<double>& values,
            std::vector<int>& variables,
            std::vector<int>& parents,
            int parent_id)
    {
        int node_id = (int)operators.size();
        operators.push_back(op);
        values.push_back(0.0);
        variables.push_back(-1);
        parents.push_back(parent_id);
        for (const NonlinearExpression& operand: operands)
            operand(operators, values, variables, parents, node_id);
    };
}

}

NonlinearExpression mathoptsolverscmake::nonlinear_variable(int variable_id)
{
    return leaf('v', 0.0, variable_id);
}

NonlinearExpression mathoptsolverscmake::nonlinear_constant(double value)
{
    return leaf('k', value, -1);
}

NonlinearExpression mathoptsolverscmake::nonlinear_negate(NonlinearExpression operand)
{
    return unary('n', std::move(operand));
}

NonlinearExpression mathoptsolverscmake::nonlinear_exp(NonlinearExpression operand)
{
    return unary('e', std::move(operand));
}

NonlinearExpression mathoptsolverscmake::nonlinear_log(NonlinearExpression operand)
{
    return unary('l', std::move(operand));
}

NonlinearExpression mathoptsolverscmake::nonlinear_sqrt(NonlinearExpression operand)
{
    return unary('q', std::move(operand));
}

NonlinearExpression mathoptsolverscmake::nonlinear_sin(NonlinearExpression operand)
{
    return unary('s', std::move(operand));
}

NonlinearExpression mathoptsolverscmake::nonlinear_cos(NonlinearExpression operand)
{
    return unary('c', std::move(operand));
}

NonlinearExpression mathoptsolverscmake::nonlinear_tan(NonlinearExpression operand)
{
    return unary('t', std::move(operand));
}

NonlinearExpression mathoptsolverscmake::nonlinear_sum(std::vector<NonlinearExpression> terms)
{
    return nary('+', std::move(terms));
}

NonlinearExpression mathoptsolverscmake::nonlinear_product(std::vector<NonlinearExpression> factors)
{
    return nary('*', std::move(factors));
}

NonlinearExpression mathoptsolverscmake::nonlinear_subtract(
        NonlinearExpression left,
        NonlinearExpression right)
{
    return nary('-', {std::move(left), std::move(right)});
}

NonlinearExpression mathoptsolverscmake::nonlinear_divide(
        NonlinearExpression left,
        NonlinearExpression right)
{
    return nary('/', {std::move(left), std::move(right)});
}

NonlinearExpression mathoptsolverscmake::nonlinear_power(
        NonlinearExpression base,
        NonlinearExpression exponent)
{
    return nary('p', {std::move(base), std::move(exponent)});
}

void mathoptsolverscmake::set_objective_nonlinear_expression(
        MathOptModel& model,
        const NonlinearExpression& expression)
{
    expression(
            model.objective_nonlinear_elements_operators,
            model.objective_nonlinear_elements_values,
            model.objective_nonlinear_elements_variables,
            model.objective_nonlinear_elements_parent,
            -1);
}

void mathoptsolverscmake::add_constraint_nonlinear_expression(
        MathOptModel& model,
        int constraint_id,
        const NonlinearExpression& expression)
{
    if ((int)model.nonlinear_elements_constraints_starts.size() < model.number_of_constraints())
        model.nonlinear_elements_constraints_starts.resize(model.number_of_constraints(), 0);
    model.nonlinear_elements_constraints_starts[constraint_id] =
        (int)model.nonlinear_elements_operators.size();
    expression(
            model.nonlinear_elements_operators,
            model.nonlinear_elements_values,
            model.nonlinear_elements_variables,
            model.nonlinear_elements_parent,
            -1);
}

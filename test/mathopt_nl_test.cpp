#include "mathoptsolverscmake/mathopt_nl.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace mathoptsolverscmake;

namespace
{

const double inf = std::numeric_limits<double>::infinity();

MathOptModel read(const std::string& nl)
{
    std::istringstream in(nl);
    return read_nl(in);
}

std::string write(const MathOptModel& model)
{
    std::ostringstream out;
    write_nl(model, out);
    return out.str();
}

/**
 * Builds the body of a binary ("b") .nl file: one ASCII byte per segment key
 * or expression node kind, followed by its values as raw binary numbers.
 */
class BinaryNlBody
{

public:

    explicit BinaryNlBody(bool big_endian): big_endian_(big_endian) { }

    BinaryNlBody& key(char c) { body_ += c; return *this; }
    BinaryNlBody& key(char c, int32_t value) { body_ += c; return put(value); }
    BinaryNlBody& integer(int32_t value) { return put(value); }
    BinaryNlBody& real(double value) { return put(value); }
    BinaryNlBody& constant(double value) { body_ += 'n'; return put(value); }
    BinaryNlBody& short_constant(int16_t value) { body_ += 's'; return put(value); }
    BinaryNlBody& bound(char type) { body_ += type; return *this; }
    BinaryNlBody& string(const std::string& value) { put((int32_t)value.size()); body_ += value; return *this; }

    const std::string& str() const { return body_; }

private:

    template <typename T>
    BinaryNlBody& put(T value)
    {
        char bytes[sizeof(T)];
        std::memcpy(bytes, &value, sizeof(T));
        uint16_t probe = 1;
        bool host_little_endian = (*reinterpret_cast<unsigned char*>(&probe) == 1);
        if (big_endian_ == host_little_endian)
            std::reverse(bytes, bytes + sizeof(T));
        body_.append(bytes, sizeof(T));
        return *this;
    }

    bool big_endian_;
    std::string body_;

};

/** The binary encoding of the LinearBoundTypesAndSenses model below. */
std::string binary_linear_model(bool big_endian)
{
    std::string header =
            "b0\n"
            "2 2 1 0 1\n"
            "0 0\n"
            "0 0\n"
            "0 0 0\n"
            "0 0 " + std::string(big_endian? "2": "1") + " 0\n"
            "0 0 0 0 0\n"
            "4 2\n"
            "0 0\n"
            "0 0 0 0 0\n";
    BinaryNlBody body(big_endian);
    body.key('b').bound('0').real(0).real(10).bound('3');
    body.key('r').bound('1').real(5).bound('4').real(0);
    body.key('C', 0).constant(0);
    body.key('C', 1).constant(0);
    body.key('O', 0).integer(0).constant(0);
    body.key('J', 0).integer(2).integer(0).real(1).integer(1).real(1);
    body.key('J', 1).integer(2).integer(0).real(1).integer(1).real(-1);
    body.key('G', 0).integer(2).integer(0).real(1).integer(1).real(2);
    return header + body.str();
}

void expect_linear_model(const MathOptModel& model)
{
    EXPECT_EQ(model.objective_direction, ObjectiveDirection::Minimize);
    ASSERT_EQ(model.number_of_variables(), 2);
    ASSERT_EQ(model.number_of_constraints(), 2);
    EXPECT_EQ(model.variables_lower_bounds[0], 0.0);
    EXPECT_EQ(model.variables_upper_bounds[0], 10.0);
    EXPECT_EQ(model.variables_lower_bounds[1], -inf);
    EXPECT_EQ(model.variables_upper_bounds[1], inf);
    EXPECT_EQ(model.constraints_lower_bounds[0], -inf);
    EXPECT_EQ(model.constraints_upper_bounds[0], 5.0);
    EXPECT_EQ(model.constraints_lower_bounds[1], 0.0);
    EXPECT_EQ(model.constraints_upper_bounds[1], 0.0);
    EXPECT_EQ(model.objective_coefficients[0], 1.0);
    EXPECT_EQ(model.objective_coefficients[1], 2.0);
    ASSERT_EQ(model.number_of_elements(), 4);
    EXPECT_EQ(model.elements_coefficients[3], -1.0);
}

}

TEST(MathOptNl, LinearBoundTypesAndSenses)
{
    // minimize x0 + 2 x1
    // s.t. c0: x0 + x1 <= 5       (bound type 1: upper only)
    //      c1: x0 - x1 = 0        (bound type 4: equality)
    // x0 in [0, 10] (type 0: range), x1 free (type 3).
    MathOptModel model = read(
            "g0\n"
            "2 2 1 0 1\n"
            "0 0\n"
            "0 0\n"
            "0 0 0\n"
            "0 0\n"
            "0 0 0 0 0\n"
            "4 2\n"
            "0 0\n"
            "0 0 0 0 0\n"
            "b\n"
            "0 0 10\n"
            "3\n"
            "r\n"
            "1 5\n"
            "4 0\n"
            "C0\n"
            "n0\n"
            "C1\n"
            "n0\n"
            "O0 0\n"
            "n0\n"
            "J0 2\n"
            "0 1\n"
            "1 1\n"
            "J1 2\n"
            "0 1\n"
            "1 -1\n"
            "G0 2\n"
            "0 1\n"
            "1 2\n");

    EXPECT_EQ(model.objective_direction, ObjectiveDirection::Minimize);
    ASSERT_EQ(model.number_of_variables(), 2);
    ASSERT_EQ(model.number_of_constraints(), 2);

    EXPECT_EQ(model.variables_lower_bounds[0], 0.0);
    EXPECT_EQ(model.variables_upper_bounds[0], 10.0);
    EXPECT_EQ(model.variables_lower_bounds[1], -inf);
    EXPECT_EQ(model.variables_upper_bounds[1], inf);

    EXPECT_EQ(model.objective_coefficients[0], 1.0);
    EXPECT_EQ(model.objective_coefficients[1], 2.0);

    EXPECT_EQ(model.constraints_lower_bounds[0], -inf);
    EXPECT_EQ(model.constraints_upper_bounds[0], 5.0);
    EXPECT_EQ(model.constraints_lower_bounds[1], 0.0);
    EXPECT_EQ(model.constraints_upper_bounds[1], 0.0);

    std::vector<double> x = {1.0, 2.0};
    EXPECT_EQ(model.evaluate_objective(x), 1.0 + 2.0 * 2.0);
    EXPECT_EQ(model.evaluate_constraint(x, 0), 3.0);
    EXPECT_EQ(model.evaluate_constraint(x, 1), -1.0);
}

TEST(MathOptNl, RangeConstraintAndFixedVariable)
{
    // c0: 2 <= x0 + x1 <= 8 (bound type 0), x0 fixed at 3 (bound type 4).
    MathOptModel model = read(
            "g0\n"
            "2 1 1 0 0\n"
            "0 0\n"
            "0 0\n"
            "0 0 0\n"
            "0 0\n"
            "0 0 0 0 0\n"
            "2 0\n"
            "0 0\n"
            "0 0 0 0 0\n"
            "b\n"
            "4 3\n"
            "2 0\n"
            "r\n"
            "0 2 8\n"
            "C0\n"
            "n0\n"
            "O0 0\n"
            "n0\n"
            "J0 2\n"
            "0 1\n"
            "1 1\n");

    ASSERT_EQ(model.number_of_variables(), 2);
    EXPECT_EQ(model.variables_lower_bounds[0], 3.0);
    EXPECT_EQ(model.variables_upper_bounds[0], 3.0);
    EXPECT_EQ(model.variables_lower_bounds[1], 0.0);
    EXPECT_EQ(model.variables_upper_bounds[1], inf);
    EXPECT_EQ(model.constraints_lower_bounds[0], 2.0);
    EXPECT_EQ(model.constraints_upper_bounds[0], 8.0);
}

TEST(MathOptNl, NonlinearObjective)
{
    // minimize exp(x0) + x1
    // s.t. c0: x0 + x1 <= 5
    // x0 in [0, 10], x1 in [0, inf).
    MathOptModel model = read(
            "g0\n"
            "2 1 1 0 0\n"
            "0 1\n"
            "0 0\n"
            "0 1 0\n"
            "0 0\n"
            "0 0 0 0 0\n"
            "2 1\n"
            "0 0\n"
            "0 0 0 0 0\n"
            "b\n"
            "0 0 10\n"
            "2 0\n"
            "r\n"
            "1 5\n"
            "C0\n"
            "n0\n"
            "O0 0\n"
            "o44\n"
            "v0\n"
            "J0 2\n"
            "0 1\n"
            "1 1\n"
            "G0 1\n"
            "1 1\n");

    ASSERT_EQ(model.number_of_variables(), 2);
    EXPECT_TRUE(model.has_nonlinear());
    EXPECT_FALSE(model.has_quadratic());

    std::vector<double> x = {0.0, 3.0};
    EXPECT_NEAR(model.evaluate_objective(x), std::exp(0.0) + 3.0, 1e-12);
    x = {2.0, 3.0};
    EXPECT_NEAR(model.evaluate_objective(x), std::exp(2.0) + 3.0, 1e-12);
}

TEST(MathOptNl, LinearBinaryFromHeaderNbv)
{
    // A variable declared "binary" via the header's nbv count, but whose
    // own 'b' segment bound is [0, +infinity): read_nl() must trust the
    // header's explicit nbv/niv split for the linear portion, not infer
    // Binary-vs-Integer from bounds (that heuristic only applies to the
    // nonlinear-variable positions, which have no such explicit signal).
    MathOptModel model = read(
            "g0\n"
            "1 0 1 0 0\n"
            "0 0\n"
            "0 0\n"
            "0 0 0\n"
            "0 0\n"
            "1 0 0 0 0\n"
            "0 1\n"
            "0 0\n"
            "0 0 0 0 0\n"
            "b\n"
            "2 0\n"
            "r\n"
            "O0 0\n"
            "n0\n"
            "G0 1\n"
            "0 1\n");

    ASSERT_EQ(model.number_of_variables(), 1);
    EXPECT_EQ(model.variables_types[0], VariableType::Binary);
    EXPECT_EQ(model.variables_lower_bounds[0], 0.0);
    EXPECT_EQ(model.variables_upper_bounds[0], inf);
}

TEST(MathOptNl, LinearContinuousBeforeBinaryAndInteger)
{
    // With no nonlinear variables, .nl orders the linear ones as: other
    // continuous, then the nbv binary, then the niv integer variables --
    // so with nbv = niv = 1, x0 and x1 are continuous, x2 binary, x3 integer.
    MathOptModel model = read(
            "g0\n"
            "4 0 1 0 0\n"
            "0 0\n"
            "0 0\n"
            "0 0 0\n"
            "0 0\n"
            "1 1 0 0 0\n"
            "0 4\n"
            "0 0\n"
            "0 0 0 0 0\n"
            "b\n"
            "0 2.5 7\n"
            "2 3\n"
            "0 0 1\n"
            "0 0 9\n"
            "r\n"
            "O0 0\n"
            "n0\n"
            "G0 4\n"
            "0 1\n"
            "1 1\n"
            "2 1\n"
            "3 1\n");

    ASSERT_EQ(model.number_of_variables(), 4);
    EXPECT_EQ(model.variables_types[0], VariableType::Continuous);
    EXPECT_EQ(model.variables_types[1], VariableType::Continuous);
    EXPECT_EQ(model.variables_types[2], VariableType::Binary);
    EXPECT_EQ(model.variables_types[3], VariableType::Integer);
}

TEST(MathOptNl, BinaryFormatLittleEndian)
{
    expect_linear_model(read(binary_linear_model(false)));
}

TEST(MathOptNl, BinaryFormatBigEndian)
{
    expect_linear_model(read(binary_linear_model(true)));
}

TEST(MathOptNl, BinaryFormatNonlinearAndSuffix)
{
    // minimize x0^2 + x0 + x1 (o54 sum of a power, with the exponent as a
    // binary-only 's' short constant, and the linear part from 'G'), with
    // an integer-valued suffix and an initial point.
    BinaryNlBody body(false);
    body.key('S', 0).integer(1).string("sosno").integer(0).integer(1);
    body.key('b').bound('3').bound('3');
    body.key('x', 1).integer(1).real(0.5);
    body.key('O', 0).integer(0).key('o', 54).integer(2)
        .key('o', 5).key('v', 0).short_constant(2)
        .constant(0);
    body.key('G', 0).integer(2).integer(0).real(1).integer(1).real(1);
    MathOptModel model = read(
            "b0\n"
            "2 0 1 0 0\n"
            "0 1\n"
            "0 0\n"
            "0 1 0\n"
            "0 0 1 0\n"
            "0 0 0 0 0\n"
            "0 2\n"
            "0 0\n"
            "0 0 0 0 0\n"
            + body.str());

    ASSERT_EQ(model.number_of_variables(), 2);
    EXPECT_EQ(model.objective_coefficients[0], 1.0);
    EXPECT_EQ(model.objective_coefficients[1], 1.0);
    ASSERT_EQ(model.variables_initial_values.size(), 2u);
    EXPECT_EQ(model.variables_initial_values[1], 0.5);
    // Objective tree: + (pow(x0, 2), 0), in DFS pre-order.
    std::string operators(
            model.objective_nonlinear_elements_operators.begin(),
            model.objective_nonlinear_elements_operators.end());
    EXPECT_EQ(operators, "+pvkk");
    EXPECT_EQ(model.objective_nonlinear_elements_values[3], 2.0);
    EXPECT_DOUBLE_EQ(model.evaluate_objective({3.0, 1.0}), 9.0 + 3.0 + 1.0);
}

TEST(MathOptNl, BinaryFormatTruncatedThrows)
{
    std::string full = binary_linear_model(false);
    EXPECT_THROW(read(full.substr(0, full.size() - 3)), std::invalid_argument);
}

TEST(MathOptNl, ObsoleteHFormatThrows)
{
    EXPECT_THROW(read("h0\n"), std::invalid_argument);
}

namespace
{

/**
 * Two common subexpressions: V2 = 2 x0 + x0^2 (a linear and a nonlinear
 * part) and V3 = V2 * x1 (referencing V2). minimize V2 s.t. c0: V3 <= 100.
 */
const char* kCommonSubexpressionsNl =
        "g0\n"
        "2 1 1 0 0\n"
        "1 1\n"
        "0 0\n"
        "2 1 1\n"
        "0 0\n"
        "0 0 0 0 0\n"
        "0 0\n"
        "0 0\n"
        "2 0 0 0 0\n"
        "b\n"
        "3\n"
        "3\n"
        "r\n"
        "1 100\n"
        "V2 1 0\n"
        "0 2\n"
        "o5\n"
        "v0\n"
        "n2\n"
        "V3 0 0\n"
        "o2\n"
        "v2\n"
        "v1\n"
        "C0\n"
        "v3\n"
        "O0 0\n"
        "v2\n";

}

TEST(MathOptNl, CommonSubexpressionsAsDefinedVariables)
{
    MathOptModel model = read(kCommonSubexpressionsNl);
    ASSERT_EQ(model.number_of_defined_variables(), 2);
    EXPECT_TRUE(model.check());
    // V3 references V2, which is defined variable 0.
    std::string operators(
            model.defined_variables_elements_operators.begin(),
            model.defined_variables_elements_operators.end());
    EXPECT_EQ(operators, "+*kvpvk" "*dv");
    EXPECT_EQ(model.nonlinear_elements_operators, std::vector<char>{'d'});
    EXPECT_EQ(model.nonlinear_elements_variables, std::vector<int>{1});

    std::vector<double> x = {3.0, 2.0};
    EXPECT_DOUBLE_EQ(model.evaluate_objective(x), 2 * 3 + 9);
    EXPECT_DOUBLE_EQ(model.evaluate_constraint(x, 0), (2 * 3 + 9) * 2);
}

TEST(MathOptNl, ExpandDefinedVariables)
{
    MathOptModel model = read(kCommonSubexpressionsNl);
    MathOptModel expanded = expand_defined_variables(model);
    EXPECT_FALSE(expanded.has_defined_variables());
    EXPECT_TRUE(expanded.check());
    for (std::vector<double> x: {std::vector<double>{3.0, 2.0}, std::vector<double>{-1.5, 0.5}}) {
        EXPECT_DOUBLE_EQ(expanded.evaluate_objective(x), model.evaluate_objective(x));
        EXPECT_DOUBLE_EQ(expanded.evaluate_constraint(x, 0), model.evaluate_constraint(x, 0));
    }
}

TEST(MathOptNl, WriteThenReadDefinedVariables)
{
    // Written expanded, so it comes back without defined variables but
    // with the same values.
    MathOptModel model = read(kCommonSubexpressionsNl);
    MathOptModel again = read(write(model));
    EXPECT_FALSE(again.has_defined_variables());
    std::vector<double> x = {3.0, 2.0};
    EXPECT_DOUBLE_EQ(again.evaluate_objective(x), model.evaluate_objective(x));
    EXPECT_DOUBLE_EQ(again.evaluate_constraint(x, 0), model.evaluate_constraint(x, 0));
}

TEST(MathOptNl, UndefinedCommonSubexpressionThrows)
{
    // The objective references V2, which is never defined.
    EXPECT_THROW(
            read(
                    "g0\n"
                    "2 0 1 0 0\n"
                    "0 1\n"
                    "0 0\n"
                    "0 1 0\n"
                    "0 0\n"
                    "0 0 0 0 0\n"
                    "0 0\n"
                    "0 0\n"
                    "0 0 1 0 0\n"
                    "b\n"
                    "3\n"
                    "3\n"
                    "O0 0\n"
                    "v2\n"),
            std::invalid_argument);
}

TEST(MathOptNl, AbsTanhLog10MinMax)
{
    // minimize abs(x0) + tanh(x1) + log10(x1) + min(x0, x1, 1) + max(x0, x1)
    MathOptModel model = read(
            "g0\n"
            "2 0 1 0 0\n"
            "0 1\n"
            "0 0\n"
            "0 2 0\n"
            "0 0\n"
            "0 0 0 0 0\n"
            "0 0\n"
            "0 0\n"
            "0 0 0 0 0\n"
            "b\n"
            "3\n"
            "3\n"
            "O0 0\n"
            "o54\n"
            "5\n"
            "o15\n"
            "v0\n"
            "o37\n"
            "v1\n"
            "o42\n"
            "v1\n"
            "o11\n"
            "3\n"
            "v0\n"
            "v1\n"
            "n1\n"
            "o12\n"
            "2\n"
            "v0\n"
            "v1\n");
    std::string operators(
            model.objective_nonlinear_elements_operators.begin(),
            model.objective_nonlinear_elements_operators.end());
    EXPECT_EQ(operators, "+avhvgvmvvkMvv");

    auto expected = [](double x0, double x1) {
        return std::abs(x0) + std::tanh(x1) + std::log10(x1)
            + (std::min)({x0, x1, 1.0}) + (std::max)(x0, x1);
    };
    EXPECT_DOUBLE_EQ(model.evaluate_objective({-2.0, 5.0}), expected(-2.0, 5.0));
    EXPECT_DOUBLE_EQ(model.evaluate_objective({3.0, 0.5}), expected(3.0, 0.5));

    // And they are written back with the same opcodes.
    MathOptModel again = read(write(model));
    EXPECT_DOUBLE_EQ(again.evaluate_objective({-2.0, 5.0}), expected(-2.0, 5.0));
}

TEST(MathOptNl, UnsupportedOpcodeThrows)
{
    // o13 is FLOOR, not one of MathOptModel's supported operators.
    EXPECT_THROW(
            read(
                    "g0\n"
                    "1 0 1 0 0\n"
                    "0 1\n"
                    "0 0\n"
                    "0 1 0\n"
                    "0 0\n"
                    "0 0 0 0 0\n"
                    "0 0\n"
                    "0 0\n"
                    "0 0 0 0 0\n"
                    "b\n"
                    "3\n"
                    "r\n"
                    "O0 0\n"
                    "o13\n"
                    "v0\n"),
            std::invalid_argument);
}

TEST(MathOptNl, MissingFileThrows)
{
    EXPECT_THROW(read_nl("/no/such/file.nl"), std::invalid_argument);
}

TEST(MathOptNl, WriteThenReadLinearRoundTrip)
{
    // All-continuous linear model: write_nl()'s variable partition is the
    // identity in this case (no nonlinear content to sort to the front),
    // so an exact field-by-field comparison is meaningful, not just a
    // semantic (evaluate_*) one.
    MathOptModel original(2, 2, 4);
    original.objective_direction = ObjectiveDirection::Maximize;
    original.variables_lower_bounds = {0.0, -5.0};
    original.variables_upper_bounds = {10.0, 5.0};
    original.objective_coefficients = {3.0, 2.0};
    original.constraints_starts = {0, 2};
    original.elements_variables = {0, 1, 0, 1};
    original.elements_coefficients = {1.0, 1.0, 1.0, -1.0};
    original.constraints_lower_bounds = {-inf, 0.0};
    original.constraints_upper_bounds = {4.0, 0.0};

    MathOptModel round_tripped = read(write(original));

    ASSERT_EQ(round_tripped.number_of_variables(), 2);
    ASSERT_EQ(round_tripped.number_of_constraints(), 2);
    EXPECT_EQ(round_tripped.objective_direction, original.objective_direction);
    for (int v = 0; v < 2; ++v) {
        EXPECT_EQ(round_tripped.variables_lower_bounds[v], original.variables_lower_bounds[v]);
        EXPECT_EQ(round_tripped.variables_upper_bounds[v], original.variables_upper_bounds[v]);
        EXPECT_EQ(round_tripped.objective_coefficients[v], original.objective_coefficients[v]);
    }
    for (int c = 0; c < 2; ++c) {
        EXPECT_EQ(round_tripped.constraints_lower_bounds[c], original.constraints_lower_bounds[c]);
        EXPECT_EQ(round_tripped.constraints_upper_bounds[c], original.constraints_upper_bounds[c]);
    }
}

TEST(MathOptNl, WriteThenReadNonlinearAndQuadraticRoundTrip)
{
    // minimize exp(x0) + x1 * x2
    // s.t. c0: x0 + x1 + x2 <= 10
    // .nl reorders variables (nonlinear ones first), so compare
    // semantically via evaluate_objective/evaluate_constraint instead of
    // by index.
    MathOptModel original(3, 1, 3);
    original.objective_direction = ObjectiveDirection::Minimize;
    original.variables_lower_bounds = {0.0, 0.0, 0.0};
    original.variables_upper_bounds = {5.0, 5.0, 5.0};
    original.objective_quadratic_elements_variables_1 = {1};
    original.objective_quadratic_elements_variables_2 = {2};
    original.objective_quadratic_elements_coefficients = {1.0};
    original.objective_nonlinear_elements_operators = {'e', 'v'};
    original.objective_nonlinear_elements_values = {0.0, 0.0};
    original.objective_nonlinear_elements_variables = {-1, 0};
    original.objective_nonlinear_elements_number_of_children = {1, 0};
    original.constraints_starts = {0};
    original.elements_variables = {0, 1, 2};
    original.elements_coefficients = {1.0, 1.0, 1.0};
    original.constraints_lower_bounds = {-inf};
    original.constraints_upper_bounds = {10.0};

    MathOptModel round_tripped = read(write(original));

    ASSERT_EQ(round_tripped.number_of_variables(), 3);
    ASSERT_EQ(round_tripped.number_of_constraints(), 1);
    EXPECT_TRUE(round_tripped.check_solution({1.0, 2.0, 3.0}));

    std::vector<double> values = {0.5, 1.5, 2.5};
    EXPECT_NEAR(
            round_tripped.evaluate_objective(values),
            original.evaluate_objective(values),
            1e-9);
    EXPECT_NEAR(
            round_tripped.evaluate_constraint(values, 0),
            original.evaluate_constraint(values, 0),
            1e-9);
}

TEST(MathOptNl, WriteBlackBoxThrows)
{
    MathOptModel model(1, 0, 0);
    model.objective_direction = ObjectiveDirection::Minimize;
    model.objective_function = [](const std::vector<double>& x) -> BlackBoxFunctionOutput
    {
        BlackBoxFunctionOutput output;
        output.objective_value = x[0];
        return output;
    };

    std::ostringstream out;
    EXPECT_THROW(write_nl(model, out), std::invalid_argument);
}

TEST(MathOptNl, NanValueThrows)
{
    std::string header =
        "g0\n"
        "2 1 1 0 0\n"
        "0 0\n"
        "0 0\n"
        "0 0 0\n"
        "0 0\n"
        "0 0 0 0 0\n"
        "2 2\n"
        "0 0\n"
        "0 0 0 0 0\n";
    std::string body =
        "b\n"
        "0 0 2.5\n"
        "0 0 4.125\n"
        "r\n"
        "1 RHS\n"
        "C0\n"
        "nVALUE\n"
        "O0 0\n"
        "n0\n"
        "k1\n"
        "1\n"
        "J0 2\n"
        "0 0.5\n"
        "1 1.75\n"
        "G0 2\n"
        "0 1.5\n"
        "1 2.25\n";
    auto with = [&](const std::string& rhs, const std::string& value) {
        std::string text = body;
        text.replace(text.find("RHS"), 3, rhs);
        text.replace(text.find("VALUE"), 5, value);
        return header + text;
    };

    MathOptModel model = read(with("3.5", "0"));
    EXPECT_EQ(model.constraints_upper_bounds[0], 3.5);
    EXPECT_THROW(read(with("nan", "0")), std::invalid_argument);
    EXPECT_THROW(read(with("3.5abc", "0")), std::invalid_argument);
    EXPECT_THROW(read(with("3.5", "nan")), std::invalid_argument);
}

namespace
{

/** max 1.5 x + 2.25 y  s.t.  c: 0.5 x + 1.75 y <= 3.5, as an .nl file. */
std::string two_variable_nl()
{
    MathOptModel model;
    model.objective_direction = ObjectiveDirection::Maximize;
    model.variables_lower_bounds = {0.0, 0.0};
    model.variables_upper_bounds = {2.5, 4.125};
    model.variables_types = {VariableType::Continuous, VariableType::Continuous};
    model.objective_coefficients = {1.5, 2.25};
    model.constraints_lower_bounds = {-inf};
    model.constraints_upper_bounds = {3.5};
    model.constraints_starts = {0};
    model.elements_variables = {0, 1};
    model.elements_coefficients = {0.5, 1.75};
    return write(model);
}

}

TEST(MathOptNl, NamesFromColAndRowStreams)
{
    MathOptModel model = read(two_variable_nl());
    EXPECT_TRUE(model.variables_names.empty());

    // The .row file lists the constraints' names, then the objective's.
    std::istringstream col("x[1]\ny[1]  \r\n");
    std::istringstream row("capacity\nprofit\n");
    read_nl_names(model, &col, &row);
    EXPECT_EQ(model.variables_names, (std::vector<std::string>{"x[1]", "y[1]"}));
    EXPECT_EQ(model.constraints_names, (std::vector<std::string>{"capacity"}));
}

TEST(MathOptNl, MismatchedNameFilesAreIgnored)
{
    MathOptModel model = read(two_variable_nl());
    std::istringstream col("x\ny\nz\n");
    std::istringstream row("\n");
    read_nl_names(model, &col, &row);
    EXPECT_TRUE(model.variables_names.empty());
    EXPECT_TRUE(model.constraints_names.empty());

    std::istringstream col_only("x\ny\n");
    read_nl_names(model, &col_only, nullptr);
    EXPECT_EQ(model.variables_names, (std::vector<std::string>{"x", "y"}));
    EXPECT_TRUE(model.constraints_names.empty());
}

TEST(MathOptNl, NamesFilePaths)
{
    EXPECT_EQ(nl_col_path("dir/model.nl"), "dir/model.col");
    EXPECT_EQ(nl_row_path("dir/model.nl"), "dir/model.row");
    EXPECT_EQ(nl_col_path("model.NL"), "model.col");
    EXPECT_EQ(nl_col_path("model.nl.gz"), "model.col");
    EXPECT_EQ(nl_row_path("model.nl.bz2"), "model.row");
    EXPECT_EQ(nl_col_path("model"), "model.col");
}

TEST(MathOptNl, ReadFileWithNames)
{
    std::filesystem::path directory = std::filesystem::temp_directory_path() / "mathopt_nl_test_names";
    std::filesystem::create_directories(directory);
    std::string nl_path = (directory / "model.nl").string();
    { std::ofstream(nl_path, std::ios::binary) << two_variable_nl(); }
    std::filesystem::remove(directory / "model.col");
    std::filesystem::remove(directory / "model.row");

    EXPECT_TRUE(read_nl(nl_path).variables_names.empty());

    { std::ofstream(directory / "model.col") << "x\ny\n"; }
    { std::ofstream(directory / "model.row") << "c\nobj\n"; }
    MathOptModel model = read_nl(nl_path);
    EXPECT_EQ(model.variable_name(1), "y");
    EXPECT_EQ(model.constraint_name(0), "c");
    std::filesystem::remove_all(directory);
}

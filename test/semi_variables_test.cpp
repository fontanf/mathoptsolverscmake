#include "semi_test_models.hpp"

#include "mathoptsolverscmake/mathopt_json.hpp"
#include "mathoptsolverscmake/mathopt_lp.hpp"
#include "mathoptsolverscmake/mathopt_mps.hpp"
#include "mathoptsolverscmake/mathopt_nl.hpp"

#include <gtest/gtest.h>

#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace mathoptsolverscmake;
using namespace mathoptsolverscmaketest;

namespace
{

MathOptModel read_mps_string(const std::string& text)
{
    std::istringstream in(text);
    return read_mps(in);
}

MathOptModel read_lp_string(const std::string& text)
{
    std::istringstream in(text);
    return read_lp(in);
}

void expect_semi_model(const MathOptModel& model, VariableType x3_type)
{
    ASSERT_EQ(model.number_of_variables(), 4);
    EXPECT_EQ(model.variables_types[0], VariableType::Continuous);
    EXPECT_EQ(model.variables_types[2], x3_type);
    EXPECT_EQ(model.variables_types[3], VariableType::Continuous);
    EXPECT_EQ(model.variables_lower_bounds[2], 1.1);
    EXPECT_EQ(model.variables_upper_bounds[2], 10.0);
    EXPECT_TRUE(model.has_semi_variables());
    EXPECT_TRUE(model.is_milp());
}

template <typename Writer, typename Reader>
void expect_round_trip(const MathOptModel& original, Writer writer, Reader reader)
{
    std::ostringstream out;
    writer(original, out);
    std::istringstream in(out.str());
    MathOptModel model = reader(in);
    EXPECT_EQ(model.variables_types, original.variables_types) << out.str();
    EXPECT_EQ(model.variables_lower_bounds, original.variables_lower_bounds) << out.str();
    EXPECT_EQ(model.variables_upper_bounds, original.variables_upper_bounds) << out.str();
}

}

TEST(SemiVariables, TypeNames)
{
    for (VariableType type: {VariableType::SemiContinuous, VariableType::SemiInteger}) {
        std::ostringstream out;
        out << type;
        std::istringstream in(out.str());
        VariableType read_type = VariableType::Continuous;
        in >> read_type;
        EXPECT_EQ(read_type, type) << out.str();
    }
    EXPECT_TRUE(is_semi(VariableType::SemiContinuous));
    EXPECT_FALSE(is_integer(VariableType::SemiContinuous));
    EXPECT_TRUE(is_integer(VariableType::SemiInteger));
    EXPECT_FALSE(is_semi(VariableType::Integer));
}

TEST(SemiVariables, ReadMps)
{
    expect_semi_model(read_mps_string(semi_mps(false)), VariableType::SemiContinuous);
    expect_semi_model(read_mps_string(semi_mps(true)), VariableType::SemiInteger);
}

TEST(SemiVariables, ReadLp)
{
    expect_semi_model(read_lp_string(semi_lp(false)), VariableType::SemiContinuous);
    expect_semi_model(read_lp_string(semi_lp(true)), VariableType::SemiInteger);

    // "semi" / "semis" headers, and General before Semi-continuous.
    std::string text = semi_lp(false);
    text.replace(text.find("Semi-continuous"), 15, "General\nx3\nsemis");
    expect_semi_model(read_lp_string(text), VariableType::SemiInteger);
}

TEST(SemiVariables, RoundTrips)
{
    for (bool semi_integer: {false, true}) {
        MathOptModel original = read_mps_string(semi_mps(semi_integer));
        expect_round_trip(original,
                [](const MathOptModel& m, std::ostream& out) { write_mps(m, out); },
                [](std::istream& in) { return read_mps(in); });
        expect_round_trip(original,
                [](const MathOptModel& m, std::ostream& out) { write_lp(m, out); },
                [](std::istream& in) { return read_lp(in); });
        expect_round_trip(original,
                [](const MathOptModel& m, std::ostream& out) { write_json(m, out); },
                [](std::istream& in) { return read_json(in); });
    }

    // A semi-continuous variable with the default bounds [0, +inf), and a
    // semi-integer one with a negative lower bound.
    MathOptModel original = read_mps_string(semi_mps(false));
    original.variables_lower_bounds[2] = 0.0;
    original.variables_upper_bounds[2] = std::numeric_limits<double>::infinity();
    original.variables_types[3] = VariableType::SemiInteger;
    original.variables_lower_bounds[3] = -5.0;
    original.variables_upper_bounds[3] = 7.0;
    expect_round_trip(original,
            [](const MathOptModel& m, std::ostream& out) { write_mps(m, out); },
            [](std::istream& in) { return read_mps(in); });
    expect_round_trip(original,
            [](const MathOptModel& m, std::ostream& out) { write_lp(m, out); },
            [](std::istream& in) { return read_lp(in); });
}

TEST(SemiVariables, NlWriterThrows)
{
    MathOptModel model = read_mps_string(semi_mps(false));
    std::ostringstream out;
    EXPECT_THROW(write_nl(model, out), std::invalid_argument);
}

TEST(SemiVariables, CheckSolution)
{
    MathOptModel model = read_mps_string(semi_mps(true));
    model.feasibility_tolerance = 1e-9;
    model.integrality_tolerance = 1e-9;
    // x3: 0, or an integer in [1.1, 10].
    for (double x3: {0.0, 2.0, 10.0})
        EXPECT_TRUE(model.check_solution_variable({0, 0, x3, 0}, 2, 0)) << x3;
    for (double x3: {0.5, 1.0, 2.5, 11.0, -1.0})
        EXPECT_FALSE(model.check_solution_variable({0, 0, x3, 0}, 2, 0)) << x3;

    model.variables_types[2] = VariableType::SemiContinuous;
    for (double x3: {0.0, 1.1, 2.5, 10.0})
        EXPECT_TRUE(model.check_solution_variable({0, 0, x3, 0}, 2, 0)) << x3;
    for (double x3: {0.5, 1.0, 10.5})
        EXPECT_FALSE(model.check_solution_variable({0, 0, x3, 0}, 2, 0)) << x3;
}

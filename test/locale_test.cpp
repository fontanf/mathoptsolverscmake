#include "mathoptsolverscmake/mathopt_json.hpp"
#include "mathoptsolverscmake/mathopt_lp.hpp"
#include "mathoptsolverscmake/mathopt_mps.hpp"
#include "mathoptsolverscmake/mathopt_nl.hpp"

#include <gtest/gtest.h>

#include <clocale>
#include <cstring>
#include <functional>
#include <limits>
#include <sstream>
#include <string>

using namespace mathoptsolverscmake;

namespace
{

const double inf = std::numeric_limits<double>::infinity();

/** max 1.5 x + 2.25 y  s.t.  0.5 x + 1.75 y <= 3.5,  0 <= x <= 2.5,  0 <= y <= 4.125. */
MathOptModel fractional_model()
{
    MathOptModel model;
    model.objective_direction = ObjectiveDirection::Maximize;
    model.variables_names = {"x", "y"};
    model.variables_lower_bounds = {0.0, 0.0};
    model.variables_upper_bounds = {2.5, 4.125};
    model.variables_types = {VariableType::Continuous, VariableType::Continuous};
    model.objective_coefficients = {1.5, 2.25};
    model.constraints_names = {"c"};
    model.constraints_lower_bounds = {-inf};
    model.constraints_upper_bounds = {3.5};
    model.constraints_starts = {0};
    model.elements_variables = {0, 1};
    model.elements_coefficients = {0.5, 1.75};
    return model;
}

/** Sets the C locale to one with a decimal comma for the scope, if one is installed. */
class CommaLocale
{

public:

    CommaLocale()
    {
        const char* current = std::setlocale(LC_ALL, nullptr);
        saved_ = current? current: "C";
        for (const char* name: {"fr_FR.UTF-8", "fr_FR.utf8", "de_DE.UTF-8", "de_DE.utf8", "French_France.1252"}) {
            if (std::setlocale(LC_ALL, name) != nullptr
                    && std::strcmp(std::localeconv()->decimal_point, ",") == 0) {
                active_ = true;
                return;
            }
        }
        std::setlocale(LC_ALL, saved_.c_str());
    }

    ~CommaLocale() { std::setlocale(LC_ALL, saved_.c_str()); }

    bool active() const { return active_; }

private:

    std::string saved_;
    bool active_ = false;

};

void check_read_back(
        const std::function<void(const MathOptModel&, std::ostream&)>& writer,
        const std::function<MathOptModel(std::istream&)>& reader)
{
    MathOptModel original = fractional_model();
    std::ostringstream out;
    writer(original, out);

    CommaLocale locale;
    if (!locale.active())
        GTEST_SKIP() << "no locale with a decimal comma installed";
    std::istringstream in(out.str());
    MathOptModel model = reader(in);

    EXPECT_EQ(model.objective_coefficients, original.objective_coefficients);
    EXPECT_EQ(model.variables_upper_bounds, original.variables_upper_bounds);
    EXPECT_EQ(model.constraints_upper_bounds, original.constraints_upper_bounds);
    EXPECT_EQ(model.elements_coefficients, original.elements_coefficients);
}

}

TEST(Locale, MpsReaderIgnoresLocale)
{
    check_read_back(
            [](const MathOptModel& model, std::ostream& out) { write_mps(model, out); },
            [](std::istream& in) { return read_mps(in); });
}

TEST(Locale, LpReaderIgnoresLocale)
{
    check_read_back(
            [](const MathOptModel& model, std::ostream& out) { write_lp(model, out); },
            [](std::istream& in) { return read_lp(in); });
}

TEST(Locale, JsonReaderIgnoresLocale)
{
    check_read_back(
            [](const MathOptModel& model, std::ostream& out) { write_json(model, out); },
            [](std::istream& in) { return read_json(in); });
}

TEST(Locale, NlReaderIgnoresLocale)
{
    check_read_back(
            [](const MathOptModel& model, std::ostream& out) { write_nl(model, out); },
            [](std::istream& in) { return read_nl(in); });
}

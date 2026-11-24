/**
 * Python bindings for the solver-independent part of MathOptSolversCMake:
 * MathOptModel and the MPS / LP / JSON / .nl readers and writers.
 *
 * MathOptModel's arrays are exposed as read-only NumPy views onto the C++
 * vectors (no copy); nanobind's reference_internal policy (the default for
 * properties) makes each view keep its model alive. Names are exposed as
 * lists of str, and expression-tree operators as bytes. Readers take the
 * file contents as bytes (see read() in __init__.py, which also handles
 * .gz/.bz2 decompression) and release the GIL while parsing.
 */
#include "mathoptsolverscmake/mathopt.hpp"
#include "mathoptsolverscmake/mathopt_mps.hpp"
#include "mathoptsolverscmake/mathopt_lp.hpp"
#include "mathoptsolverscmake/mathopt_json.hpp"
#include "mathoptsolverscmake/mathopt_nl.hpp"
#include "mathoptsolverscmake/dual.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <sstream>

namespace nb = nanobind;
using namespace nb::literals;
using namespace mathoptsolverscmake;

namespace
{

template <typename T>
using ReadOnlyArray = nb::ndarray<nb::numpy, const T, nb::ndim<1>>;

/** A read-only NumPy view onto one of the model's vectors. */
template <typename T, std::vector<T> MathOptModel::*Member>
ReadOnlyArray<T> view_member(const MathOptModel& model)
{
    const std::vector<T>& vector = model.*Member;
    return ReadOnlyArray<T>(vector.data(), {vector.size()});
}

/** Variable types as a NumPy array of their integer values (see VariableType). */
nb::ndarray<nb::numpy, int8_t, nb::ndim<1>> variables_types_array(const MathOptModel& model)
{
    size_t n = model.variables_types.size();
    int8_t* data = new int8_t[n];
    for (size_t i = 0; i < n; ++i)
        data[i] = static_cast<int8_t>(model.variables_types[i]);
    nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<int8_t*>(p); });
    return nb::ndarray<nb::numpy, int8_t, nb::ndim<1>>(data, {n}, owner);
}

/** Constraint classes as a NumPy array of their integer values (see ConstraintClass). */
nb::ndarray<nb::numpy, int8_t, nb::ndim<1>> constraints_classes_array(const MathOptModel& model)
{
    size_t n = model.number_of_constraints();
    int8_t* data = new int8_t[n];
    for (size_t i = 0; i < n; ++i)
        data[i] = static_cast<int8_t>(model.constraint_class(i));
    nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<int8_t*>(p); });
    return nb::ndarray<nb::numpy, int8_t, nb::ndim<1>>(data, {n}, owner);
}

nb::bytes to_bytes(const std::vector<char>& operators)
{
    return nb::bytes(operators.data(), operators.size());
}

/** Parse `data` with `reader`, without holding the GIL. */
MathOptModel read_from_bytes(
        MathOptModel (*reader)(std::istream&),
        nb::bytes data)
{
    std::string contents(data.c_str(), data.size());
    nb::gil_scoped_release release;
    std::istringstream in(contents, std::ios::in | std::ios::binary);
    return reader(in);
}

std::string write_to_string(
        void (*writer)(const MathOptModel&, std::ostream&),
        const MathOptModel& model)
{
    std::ostringstream out;
    writer(model, out);
    return out.str();
}

}

NB_MODULE(_core, m)
{
    m.doc() = "Python bindings for MathOptSolversCMake's model and file formats.";

    nb::enum_<ObjectiveDirection>(m, "ObjectiveDirection", nb::is_arithmetic())
        .value("Minimize", ObjectiveDirection::Minimize)
        .value("Maximize", ObjectiveDirection::Maximize);

    nb::enum_<VariableType>(m, "VariableType", nb::is_arithmetic())
        .value("Continuous", VariableType::Continuous)
        .value("Binary", VariableType::Binary)
        .value("Integer", VariableType::Integer);

    nb::enum_<ConstraintSense>(m, "ConstraintSense", nb::is_arithmetic())
        .value("LessThanOrEqualTo", ConstraintSense::LessThanOrEqualTo)
        .value("GreaterThanOrEqualTo", ConstraintSense::GreaterThanOrEqualTo)
        .value("Equality", ConstraintSense::Equality)
        .value("Range", ConstraintSense::Range)
        .value("Free", ConstraintSense::Free);

    nb::enum_<ConstraintClass>(m, "ConstraintClass", nb::is_arithmetic())
        .value("Empty", ConstraintClass::Empty)
        .value("Free", ConstraintClass::Free)
        .value("Singleton", ConstraintClass::Singleton)
        .value("Aggregation", ConstraintClass::Aggregation)
        .value("Precedence", ConstraintClass::Precedence)
        .value("VariableBound", ConstraintClass::VariableBound)
        .value("SetPartitioning", ConstraintClass::SetPartitioning)
        .value("SetPacking", ConstraintClass::SetPacking)
        .value("SetCovering", ConstraintClass::SetCovering)
        .value("Cardinality", ConstraintClass::Cardinality)
        .value("InvariantKnapsack", ConstraintClass::InvariantKnapsack)
        .value("EquationKnapsack", ConstraintClass::EquationKnapsack)
        .value("Binpacking", ConstraintClass::Binpacking)
        .value("Knapsack", ConstraintClass::Knapsack)
        .value("IntegerKnapsack", ConstraintClass::IntegerKnapsack)
        .value("MixedBinary", ConstraintClass::MixedBinary)
        .value("GeneralLinear", ConstraintClass::GeneralLinear)
        .value("GeneralQuadratic", ConstraintClass::GeneralQuadratic)
        .value("GeneralBlackBox", ConstraintClass::GeneralBlackBox)
        .value("GeneralNonlinear", ConstraintClass::GeneralNonlinear);

    nb::class_<MathOptModel>(m, "Model",
            "An optimization model. Its arrays are read-only views; see mathopt.hpp "
            "for the meaning of each field.")
        .def_prop_ro("objective_direction",
                [](const MathOptModel& model) { return model.objective_direction; })
        .def_prop_ro("number_of_variables",
                [](const MathOptModel& model) { return model.number_of_variables(); })
        .def_prop_ro("number_of_constraints",
                [](const MathOptModel& model) { return model.number_of_constraints(); })
        .def_prop_ro("number_of_elements",
                [](const MathOptModel& model) { return model.number_of_elements(); })

        // Variables.
        .def_prop_ro("variables_lower_bounds",
                view_member<double, &MathOptModel::variables_lower_bounds>)
        .def_prop_ro("variables_upper_bounds",
                view_member<double, &MathOptModel::variables_upper_bounds>)
        .def_prop_ro("variables_types", variables_types_array, nb::rv_policy::move,
                "Variable types as an int8 array of VariableType values.")
        .def_prop_ro("variables_names",
                [](const MathOptModel& model) { return model.variables_names; },
                "Variable names (empty if the file format doesn't store them).")
        .def_prop_ro("variables_initial_values",
                view_member<double, &MathOptModel::variables_initial_values>)
        .def("variable_name", &MathOptModel::variable_name, "variable_id"_a)

        // Constraints.
        .def_prop_ro("constraints_lower_bounds",
                view_member<double, &MathOptModel::constraints_lower_bounds>)
        .def_prop_ro("constraints_upper_bounds",
                view_member<double, &MathOptModel::constraints_upper_bounds>)
        .def_prop_ro("constraints_names",
                [](const MathOptModel& model) { return model.constraints_names; },
                "Constraint names (empty if the file format doesn't store them).")
        .def("constraint_name", &MathOptModel::constraint_name, "constraint_id"_a)
        .def("constraint_sense", &MathOptModel::constraint_sense, "constraint_id"_a)
        .def("constraint_class", &MathOptModel::constraint_class, "constraint_id"_a)
        .def_prop_ro("constraints_classes", constraints_classes_array, nb::rv_policy::move,
                "Constraint classes as an int8 array of ConstraintClass values.")
        .def("constraint_end", &MathOptModel::constraint_end, "constraint_id"_a)
        .def("quadratic_constraint_end", &MathOptModel::quadratic_constraint_end, "constraint_id"_a)
        .def("nonlinear_constraint_end", &MathOptModel::nonlinear_constraint_end, "constraint_id"_a)
        .def("format_constraint",
                [](const MathOptModel& model, int constraint_id) {
                    std::ostringstream os;
                    model.format_constraint(os, constraint_id);
                    return os.str();
                },
                "constraint_id"_a)

        // Linear structures.
        .def_prop_ro("objective_coefficients",
                view_member<double, &MathOptModel::objective_coefficients>)
        .def_prop_ro("constraints_starts",
                view_member<int, &MathOptModel::constraints_starts>)
        .def_prop_ro("elements_variables",
                view_member<int, &MathOptModel::elements_variables>)
        .def_prop_ro("elements_coefficients",
                view_member<double, &MathOptModel::elements_coefficients>)

        // Quadratic structures.
        .def_prop_ro("objective_quadratic_elements_variables_1",
                view_member<int, &MathOptModel::objective_quadratic_elements_variables_1>)
        .def_prop_ro("objective_quadratic_elements_variables_2",
                view_member<int, &MathOptModel::objective_quadratic_elements_variables_2>)
        .def_prop_ro("objective_quadratic_elements_coefficients",
                view_member<double, &MathOptModel::objective_quadratic_elements_coefficients>)
        .def_prop_ro("quadratic_elements_constraints_starts",
                view_member<int, &MathOptModel::quadratic_elements_constraints_starts>)
        .def_prop_ro("quadratic_elements_variables_1",
                view_member<int, &MathOptModel::quadratic_elements_variables_1>)
        .def_prop_ro("quadratic_elements_variables_2",
                view_member<int, &MathOptModel::quadratic_elements_variables_2>)
        .def_prop_ro("quadratic_elements_coefficients",
                view_member<double, &MathOptModel::quadratic_elements_coefficients>)

        // Nonlinear structures.
        .def_prop_ro("objective_nonlinear_elements_operators",
                [](const MathOptModel& model) { return to_bytes(model.objective_nonlinear_elements_operators); },
                "Operator char codes of the objective's expression tree (see mathopt.hpp), as bytes.")
        .def_prop_ro("objective_nonlinear_elements_values",
                view_member<double, &MathOptModel::objective_nonlinear_elements_values>)
        .def_prop_ro("objective_nonlinear_elements_variables",
                view_member<int, &MathOptModel::objective_nonlinear_elements_variables>)
        .def_prop_ro("objective_nonlinear_elements_number_of_children",
                view_member<int, &MathOptModel::objective_nonlinear_elements_number_of_children>)
        .def_prop_ro("nonlinear_elements_constraints_starts",
                view_member<int, &MathOptModel::nonlinear_elements_constraints_starts>)
        .def_prop_ro("nonlinear_elements_operators",
                [](const MathOptModel& model) { return to_bytes(model.nonlinear_elements_operators); },
                "Operator char codes of the constraints' expression trees (see mathopt.hpp), as bytes.")
        .def_prop_ro("nonlinear_elements_values",
                view_member<double, &MathOptModel::nonlinear_elements_values>)
        .def_prop_ro("nonlinear_elements_variables",
                view_member<int, &MathOptModel::nonlinear_elements_variables>)
        .def_prop_ro("nonlinear_elements_number_of_children",
                view_member<int, &MathOptModel::nonlinear_elements_number_of_children>)

        // Defined variables (common subexpressions).
        .def_prop_ro("number_of_defined_variables",
                [](const MathOptModel& model) { return model.number_of_defined_variables(); })
        .def("defined_variable_end", &MathOptModel::defined_variable_end, "defined_variable_id"_a)
        .def_prop_ro("defined_variables_starts",
                view_member<int, &MathOptModel::defined_variables_starts>)
        .def_prop_ro("defined_variables_elements_operators",
                [](const MathOptModel& model) { return to_bytes(model.defined_variables_elements_operators); },
                "Operator char codes of the defined variables' expression trees (see mathopt.hpp), as bytes.")
        .def_prop_ro("defined_variables_elements_values",
                view_member<double, &MathOptModel::defined_variables_elements_values>)
        .def_prop_ro("defined_variables_elements_variables",
                view_member<int, &MathOptModel::defined_variables_elements_variables>)
        .def_prop_ro("defined_variables_elements_number_of_children",
                view_member<int, &MathOptModel::defined_variables_elements_number_of_children>)

        // Model kind.
        .def("has_non_continuous_variables", &MathOptModel::has_non_continuous_variables)
        .def("has_defined_variables", &MathOptModel::has_defined_variables)
        .def("has_quadratic", &MathOptModel::has_quadratic)
        .def("has_nonlinear", &MathOptModel::has_nonlinear)
        .def("is_lp", &MathOptModel::is_lp)
        .def("is_milp", &MathOptModel::is_milp)

        // Evaluations and checks.
        .def("evaluate_objective", &MathOptModel::evaluate_objective, "solution"_a)
        .def("evaluate_constraint", &MathOptModel::evaluate_constraint, "solution"_a, "constraint_id"_a)
        .def("check", &MathOptModel::check, "verbosity_level"_a = 0,
                "Check that the model's arrays are consistent.")
        .def("check_solution", &MathOptModel::check_solution, "solution"_a, "verbosity_level"_a = 0)
        .def("format",
                [](const MathOptModel& model, int verbosity_level) {
                    std::ostringstream os;
                    model.format(os, verbosity_level);
                    return os.str();
                },
                "verbosity_level"_a = 1)
        .def("__repr__",
                [](const MathOptModel& model) {
                    return "<mathoptmodel.Model: "
                        + std::to_string(model.number_of_variables()) + " variables, "
                        + std::to_string(model.number_of_constraints()) + " constraints>";
                });

    // Readers, from the raw file contents.
    m.def("read_mps_bytes",
            [](nb::bytes data) { return read_from_bytes(read_mps, data); },
            "data"_a, "Read a model from the contents of a free-format MPS file.");
    m.def("read_lp_bytes",
            [](nb::bytes data) { return read_from_bytes(read_lp, data); },
            "data"_a, "Read a model from the contents of a CPLEX-style .lp file.");
    m.def("read_json_bytes",
            [](nb::bytes data) { return read_from_bytes(read_json, data); },
            "data"_a, "Read a model from the contents of a MathOptSolversCMake JSON file.");
    m.def("read_nl_bytes",
            [](nb::bytes data) { return read_from_bytes(read_nl, data); },
            "data"_a, "Read a model from the contents of an AMPL .nl file.");

    // Writers, to a str.
    m.def("write_mps_string",
            [](const MathOptModel& model) { return write_to_string(write_mps, model); },
            "model"_a, "Write a model as a free-format MPS file, returned as a str.");
    m.def("write_lp_string",
            [](const MathOptModel& model) { return write_to_string(write_lp, model); },
            "model"_a, "Write a model as a CPLEX-style .lp file, returned as a str.");
    m.def("write_json_string",
            [](const MathOptModel& model) { return write_to_string(write_json, model); },
            "model"_a, "Write a model as a MathOptSolversCMake JSON file, returned as a str.");
    m.def("write_nl_string",
            [](const MathOptModel& model) { return write_to_string(write_nl, model); },
            "model"_a, "Write a model as an AMPL .nl file (text format), returned as a str.");

    m.def("dual", &dual, "model"_a, "Build the dual of a linear program.");
    m.def("expand_defined_variables", &expand_defined_variables, "model"_a,
            "Return an equivalent model whose expression trees reference no defined variable.");
}

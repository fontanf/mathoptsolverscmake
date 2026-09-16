#pragma once

#include "mathoptsolverscmake/mathopt.hpp"

#include <istream>
#include <ostream>
#include <string>

namespace mathoptsolverscmake
{

/**
 * Build a MathOptModel from JSON that directly mirrors its fields: one
 * top-level object whose keys are (a subset of) MathOptModel's member
 * names, each holding that member's value straight (a double/int array
 * stays an array, ObjectiveDirection/VariableType are written as the
 * strings produced by their operator<<, e.g. "Minimize"/"Continuous", and
 * a std::vector<char> of nonlinear-expression operators is an array of
 * single-character strings). A field absent from the JSON is left at its
 * default (an empty vector, or 0.0 for the tolerances); "objective_direction"
 * is the only required field. Black-box objective/constraint functions
 * (std::function) have no JSON representation and are never read.
 *
 * Infinite bounds are written/read as the bare (non-standard, but widely
 * accepted) JSON tokens Infinity / -Infinity / NaN rather than as quoted
 * strings.
 *
 * This is a structural dump, not a solver interchange format like MPS:
 * fields are copied as given with no cross-field consistency checks, so
 * malformed combinations (e.g. mismatched array sizes) simply produce a
 * MathOptModel that fails MathOptModel::check().
 */
MathOptModel read_json(const std::string& file_path);

/** Same as read_json(file_path), reading from an already-open stream. */
MathOptModel read_json(std::istream& in);

/**
 * Write a MathOptModel as JSON, field by field (see read_json()'s doc
 * comment for the mapping). A field that is empty/at its default (name
 * arrays, quadratic/nonlinear structures, ...) is omitted rather than
 * written as an empty array. Throws std::invalid_argument if the model has
 * a black-box objective or constraint function, since those can't be
 * represented.
 *
 * write_json() followed by read_json() round-trips exactly.
 */
void write_json(
        const MathOptModel& model,
        const std::string& file_path);

/** Same as write_json(model, file_path), writing to an already-open stream. */
void write_json(
        const MathOptModel& model,
        std::ostream& out);

}

#pragma once

#include "mathoptsolverscmake/mathopt.hpp"

#include <istream>
#include <ostream>
#include <string>

namespace mathoptsolverscmake
{

/**
 * Build a MathOptModel from an AMPL .nl file (ASCII "g" format only; the
 * "b"/"h" binary encodings are not supported and cause an
 * std::invalid_argument to be thrown).
 *
 * Unlike MPS/LP, .nl is expression-tree-native, so linear, quadratic and
 * general nonlinear objectives/constraints are all supported: its
 * "o<opcode>"/"v<index>"/"n<value>" prefix-notation opcodes are translated
 * to/from MathOptModel's own char-opcode tree (see mathopt.hpp), for
 * exactly the operators MathOptModel itself supports (+, -, *, /, unary
 * negate, exp, log, sqrt, sin, cos, tan, pow -- both o0/o54 forms of +).
 * Any other opcode (min/max, floor/ceil, logical/comparison operators,
 * piecewise-linear terms, ...) causes an std::invalid_argument to be
 * thrown, as does a nonzero header count for anything else unsupported:
 * user-defined functions, common subexpressions, logical or
 * complementarity constraints, and network ("arc") variables. Suffix
 * ('S') segments are skipped, not stored, and so are initial dual values
 * ('d'); initial primal values ('x') populate variables_initial_values
 * (a variable not covered by a sparse 'x' segment defaults to 0).
 *
 * Variable/constraint names are not part of .nl's core structure (AMPL
 * writes those to separate .col/.row files) and are never read; the
 * model's variables_names/constraints_names are left empty. Black-box
 * functions have no .nl representation and are of course never read
 * either.
 */
MathOptModel read_nl(const std::string& file_path);

/** Same as read_nl(file_path), reading from an already-open stream. */
MathOptModel read_nl(std::istream& in);

/**
 * Write a MathOptModel as an AMPL .nl file (ASCII "g" format).
 *
 * Unlike MPS/LP, quadratic and nonlinear terms are supported (they are
 * .nl's native structure); only a black-box objective/constraint function
 * is unsupported (throws std::invalid_argument, i.e. unless
 * !model.has_black_box()). A quadratic term is folded into the written
 * nonlinear expression as a product (.nl has no separate quadratic
 * segment), so it comes back from read_nl() as a nonlinear term rather
 * than a quadratic one -- semantically identical, structurally different.
 * Variable/constraint names are not written (see read_nl()'s doc
 * comment); a variables_initial_values is written as a dense 'x' segment
 * when present.
 *
 * .nl requires variables partitioned into a specific order (nonlinear in
 * both the objective and constraints, then constraints-only, then
 * objective-only, then linear -- continuous before integer/binary within
 * each group), so write_nl() reindexes variables (and puts nonlinear
 * constraints first) to satisfy this. The model that comes back from
 * read_nl() is therefore equivalent -- same variables, bounds, types and
 * coefficients -- but not necessarily in the same variable/constraint
 * order as the original.
 */
void write_nl(
        const MathOptModel& model,
        const std::string& file_path);

/** Same as write_nl(model, file_path), writing to an already-open stream. */
void write_nl(
        const MathOptModel& model,
        std::ostream& out);

}

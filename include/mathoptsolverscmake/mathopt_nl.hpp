#pragma once

#include "mathoptsolverscmake/mathopt.hpp"

#include <istream>
#include <ostream>
#include <string>

namespace mathoptsolverscmake
{

/**
 * Build a MathOptModel from an AMPL .nl file, in either the text ("g") or
 * the binary ("b") encoding (both little- and big-endian, per the header's
 * "arith" field). The obsolete "h" encoding is not supported and causes an
 * std::invalid_argument to be thrown.
 *
 * Unlike MPS/LP, .nl is expression-tree-native, so linear, quadratic and
 * general nonlinear objectives/constraints are all supported: its
 * "o<opcode>"/"v<index>"/"n<value>" prefix-notation opcodes are translated
 * to/from MathOptModel's own char-opcode tree (see mathopt.hpp), for
 * exactly the operators MathOptModel itself supports (+, -, *, /, unary
 * negate, exp, log, log10, sqrt, sin, cos, tan, tanh, abs, pow, n-ary min
 * and max -- both o0/o54 forms of +). Any other opcode (floor/ceil,
 * logical/comparison operators, piecewise-linear terms, ...) causes an
 * std::invalid_argument to be thrown, as does a nonzero header count for
 * anything else unsupported: user-defined functions, logical or
 * complementarity constraints, and network ("arc") variables. Common
 * subexpressions ('V' segments) are read as MathOptModel defined variables,
 * numbered in the order they are defined in the file. Suffix
 * ('S') segments are skipped, not stored, and so are initial dual values
 * ('d'); initial primal values ('x') populate variables_initial_values
 * (a variable not covered by a sparse 'x' segment defaults to 0).
 *
 * Variable/constraint names are not part of .nl's core structure: AMPL
 * writes them to separate .col/.row files, next to the .nl file. read_nl()
 * reads them if they exist (see read_nl_names()); read_nl(std::istream&)
 * doesn't, and leaves the model's variables_names/constraints_names empty.
 * Black-box functions have no .nl representation and are of course never
 * read either.
 */
MathOptModel read_nl(const std::string& file_path);

/** Same as read_nl(file_path), reading from an already-open stream (without names). */
MathOptModel read_nl(std::istream& in);

/**
 * Set the names of a model read from an .nl file from the .col / .row
 * files AMPL writes next to it ("model.nl", "model.col", "model.row"):
 * .col has one variable name per line, .row one constraint name per line,
 * followed by the objectives' names. Either stream may be null.
 *
 * A file whose number of names doesn't match the model (fewer constraint
 * names than constraints, or a number of variable names other than the
 * number of variables) is ignored, as names don't change the model: it may
 * be left over from another version of the model.
 */
void read_nl_names(
        MathOptModel& model,
        std::istream* col,
        std::istream* row);

/**
 * The paths of the .col and .row files of an .nl file: its path with its
 * ".nl" extension (in any case; ignoring a trailing ".gz" or ".bz2")
 * replaced by ".col" and ".row".
 */
std::string nl_col_path(const std::string& nl_path);
std::string nl_row_path(const std::string& nl_path);

/**
 * Write a MathOptModel as an AMPL .nl file (ASCII "g" format).
 *
 * Unlike MPS/LP, quadratic and nonlinear terms are supported (they are
 * .nl's native structure); defined variables are expanded (see
 * expand_defined_variables()) rather than written as 'V' segments; only a black-box objective/constraint function
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

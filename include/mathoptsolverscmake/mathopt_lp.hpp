#pragma once

#include "mathoptsolverscmake/mathopt.hpp"

#include <istream>
#include <ostream>
#include <string>

namespace mathoptsolverscmake
{

/**
 * Build a MathOptModel from a CPLEX-style .lp file.
 *
 * Supports the objective section (Minimize/Maximize, with the usual
 * Min/Max/Minimum/Maximum spelling variants), Subject To / Such That / st /
 * s.t., Bounds (including "x free", "x = v", and both the "v <= x <= v" and
 * "v >= x >= v" range forms), General/Generals/Gen/Integer/Integers and
 * Binary/Binaries/Bin sections, terminated by an optional End. '\' starts a
 * line comment, and a C-style slash-star block comment is also supported.
 * Unlike MPS's
 * INTORG/INTEND marker, a variable's section placement directly gives its
 * VariableType (no bounds-based Binary/Integer inference): a variable
 * listed under Binary is Binary, one listed under General/Integer is
 * Integer, regardless of its bounds. Quadratic terms ('[' ... ']') and SOS
 * / semi-continuous / semi-integer sections are not supported and cause an
 * std::invalid_argument to be thrown.
 *
 * A variable/constraint name may contain letters, digits, and any of
 * "_.()#&%!?@{}$,;'" (a stricter, unambiguous subset of what real .lp
 * files may contain: e.g. names with an embedded '-' aren't supported,
 * since '-' is also the subtraction operator); such a name causes
 * read_lp() to misinterpret it, since the reader has no way to tell a
 * name apart from an expression. A constraint name may additionally be a
 * plain number (e.g. row "2"), matching HiGHS's own .lp reader; a
 * variable name may not, since a bare number in an expression is always
 * read as a coefficient.
 *
 * Infinite bounds are read from (and, in write_lp(), written as) the
 * bare, case-insensitive keyword "infinity" (or "inf"), following the
 * same convention as e.g. HiGHS's own .lp reader.
 */
MathOptModel read_lp(const std::string& file_path);

/** Same as read_lp(file_path), reading from an already-open stream. */
MathOptModel read_lp(std::istream& in);

/**
 * Write a MathOptModel as a .lp file (see read_lp()'s doc comment for the
 * supported subset). Only linear / mixed-integer-linear models are
 * supported (throws std::invalid_argument otherwise, i.e. unless
 * model.is_milp()), and variable/constraint names are restricted as
 * described above (throws otherwise).
 *
 * A Range constraint is written using the two-sided "lower <= expr <=
 * upper" form, and a Free constraint (both bounds infinite, which .lp has
 * no direct syntax for) as "expr >= -infinity". Both round-trip exactly
 * through read_lp(), though the Free-constraint encoding in particular is
 * specific to this reader/writer pair and not a general .lp convention.
 *
 * write_lp() followed by read_lp() round-trips exactly.
 */
void write_lp(
        const MathOptModel& model,
        const std::string& file_path);

/** Same as write_lp(model, file_path), writing to an already-open stream. */
void write_lp(
        const MathOptModel& model,
        std::ostream& out);

}

#pragma once

#include "mathoptsolverscmake/mathopt.hpp"

#include <istream>
#include <ostream>
#include <string>

namespace mathoptsolverscmake
{

/**
 * Build a MathOptModel from a free-format MPS file.
 *
 * Supports the NAME, OBJSENSE, ROWS, COLUMNS (including 'MARKER'
 * INTORG/INTEND integer sections), RHS, RANGES and BOUNDS (UP, LO, FX, MI,
 * PL, BV, FR, LI, UI) sections of a linear / mixed-integer-linear model.
 * Parsing stops at ENDATA. Rows other than the first 'N' row (the
 * objective) are free rows and are dropped, matching how solvers such as
 * HiGHS treat them. Quadratic and other extension sections (QUADOBJ,
 * QCMATRIX, SOS, SC, SI, ...) are not supported and cause an
 * std::invalid_argument to be thrown.
 *
 * Row/bound/marker conventions (default variable bounds, integer columns
 * defaulting to [0, 1] unless a BOUNDS entry overrides them, RANGES
 * translation, ...) follow HiGHS's free-format MPS reader
 * (io/HMpsFF.cpp), so a model built with read_mps() from a given file is
 * equivalent to the one obtained via Highs::readModel() on the same file
 * followed by mathoptsolverscmake::to_mathopt().
 */
MathOptModel read_mps(const std::string& file_path);

/** Same as read_mps(file_path), reading from an already-open stream. */
MathOptModel read_mps(std::istream& in);

/**
 * Write a MathOptModel to a free-format MPS file.
 *
 * Only linear / mixed-integer-linear models are supported (throws
 * std::invalid_argument otherwise, i.e. unless model.is_milp()). Writes
 * OBJSENSE only when maximizing, ROWS/COLUMNS/RHS/BOUNDS and, for range
 * constraints, RANGES. Every variable is guaranteed to appear in COLUMNS
 * (even with all-zero coefficients) so that read_mps() recovers exactly
 * model.number_of_variables() columns. Variable and constraint names
 * containing whitespace can't be represented in free-format MPS and cause
 * an std::invalid_argument to be thrown.
 *
 * write_mps() followed by read_mps() round-trips exactly, except: a Free
 * constraint (both bounds infinite) is written as a second 'N' row, which
 * read_mps() (like HiGHS) drops as a free row; and an Integer variable
 * whose bounds happen to be exactly [0, 1] is indistinguishable in MPS
 * from a declared Binary one and comes back as Binary (see read_mps()'s
 * doc comment).
 */
void write_mps(
        const MathOptModel& model,
        const std::string& file_path);

/** Same as write_mps(model, file_path), writing to an already-open stream. */
void write_mps(
        const MathOptModel& model,
        std::ostream& out);

}

#include "mathoptsolverscmake/mathopt_mps.hpp"

#include <cmath>
#include <cstdlib>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

using namespace mathoptsolverscmake;

namespace
{

enum class Section
{
    None,
    Name,
    Objsense,
    Rows,
    Columns,
    Rhs,
    Ranges,
    Bounds,
};

std::string to_upper(const std::string& s)
{
    std::string res = s;
    for (char& c: res)
        c = (char)std::toupper((unsigned char)c);
    return res;
}

std::vector<std::string> split_line(const std::string& line)
{
    std::vector<std::string> tokens;
    std::istringstream iss(line);
    std::string token;
    while (iss >> token)
        tokens.push_back(token);
    return tokens;
}

/** Parse an MPS number, translating the Fortran 'D' exponent to 'E'. */
double parse_number(const std::string& word)
{
    std::string w = word;
    for (char& c: w) {
        if (c == 'D' || c == 'd')
            c = 'E';
    }
    return std::strtod(w.c_str(), nullptr);
}

/**
 * Returns Section::None if 'word_upper' is not a recognized section
 * keyword, throws for a recognized-but-unsupported one.
 */
Section section_from_keyword(const std::string& word_upper)
{
    if (word_upper == "NAME") return Section::Name;
    if (word_upper == "OBJSENSE") return Section::Objsense;
    if (word_upper == "ROWS") return Section::Rows;
    if (word_upper == "COLUMNS") return Section::Columns;
    if (word_upper == "RHS") return Section::Rhs;
    if (word_upper == "RANGES") return Section::Ranges;
    if (word_upper == "BOUNDS") return Section::Bounds;
    return Section::None;
}

bool is_unsupported_section_keyword(const std::string& word_upper)
{
    return word_upper == "QSECTION"
        || word_upper == "QMATRIX"
        || word_upper == "QUADOBJ"
        || word_upper == "QCMATRIX"
        || word_upper == "CSECTION"
        || word_upper == "SETS"
        || word_upper == "SOS"
        || word_upper == "INDICATORS"
        || word_upper == "GENCONS"
        || word_upper == "PWLOBJ"
        || word_upper == "PWLNAM"
        || word_upper == "PWLCON"
        || word_upper == "DELAYEDROWS"
        || word_upper == "MODELCUTS"
        || word_upper == "USERCUTS";
}

bool has_whitespace(const std::string& s)
{
    for (char c: s)
        if (std::isspace((unsigned char)c))
            return true;
    return false;
}

void check_mps_name(const std::string& name)
{
    if (has_whitespace(name)) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "name \"" + name + "\" contains whitespace and cannot be "
                "represented in free-format MPS.");
    }
}

/**
 * All the bookkeeping needed while a MathOptModel is being built up
 * incrementally from an MPS file (row/column name lookups, dedup flags,
 * the per-row nonzero buffers COLUMNS fills in column order but that must
 * end up grouped by row, ...).
 */
struct MpsParseState
{
    const double inf = std::numeric_limits<double>::infinity();

    MathOptModel model;

    // Rows: -1 is the objective row, -2 is any other (free) 'N' row.
    std::unordered_map<std::string, int> row_ids;
    bool has_objective_row = false;
    std::vector<char> row_types; // 'L', 'G' or 'E', parallel to model constraints.
    std::vector<bool> row_has_rhs;
    std::vector<bool> row_has_range;
    std::vector<std::vector<std::pair<int, double>>> row_entries; // (variable_id, coefficient), parallel to model constraints.

    // Columns.
    std::unordered_map<std::string, int> column_ids;
    std::vector<bool> column_has_cost;
    std::vector<bool> column_is_integer;
    // True for an integer column until it receives any BOUNDS entry; such
    // columns default to [0, 1] (classic MPS convention, also HiGHS's).
    std::vector<bool> column_still_default_binary;
    std::vector<bool> column_has_lower;
    std::vector<bool> column_has_upper;

    bool in_integer_marker_section = false;

    int get_or_create_column(const std::string& name)
    {
        auto it = column_ids.find(name);
        if (it != column_ids.end())
            return it->second;

        int column_id = model.number_of_variables();
        column_ids[name] = column_id;
        model.variables_names.push_back(name);
        model.variables_lower_bounds.push_back(0.0);
        model.variables_upper_bounds.push_back(inf);
        model.objective_coefficients.push_back(0.0);
        column_has_cost.push_back(false);
        column_is_integer.push_back(in_integer_marker_section);
        column_still_default_binary.push_back(in_integer_marker_section);
        column_has_lower.push_back(false);
        column_has_upper.push_back(false);
        return column_id;
    }

    void parse_rows_line(const std::vector<std::string>& tokens)
    {
        if (tokens.size() < 2) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "malformed ROWS line.");
        }
        char row_type = tokens[0][0];
        const std::string& row_name = tokens[1];

        if (row_type == 'N') {
            if (!has_objective_row) {
                has_objective_row = true;
                row_ids[row_name] = -1;
            } else {
                row_ids[row_name] = -2;
            }
            return;
        }

        if (row_type != 'L' && row_type != 'G' && row_type != 'E') {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "unknown row type \"" + tokens[0] + "\" for row \"" + row_name + "\".");
        }

        int row_id = model.number_of_constraints();
        row_ids[row_name] = row_id;
        model.constraints_names.push_back(row_name);
        row_types.push_back(row_type);
        row_has_rhs.push_back(false);
        row_has_range.push_back(false);
        row_entries.emplace_back();
        if (row_type == 'L') {
            model.constraints_lower_bounds.push_back(-inf);
            model.constraints_upper_bounds.push_back(0.0);
        } else if (row_type == 'G') {
            model.constraints_lower_bounds.push_back(0.0);
            model.constraints_upper_bounds.push_back(inf);
        } else {
            model.constraints_lower_bounds.push_back(0.0);
            model.constraints_upper_bounds.push_back(0.0);
        }
    }

    void parse_columns_line(const std::vector<std::string>& tokens)
    {
        if (tokens.size() >= 2 && tokens[1] == "'MARKER'") {
            if (tokens.size() >= 3) {
                if (tokens[2] == "'INTORG'") {
                    in_integer_marker_section = true;
                } else if (tokens[2] == "'INTEND'") {
                    in_integer_marker_section = false;
                } else {
                    throw std::invalid_argument(
                            FUNC_SIGNATURE + ": "
                            "unknown marker \"" + tokens[2] + "\" in COLUMNS section.");
                }
            }
            return;
        }

        if (tokens.size() < 3 || tokens.size() % 2 == 0) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "malformed COLUMNS line.");
        }

        int column_id = get_or_create_column(tokens[0]);
        for (size_t token_id = 1; token_id + 1 < tokens.size(); token_id += 2) {
            const std::string& row_name = tokens[token_id];
            double value = parse_number(tokens[token_id + 1]);

            auto row_it = row_ids.find(row_name);
            if (row_it == row_ids.end() || value == 0.0)
                continue;
            int row_id = row_it->second;

            if (row_id == -1) {
                if (!column_has_cost[column_id]) {
                    model.objective_coefficients[column_id] = value;
                    column_has_cost[column_id] = true;
                }
            } else if (row_id >= 0) {
                row_entries[row_id].push_back(std::make_pair(column_id, value));
            }
            // row_id == -2: free row, dropped.
        }
    }

    void parse_rhs_line(const std::vector<std::string>& tokens)
    {
        // tokens[0] is the RHS vector name and is ignored.
        if (tokens.size() < 3 || tokens.size() % 2 == 0) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "malformed RHS line.");
        }
        for (size_t token_id = 1; token_id + 1 < tokens.size(); token_id += 2) {
            const std::string& row_name = tokens[token_id];
            double value = parse_number(tokens[token_id + 1]);

            auto row_it = row_ids.find(row_name);
            if (row_it == row_ids.end())
                continue;
            int row_id = row_it->second;
            // row_id == -1 (objective offset) and -2 (free row) are dropped:
            // MathOptModel has no field for an objective constant term.
            if (row_id < 0 || row_has_rhs[row_id])
                continue;
            row_has_rhs[row_id] = true;

            char row_type = row_types[row_id];
            if (row_type == 'E' || row_type == 'L')
                model.constraints_upper_bounds[row_id] = value;
            if (row_type == 'E' || row_type == 'G')
                model.constraints_lower_bounds[row_id] = value;
        }
    }

    void parse_ranges_line(const std::vector<std::string>& tokens)
    {
        // tokens[0] is the RANGES vector name and is ignored.
        if (tokens.size() < 3 || tokens.size() % 2 == 0) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "malformed RANGES line.");
        }
        for (size_t token_id = 1; token_id + 1 < tokens.size(); token_id += 2) {
            const std::string& row_name = tokens[token_id];
            double value = parse_number(tokens[token_id + 1]);

            auto row_it = row_ids.find(row_name);
            if (row_it == row_ids.end())
                continue;
            int row_id = row_it->second;
            if (row_id < 0 || row_has_range[row_id])
                continue;
            row_has_range[row_id] = true;

            char row_type = row_types[row_id];
            double r = std::abs(value);
            if ((row_type == 'E' && value < 0) || row_type == 'L') {
                model.constraints_lower_bounds[row_id] = model.constraints_upper_bounds[row_id] - r;
            } else if ((row_type == 'E' && value >= 0) || row_type == 'G') {
                model.constraints_upper_bounds[row_id] = model.constraints_lower_bounds[row_id] + r;
            }
        }
    }

    void parse_bounds_line(const std::vector<std::string>& tokens)
    {
        // tokens[0] is the bound type, tokens[1] the bound vector name
        // (ignored), tokens[2] the column name, tokens[3] (if present) the
        // value.
        if (tokens.size() < 3) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "malformed BOUNDS line.");
        }
        const std::string& bound_type = tokens[0];
        int column_id = get_or_create_column(tokens[2]);
        column_still_default_binary[column_id] = false;

        if (bound_type == "MI") {
            if (!column_has_lower[column_id]) {
                model.variables_lower_bounds[column_id] = -inf;
                column_has_lower[column_id] = true;
            }
            return;
        }
        if (bound_type == "PL") {
            if (!column_has_upper[column_id]) {
                model.variables_upper_bounds[column_id] = inf;
                column_has_upper[column_id] = true;
            }
            return;
        }
        if (bound_type == "FR") {
            if (!column_has_lower[column_id]) {
                model.variables_lower_bounds[column_id] = -inf;
                column_has_lower[column_id] = true;
            }
            if (!column_has_upper[column_id]) {
                model.variables_upper_bounds[column_id] = inf;
                column_has_upper[column_id] = true;
            }
            return;
        }
        if (bound_type == "BV") {
            column_is_integer[column_id] = true;
            if (!column_has_lower[column_id]) {
                model.variables_lower_bounds[column_id] = 0.0;
                column_has_lower[column_id] = true;
            }
            if (!column_has_upper[column_id]) {
                model.variables_upper_bounds[column_id] = 1.0;
                column_has_upper[column_id] = true;
            }
            return;
        }

        // UP, LO, FX, LI, UI need a value.
        if (tokens.size() < 4) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "no value given for \"" + bound_type + "\" bound of column \"" + tokens[2] + "\".");
        }
        double value = parse_number(tokens[3]);

        if (bound_type == "UP") {
            if (!column_has_upper[column_id]) {
                model.variables_upper_bounds[column_id] = value;
                column_has_upper[column_id] = true;
            }
        } else if (bound_type == "LO") {
            if (!column_has_lower[column_id]) {
                model.variables_lower_bounds[column_id] = value;
                column_has_lower[column_id] = true;
            }
        } else if (bound_type == "FX") {
            if (!column_has_lower[column_id]) {
                model.variables_lower_bounds[column_id] = value;
                column_has_lower[column_id] = true;
            }
            if (!column_has_upper[column_id]) {
                model.variables_upper_bounds[column_id] = value;
                column_has_upper[column_id] = true;
            }
        } else if (bound_type == "LI") {
            column_is_integer[column_id] = true;
            if (!column_has_lower[column_id]) {
                model.variables_lower_bounds[column_id] = value;
                column_has_lower[column_id] = true;
            }
        } else if (bound_type == "UI") {
            column_is_integer[column_id] = true;
            if (!column_has_upper[column_id]) {
                model.variables_upper_bounds[column_id] = value;
                column_has_upper[column_id] = true;
            }
        } else {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "unsupported BOUNDS type \"" + bound_type + "\".");
        }
    }

    MathOptModel finalize()
    {
        // Integer columns that never received a BOUNDS entry default to
        // [0, 1] (classic MPS convention).
        for (int column_id = 0; column_id < model.number_of_variables(); ++column_id) {
            if (column_still_default_binary[column_id]) {
                model.variables_lower_bounds[column_id] = 0.0;
                model.variables_upper_bounds[column_id] = 1.0;
            }
        }

        // Classify variable types: an integer column with final bounds
        // [0, 1] is Binary, otherwise Integer.
        model.variables_types.resize(model.number_of_variables());
        for (int column_id = 0; column_id < model.number_of_variables(); ++column_id) {
            if (!column_is_integer[column_id]) {
                model.variables_types[column_id] = VariableType::Continuous;
            } else if (model.variables_lower_bounds[column_id] == 0.0
                    && model.variables_upper_bounds[column_id] == 1.0) {
                model.variables_types[column_id] = VariableType::Binary;
            } else {
                model.variables_types[column_id] = VariableType::Integer;
            }
        }

        // Flatten the per-row nonzero buffers (filled in column order by
        // COLUMNS) into the CSR layout MathOptModel expects (grouped by
        // row).
        model.constraints_starts.resize(model.number_of_constraints());
        for (int row_id = 0; row_id < model.number_of_constraints(); ++row_id) {
            model.constraints_starts[row_id] = model.number_of_elements();
            for (const auto& entry: row_entries[row_id]) {
                model.elements_variables.push_back(entry.first);
                model.elements_coefficients.push_back(entry.second);
            }
        }

        return std::move(model);
    }
};

}

MathOptModel mathoptsolverscmake::read_mps(const std::string& file_path)
{
    std::ifstream file(file_path);
    if (!file.good()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "unable to open file \"" + file_path + "\".");
    }
    return read_mps(file);
}

MathOptModel mathoptsolverscmake::read_mps(std::istream& in)
{
    MpsParseState state;
    state.model.objective_direction = ObjectiveDirection::Minimize;

    Section section = Section::None;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        if (line.empty() || line[0] == '*')
            continue;

        std::vector<std::string> tokens = split_line(line);
        if (tokens.empty())
            continue;

        std::string first_upper = to_upper(tokens[0]);
        if (first_upper == "ENDATA")
            break;
        if (is_unsupported_section_keyword(first_upper)) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "unsupported MPS section \"" + tokens[0] + "\".");
        }
        Section keyword = section_from_keyword(first_upper);
        // A data line's leading field (e.g. the RHS/RANGES/BOUNDS vector
        // name, conventionally "RHS") commonly coincides with a section
        // keyword; only treat it as a new section header when it actually
        // changes the section.
        if (keyword != Section::None && keyword != section) {
            section = keyword;
            // Gurobi-style "OBJSENSE MAX" on the header line itself.
            if (keyword == Section::Objsense && tokens.size() >= 2) {
                std::string sense = to_upper(tokens[1]).substr(0, 3);
                if (sense == "MAX")
                    state.model.objective_direction = ObjectiveDirection::Maximize;
                else if (sense == "MIN")
                    state.model.objective_direction = ObjectiveDirection::Minimize;
            }
            continue;
        }

        switch (section) {
        case Section::None:
        case Section::Name: {
            // Nothing to do: the model name is not stored on MathOptModel.
            break;
        } case Section::Objsense: {
            std::string sense = to_upper(tokens[0]).substr(0, 3);
            if (sense == "MAX")
                state.model.objective_direction = ObjectiveDirection::Maximize;
            else if (sense == "MIN")
                state.model.objective_direction = ObjectiveDirection::Minimize;
            break;
        } case Section::Rows: {
            state.parse_rows_line(tokens);
            break;
        } case Section::Columns: {
            state.parse_columns_line(tokens);
            break;
        } case Section::Rhs: {
            state.parse_rhs_line(tokens);
            break;
        } case Section::Ranges: {
            state.parse_ranges_line(tokens);
            break;
        } case Section::Bounds: {
            state.parse_bounds_line(tokens);
            break;
        }
        }
    }

    if (!state.has_objective_row) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "no objective ('N') row found.");
    }

    return state.finalize();
}

void mathoptsolverscmake::write_mps(
        const MathOptModel& model,
        const std::string& file_path)
{
    std::ofstream file(file_path);
    if (!file.good()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "unable to open file \"" + file_path + "\" for writing.");
    }
    write_mps(model, file);
}

void mathoptsolverscmake::write_mps(
        const MathOptModel& model,
        std::ostream& out)
{
    if (!model.is_milp()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "model is not an MILP; MPS format only supports linear / "
                "mixed-integer-linear models.");
    }

    const double inf = std::numeric_limits<double>::infinity();
    int number_of_variables = model.number_of_variables();
    int number_of_constraints = model.number_of_constraints();

    std::vector<std::string> variable_names(number_of_variables);
    for (int variable_id = 0; variable_id < number_of_variables; ++variable_id) {
        variable_names[variable_id] = model.variable_name(variable_id);
        check_mps_name(variable_names[variable_id]);
    }
    std::vector<std::string> constraint_names(number_of_constraints);
    std::set<std::string> used_row_names;
    for (int constraint_id = 0; constraint_id < number_of_constraints; ++constraint_id) {
        constraint_names[constraint_id] = model.constraint_name(constraint_id);
        check_mps_name(constraint_names[constraint_id]);
        used_row_names.insert(constraint_names[constraint_id]);
    }
    // Avoid a name clash with a constraint that happens to be named "COST".
    std::string objective_row_name = "COST";
    while (used_row_names.count(objective_row_name))
        objective_row_name += "_";

    out << std::setprecision(17);

    out << "NAME\n";

    if (model.objective_direction == ObjectiveDirection::Maximize) {
        out << "OBJSENSE\n"
            << " MAX\n";
    }

    // ROWS.
    out << "ROWS\n"
        << " N  " << objective_row_name << "\n";
    std::vector<char> row_types(number_of_constraints);
    for (int constraint_id = 0; constraint_id < number_of_constraints; ++constraint_id) {
        char row_type = 'L';
        switch (model.constraint_sense(constraint_id)) {
        case ConstraintSense::LessThanOrEqualTo: row_type = 'L'; break;
        case ConstraintSense::GreaterThanOrEqualTo: row_type = 'G'; break;
        case ConstraintSense::Equality: row_type = 'E'; break;
        // A range constraint is written as an 'L' row plus a RANGES entry.
        case ConstraintSense::Range: row_type = 'L'; break;
        // A free constraint has no MPS row type of its own; write it as a
        // (non-objective) free row, matching how solvers such as HiGHS
        // themselves emit one. It is dropped again on read-back, since
        // free rows other than the objective carry no information.
        case ConstraintSense::Free: row_type = 'N'; break;
        }
        row_types[constraint_id] = row_type;
        out << " " << row_type << "  " << constraint_names[constraint_id] << "\n";
    }

    // COLUMNS: invert the row-wise matrix into per-column nonzero lists.
    std::vector<std::vector<std::pair<int, double>>> column_entries(number_of_variables);
    for (int constraint_id = 0; constraint_id < number_of_constraints; ++constraint_id) {
        for (int element_id = model.constraints_starts[constraint_id];
                element_id < model.constraint_end(constraint_id);
                ++element_id) {
            double coefficient = model.elements_coefficients[element_id];
            if (coefficient == 0.0)
                continue;
            column_entries[model.elements_variables[element_id]].push_back(
                    std::make_pair(constraint_id, coefficient));
        }
    }

    out << "COLUMNS\n";
    bool in_integer_marker_section = false;
    for (int variable_id = 0; variable_id < number_of_variables; ++variable_id) {
        bool is_integer = (model.variables_types[variable_id] != VariableType::Continuous);
        if (is_integer != in_integer_marker_section) {
            out << "    MARKER    'MARKER'                 '"
                << (is_integer? "INTORG": "INTEND") << "'\n";
            in_integer_marker_section = is_integer;
        }

        // Every variable must appear at least once in COLUMNS, even with
        // an all-zero column, or read_mps() would not create it at all.
        std::vector<std::pair<std::string, double>> pairs;
        double objective_coefficient = model.objective_coefficients[variable_id];
        if (objective_coefficient != 0.0)
            pairs.push_back(std::make_pair(objective_row_name, objective_coefficient));
        for (const auto& entry: column_entries[variable_id])
            pairs.push_back(std::make_pair(constraint_names[entry.first], entry.second));
        if (pairs.empty())
            pairs.push_back(std::make_pair(objective_row_name, 0.0));

        for (size_t pair_id = 0; pair_id < pairs.size(); pair_id += 2) {
            out << "    " << variable_names[variable_id]
                << "  " << pairs[pair_id].first << "  " << pairs[pair_id].second;
            if (pair_id + 1 < pairs.size())
                out << "  " << pairs[pair_id + 1].first << "  " << pairs[pair_id + 1].second;
            out << "\n";
        }
    }
    if (in_integer_marker_section)
        out << "    MARKER    'MARKER'                 'INTEND'\n";

    // RHS.
    out << "RHS\n";
    std::vector<std::pair<std::string, double>> rhs_pairs;
    for (int constraint_id = 0; constraint_id < number_of_constraints; ++constraint_id) {
        double lower = model.constraints_lower_bounds[constraint_id];
        double upper = model.constraints_upper_bounds[constraint_id];
        switch (row_types[constraint_id]) {
        case 'L': rhs_pairs.push_back(std::make_pair(constraint_names[constraint_id], upper)); break;
        case 'G': rhs_pairs.push_back(std::make_pair(constraint_names[constraint_id], lower)); break;
        case 'E': rhs_pairs.push_back(std::make_pair(constraint_names[constraint_id], lower)); break;
        default: break; // 'N': free row, no RHS.
        }
    }
    for (size_t pair_id = 0; pair_id < rhs_pairs.size(); pair_id += 2) {
        out << "    RHS       " << rhs_pairs[pair_id].first << "  " << rhs_pairs[pair_id].second;
        if (pair_id + 1 < rhs_pairs.size())
            out << "  " << rhs_pairs[pair_id + 1].first << "  " << rhs_pairs[pair_id + 1].second;
        out << "\n";
    }

    // RANGES: only range constraints need one.
    std::vector<std::pair<std::string, double>> range_pairs;
    for (int constraint_id = 0; constraint_id < number_of_constraints; ++constraint_id) {
        if (model.constraint_sense(constraint_id) != ConstraintSense::Range)
            continue;
        double range = model.constraints_upper_bounds[constraint_id] - model.constraints_lower_bounds[constraint_id];
        range_pairs.push_back(std::make_pair(constraint_names[constraint_id], range));
    }
    if (!range_pairs.empty()) {
        out << "RANGES\n";
        for (size_t pair_id = 0; pair_id < range_pairs.size(); pair_id += 2) {
            out << "    RNG       " << range_pairs[pair_id].first << "  " << range_pairs[pair_id].second;
            if (pair_id + 1 < range_pairs.size())
                out << "  " << range_pairs[pair_id + 1].first << "  " << range_pairs[pair_id + 1].second;
            out << "\n";
        }
    }

    // BOUNDS: only variables that deviate from the applicable MPS default
    // (continuous [0, +inf), or [0, 1] for a variable wrapped in an
    // INTORG/INTEND marker) need an entry.
    std::ostringstream bounds_stream;
    bounds_stream << std::setprecision(17);
    bool has_bounds = false;
    for (int variable_id = 0; variable_id < number_of_variables; ++variable_id) {
        double lower = model.variables_lower_bounds[variable_id];
        double upper = model.variables_upper_bounds[variable_id];
        bool is_integer = (model.variables_types[variable_id] != VariableType::Continuous);
        const std::string& name = variable_names[variable_id];

        if (!is_integer && lower == 0.0 && upper == inf)
            continue;
        has_bounds = true;

        if (is_integer && lower == 0.0 && upper == 1.0) {
            bounds_stream << " BV BND       " << name << "\n";
        } else if (lower == -inf && upper == inf) {
            bounds_stream << " FR BND       " << name << "\n";
        } else if (lower == upper) {
            bounds_stream << " FX BND       " << name << "  " << lower << "\n";
        } else {
            if (lower == -inf)
                bounds_stream << " MI BND       " << name << "\n";
            else if (lower != 0.0)
                bounds_stream << " LO BND       " << name << "  " << lower << "\n";
            if (upper == inf)
                bounds_stream << " PL BND       " << name << "\n";
            else
                bounds_stream << " UP BND       " << name << "  " << upper << "\n";
        }
    }
    if (has_bounds) {
        out << "BOUNDS\n"
            << bounds_stream.str();
    }

    out << "ENDATA\n";
}

#include "mathoptsolverscmake/mathopt_nl.hpp"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

using namespace mathoptsolverscmake;

namespace
{

const double inf = std::numeric_limits<double>::infinity();

/*
 * Opcode table: translates between .nl's numeric "o<opcode>" prefix-tree
 * opcodes and MathOptModel's own char-opcode tree (see the "Nonlinear
 * structures" comment on MathOptModel in mathopt.hpp), for exactly the
 * operators MathOptModel supports. Arity is 1 (unary), 2 (binary) or -1
 * (n-ary, with an explicit child count following the opcode -- only used
 * by OPSUMLIST, .nl's n-ary '+').
 */

struct OpcodeInfo
{
    int opcode;
    char op;
    int arity;
};

const OpcodeInfo kOpcodeTable[] = {
    {0, '+', 2},    // OPPLUS
    {1, '-', 2},    // OPMINUS
    {2, '*', 2},    // OPMULT
    {3, '/', 2},    // OPDIV
    {5, 'p', 2},    // OPPOW
    {16, 'n', 1},   // OPUMINUS
    {38, 't', 1},   // OP_tan
    {39, 'q', 1},   // OP_sqrt
    {41, 's', 1},   // OP_sin
    {43, 'l', 1},   // OP_log
    {44, 'e', 1},   // OP_exp
    {46, 'c', 1},   // OP_cos
    {54, '+', -1},  // OPSUMLIST
};

bool opcode_to_char(int opcode, char& op, int& arity)
{
    for (const OpcodeInfo& entry: kOpcodeTable) {
        if (entry.opcode == opcode) {
            op = entry.op;
            arity = entry.arity;
            return true;
        }
    }
    return false;
}

int char_to_opcode(char op)
{
    switch (op) {
    case '-': return 1;
    case '*': return 2;
    case '/': return 3;
    case 'p': return 5;
    case 'n': return 16;
    case 't': return 38;
    case 'q': return 39;
    case 's': return 41;
    case 'l': return 43;
    case 'e': return 44;
    case 'c': return 46;
    default:
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "no .nl opcode for nonlinear-expression operator '" + std::string(1, op) + "'.");
    }
}

/*
 * Header parsing: the first 10 lines of a .nl file, each a fixed but
 * sometimes-truncated (trailing fields default to 0) whitespace-separated
 * list of integers, following jac0dim.c's Sscanf(" %d %d ...", ...)
 * pattern. See mathopt_nl.hpp's doc comment for which fields this reader
 * requires to be 0 (unsupported .nl features).
 */

struct NlHeader
{
    int n_var = 0, n_con = 0, n_obj = 0, n_lcon = 0;
    int nlc = 0, nlo = 0, n_cc = 0, nlcc = 0, ndcc = 0, nzlb = 0;
    int nlnc = 0, lnc = 0;
    int nlvc = 0, nlvo = 0, nlvb = 0;
    int nwv = 0, nfunc = 0;
    int nbv = 0, niv = 0, nlvbi = 0, nlvci = 0, nlvoi = 0;
    int comb = 0, comc = 0, como = 0, comc1 = 0, como1 = 0;
};

std::string strip_comment(const std::string& line)
{
    size_t hash = line.find('#');
    return hash == std::string::npos? line: line.substr(0, hash);
}

/** Splits a (comment-stripped) line into as many integers as are present. */
std::vector<long> split_ints(const std::string& line)
{
    std::vector<long> result;
    std::istringstream iss(line);
    long value;
    while (iss >> value)
        result.push_back(value);
    return result;
}

long field(const std::vector<long>& fields, size_t index, long default_value = 0)
{
    return index < fields.size()? fields[index]: default_value;
}

void require_zero(long value, const std::string& what)
{
    if (value != 0) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "unsupported .nl file: " + what + " is not supported (got " + std::to_string(value) + ").");
    }
}

NlHeader read_nl_header(std::istream& in, std::string& first_data_line)
{
    std::string line;
    if (!std::getline(in, line)) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "empty .nl input.");
    }
    if (line.empty() || (line[0] != 'g' && line[0] != 'G')) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "unsupported .nl format: only the ASCII \"g\" format is supported "
                "(binary \"b\"/\"h\" formats are not).");
    }
    // Skip the trailing "k opt1 .. optk" AMPL-internal option list, if any:
    // 'k' (how many options follow) is glued right after 'g', the options
    // themselves are whitespace-separated on the same line.

    NlHeader header;

    std::getline(in, line);
    std::vector<long> l2 = split_ints(strip_comment(line));
    if (l2.size() < 3) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "malformed .nl header (variable/constraint/objective counts).");
    }
    header.n_var = (int)field(l2, 0);
    header.n_con = (int)field(l2, 1);
    header.n_obj = (int)field(l2, 2);
    header.n_lcon = (int)field(l2, 5);

    std::getline(in, line);
    std::vector<long> l3 = split_ints(strip_comment(line));
    header.nlc = (int)field(l3, 0);
    header.nlo = (int)field(l3, 1);
    header.n_cc = (int)field(l3, 2);
    header.nlcc = (int)field(l3, 3);
    header.ndcc = (int)field(l3, 4);
    header.nzlb = (int)field(l3, 5);

    std::getline(in, line);
    std::vector<long> l4 = split_ints(strip_comment(line));
    header.nlnc = (int)field(l4, 0);
    header.lnc = (int)field(l4, 1);

    std::getline(in, line);
    std::vector<long> l5 = split_ints(strip_comment(line));
    header.nlvc = (int)field(l5, 0);
    header.nlvo = (int)field(l5, 1);
    header.nlvb = (int)field(l5, 2, -1);

    std::getline(in, line);
    std::vector<long> l6 = split_ints(strip_comment(line));
    header.nwv = (int)field(l6, 0);
    header.nfunc = (int)field(l6, 1);

    if (header.nlvb < 0) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "unsupported .nl file: pre-1993 header format (missing nlvb).");
    }
    std::getline(in, line);
    std::vector<long> l7 = split_ints(strip_comment(line));
    if (l7.size() != 5) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "malformed .nl header (discrete variable counts).");
    }
    header.nbv = (int)l7[0];
    header.niv = (int)l7[1];
    header.nlvbi = (int)l7[2];
    header.nlvci = (int)l7[3];
    header.nlvoi = (int)l7[4];

    std::getline(in, line); // nzc nzo: nonzero counts, informational only.
    std::getline(in, line); // maxrownamelen maxcolnamelen: names aren't read.

    std::getline(in, line);
    std::vector<long> l10 = split_ints(strip_comment(line));
    header.comb = (int)field(l10, 0);
    header.comc = (int)field(l10, 1);
    header.como = (int)field(l10, 2);
    header.comc1 = (int)field(l10, 3);
    header.como1 = (int)field(l10, 4);

    require_zero(header.n_lcon, "logical constraints");
    require_zero(header.n_cc + header.nlcc, "complementarity constraints");
    require_zero(header.nlnc + header.lnc, "network constraints");
    require_zero(header.nwv, "network (\"arc\") variables");
    require_zero(header.nfunc, "user-defined functions");
    require_zero(header.comb + header.comc + header.como + header.comc1 + header.como1,
            "common subexpressions");
    if (header.n_obj > 1) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "unsupported .nl file: only a single objective is supported "
                "(got " + std::to_string(header.n_obj) + ").");
    }

    first_data_line.clear();
    return header;
}

/*
 * Body parsing: everything after the header is tokenized as a flat,
 * whitespace-separated stream (newlines carry no meaning beyond
 * separating tokens here), since every segment either gives its own
 * element count up front or has a statically-known one. A marker token
 * such as "C3", "o16" or "x7" is a single letter immediately followed
 * (no space) by an integer; segments below read that glued integer via
 * expect_marker().
 */

class NlBodyParser
{

public:

    explicit NlBodyParser(const std::string& content)
    {
        std::string stripped;
        stripped.reserve(content.size());
        std::istringstream lines(content);
        std::string line;
        while (std::getline(lines, line)) {
            stripped += strip_comment(line);
            stripped += '\n';
        }
        std::istringstream iss(stripped);
        std::string token;
        while (iss >> token)
            tokens_.push_back(token);
    }

    bool at_end() const { return pos_ >= tokens_.size(); }

    const std::string& peek() const
    {
        if (at_end()) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "unexpected end of .nl input.");
        }
        return tokens_[pos_];
    }

    std::string advance()
    {
        const std::string& token = peek();
        ++pos_;
        return token;
    }

    long advance_int()
    {
        const std::string& token = advance();
        return std::strtol(token.c_str(), nullptr, 10);
    }

    double advance_double()
    {
        const std::string& token = advance();
        return std::strtod(token.c_str(), nullptr);
    }

    /** Reads a token required to start with 'marker', returning the glued integer that follows. */
    long expect_marker(char marker)
    {
        std::string token = advance();
        if (token.empty() || token[0] != marker || token.size() < 2) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "expected a \"" + std::string(1, marker) + "<number>\" token, got \"" + token + "\".");
        }
        return std::strtol(token.c_str() + 1, nullptr, 10);
    }

    /**
     * Reads one expression node (and, recursively, its children) in
     * .nl's prefix notation directly into MathOptModel's own DFS
     * pre-order parallel arrays.
     */
    void read_expression(
            std::vector<char>& operators,
            std::vector<double>& values,
            std::vector<int>& variables,
            std::vector<int>& number_of_children)
    {
        std::string token = advance();
        if (token.empty()) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "expected an expression token.");
        }
        char kind = token[0];
        if (kind == 'n') {
            operators.push_back('k');
            values.push_back(std::strtod(token.c_str() + 1, nullptr));
            variables.push_back(-1);
            number_of_children.push_back(0);
            return;
        }
        if (kind == 'v') {
            operators.push_back('v');
            values.push_back(0.0);
            variables.push_back((int)std::strtol(token.c_str() + 1, nullptr, 10));
            number_of_children.push_back(0);
            return;
        }
        if (kind != 'o') {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "expected an 'n', 'v' or 'o' expression token, got \"" + token + "\".");
        }
        int opcode = (int)std::strtol(token.c_str() + 1, nullptr, 10);
        char op;
        int arity;
        if (!opcode_to_char(opcode, op, arity)) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "unsupported .nl opcode o" + std::to_string(opcode) + ".");
        }
        if (arity == 1) {
            operators.push_back(op);
            values.push_back(0.0);
            variables.push_back(-1);
            number_of_children.push_back(1);
            read_expression(operators, values, variables, number_of_children);
        } else if (arity == 2) {
            operators.push_back(op);
            values.push_back(0.0);
            variables.push_back(-1);
            number_of_children.push_back(2);
            read_expression(operators, values, variables, number_of_children);
            read_expression(operators, values, variables, number_of_children);
        } else {
            int count = (int)advance_int();
            operators.push_back(op);
            values.push_back(0.0);
            variables.push_back(-1);
            number_of_children.push_back(count);
            for (int i = 0; i < count; ++i)
                read_expression(operators, values, variables, number_of_children);
        }
    }

private:

    std::vector<std::string> tokens_;
    size_t pos_ = 0;

};

void read_bound(NlBodyParser& parser, double& lower, double& upper)
{
    int type = (int)parser.advance_int();
    switch (type) {
    case 0: lower = parser.advance_double(); upper = parser.advance_double(); break;
    case 1: lower = -inf; upper = parser.advance_double(); break;
    case 2: lower = parser.advance_double(); upper = inf; break;
    case 3: lower = -inf; upper = inf; break;
    case 4: lower = parser.advance_double(); upper = lower; break;
    default:
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "unsupported bound type " + std::to_string(type) + " "
                "(complementarity bounds are not supported).");
    }
}

/** True if the just-parsed expression is exactly the trivial constant 0 (.nl's placeholder for "no nonlinear part"). */
bool is_trivial_zero(
        const std::vector<char>& operators,
        const std::vector<double>& values,
        size_t start)
{
    return operators.size() == start + 1 && operators[start] == 'k' && values[start] == 0.0;
}

/*
 * Writer helpers.
 */

struct ExpressionArrays
{
    const std::vector<char>& operators;
    const std::vector<double>& values;
    const std::vector<int>& variables;
    const std::vector<int>& number_of_children;
};

size_t write_subtree(
        std::ostream& out,
        const ExpressionArrays& arrays,
        size_t index,
        const std::vector<int>& variable_remap);

size_t write_nary_multiply(
        std::ostream& out,
        const ExpressionArrays& arrays,
        size_t index,
        int remaining,
        const std::vector<int>& variable_remap)
{
    if (remaining == 1)
        return write_subtree(out, arrays, index, variable_remap);
    out << "o2\n";
    size_t next = write_subtree(out, arrays, index, variable_remap);
    return write_nary_multiply(out, arrays, next, remaining - 1, variable_remap);
}

size_t write_subtree(
        std::ostream& out,
        const ExpressionArrays& arrays,
        size_t index,
        const std::vector<int>& variable_remap)
{
    char op = arrays.operators[index];
    if (op == 'k') {
        out << "n" << arrays.values[index] << "\n";
        return index + 1;
    }
    if (op == 'v') {
        out << "v" << variable_remap[arrays.variables[index]] << "\n";
        return index + 1;
    }
    if (op == '+' && arrays.number_of_children[index] > 2) {
        int n = arrays.number_of_children[index];
        out << "o54\n" << n << "\n";
        size_t next = index + 1;
        for (int i = 0; i < n; ++i)
            next = write_subtree(out, arrays, next, variable_remap);
        return next;
    }
    if (op == '*' && arrays.number_of_children[index] > 2) {
        return write_nary_multiply(out, arrays, index + 1, arrays.number_of_children[index], variable_remap);
    }
    if (arrays.number_of_children[index] == 1) {
        out << "o" << char_to_opcode(op) << "\n";
        return write_subtree(out, arrays, index + 1, variable_remap);
    }
    // Binary: '+', '-', '*', '/', 'p'.
    out << "o" << (op == '+'? 0: char_to_opcode(op)) << "\n";
    size_t next = write_subtree(out, arrays, index + 1, variable_remap);
    return write_subtree(out, arrays, next, variable_remap);
}

/** Writes the combined (nonlinear tree + folded quadratic terms) expression for one row. */
void write_combined_expression(
        std::ostream& out,
        const ExpressionArrays& nonlinear,
        size_t nonlinear_start,
        size_t nonlinear_end,
        const std::vector<int>& quadratic_variables_1,
        const std::vector<int>& quadratic_variables_2,
        const std::vector<double>& quadratic_coefficients,
        size_t quadratic_start,
        size_t quadratic_end,
        const std::vector<int>& variable_remap)
{
    bool has_nonlinear = nonlinear_end > nonlinear_start;
    size_t quadratic_count = quadratic_end - quadratic_start;
    size_t term_count = (has_nonlinear? 1: 0) + quadratic_count;

    if (term_count == 0) {
        out << "n0\n";
        return;
    }
    if (term_count > 1)
        out << "o54\n" << term_count << "\n";
    if (has_nonlinear)
        write_subtree(out, nonlinear, nonlinear_start, variable_remap);
    for (size_t i = quadratic_start; i < quadratic_end; ++i) {
        out << "o2\nn" << quadratic_coefficients[i]
            << "\no2\nv" << variable_remap[quadratic_variables_1[i]]
            << "\nv" << variable_remap[quadratic_variables_2[i]] << "\n";
    }
}

void write_bound(std::ostream& out, double lower, double upper)
{
    if (lower == -inf && upper == inf) {
        out << "3\n";
    } else if (lower == upper) {
        out << "4 " << lower << "\n";
    } else if (lower == -inf) {
        out << "1 " << upper << "\n";
    } else if (upper == inf) {
        out << "2 " << lower << "\n";
    } else {
        out << "0 " << lower << " " << upper << "\n";
    }
}

}

MathOptModel mathoptsolverscmake::read_nl(const std::string& file_path)
{
    std::ifstream file(file_path);
    if (!file.good()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "unable to open file \"" + file_path + "\".");
    }
    return read_nl(file);
}

MathOptModel mathoptsolverscmake::read_nl(std::istream& in)
{
    std::string unused;
    NlHeader header = read_nl_header(in, unused);

    int n_var = header.n_var;
    int n_con = header.n_con;
    int nlv = std::max(header.nlvc, header.nlvo);
    if (header.nlvbi > header.nlvb
            || header.nlvci > header.nlvc - header.nlvb
            || header.nlvoi > nlv - header.nlvc
            || nlv + header.nbv + header.niv > n_var) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "inconsistent .nl header variable-partition counts.");
    }

    MathOptModel model;
    model.variables_lower_bounds.assign(n_var, 0.0);
    model.variables_upper_bounds.assign(n_var, inf);
    model.variables_types.assign(n_var, VariableType::Continuous);
    model.objective_coefficients.assign(n_var, 0.0);
    model.constraints_lower_bounds.assign(n_con, -inf);
    model.constraints_upper_bounds.assign(n_con, inf);
    model.constraints_starts.assign(n_con, 0);
    model.nonlinear_elements_constraints_starts.assign(n_con, 0);

    std::string body((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    NlBodyParser parser(body);

    std::vector<std::vector<std::pair<int, double>>> constraint_linear_terms(n_con);
    bool have_x0 = false;

    while (!parser.at_end()) {
        const std::string& token = parser.peek();
        char kind = token[0];
        switch (kind) {
        case 'C': {
            int constraint_id = (int)parser.expect_marker('C');
            if (constraint_id < 0 || constraint_id >= n_con) {
                throw std::invalid_argument(
                        FUNC_SIGNATURE + ": "
                        "constraint index " + std::to_string(constraint_id) + " out of range in 'C' segment.");
            }
            size_t start = model.nonlinear_elements_operators.size();
            model.nonlinear_elements_constraints_starts[constraint_id] = (int)start;
            parser.read_expression(
                    model.nonlinear_elements_operators,
                    model.nonlinear_elements_values,
                    model.nonlinear_elements_variables,
                    model.nonlinear_elements_number_of_children);
            if (is_trivial_zero(model.nonlinear_elements_operators, model.nonlinear_elements_values, start)) {
                model.nonlinear_elements_operators.resize(start);
                model.nonlinear_elements_values.resize(start);
                model.nonlinear_elements_variables.resize(start);
                model.nonlinear_elements_number_of_children.resize(start);
                model.nonlinear_elements_constraints_starts[constraint_id] = (int)start;
            }
            break;
        } case 'O': {
            int objective_id = (int)parser.expect_marker('O');
            if (objective_id != 0) {
                throw std::invalid_argument(
                        FUNC_SIGNATURE + ": "
                        "unsupported .nl file: objective index " + std::to_string(objective_id) + " != 0.");
            }
            long sense = parser.advance_int();
            model.objective_direction = (sense == 1)? ObjectiveDirection::Maximize: ObjectiveDirection::Minimize;
            size_t start = model.objective_nonlinear_elements_operators.size();
            parser.read_expression(
                    model.objective_nonlinear_elements_operators,
                    model.objective_nonlinear_elements_values,
                    model.objective_nonlinear_elements_variables,
                    model.objective_nonlinear_elements_number_of_children);
            if (is_trivial_zero(model.objective_nonlinear_elements_operators, model.objective_nonlinear_elements_values, start)) {
                model.objective_nonlinear_elements_operators.resize(start);
                model.objective_nonlinear_elements_values.resize(start);
                model.objective_nonlinear_elements_variables.resize(start);
                model.objective_nonlinear_elements_number_of_children.resize(start);
            }
            break;
        } case 'J': {
            int constraint_id = (int)parser.expect_marker('J');
            long count = parser.advance_int();
            if (constraint_id < 0 || constraint_id >= n_con) {
                throw std::invalid_argument(
                        FUNC_SIGNATURE + ": "
                        "constraint index " + std::to_string(constraint_id) + " out of range in 'J' segment.");
            }
            for (long i = 0; i < count; ++i) {
                int variable_id = (int)parser.advance_int();
                double coefficient = parser.advance_double();
                constraint_linear_terms[constraint_id].push_back(std::make_pair(variable_id, coefficient));
            }
            break;
        } case 'G': {
            int objective_id = (int)parser.expect_marker('G');
            (void)objective_id;
            long count = parser.advance_int();
            for (long i = 0; i < count; ++i) {
                int variable_id = (int)parser.advance_int();
                double coefficient = parser.advance_double();
                model.objective_coefficients[variable_id] += coefficient;
            }
            break;
        } case 'b': {
            parser.advance(); // marker, no glued value.
            for (int variable_id = 0; variable_id < n_var; ++variable_id)
                read_bound(parser, model.variables_lower_bounds[variable_id], model.variables_upper_bounds[variable_id]);
            break;
        } case 'r': {
            parser.advance();
            for (int constraint_id = 0; constraint_id < n_con; ++constraint_id)
                read_bound(parser, model.constraints_lower_bounds[constraint_id], model.constraints_upper_bounds[constraint_id]);
            break;
        } case 'k': case 'K': {
            long count = parser.expect_marker(kind);
            for (long i = 0; i < count; ++i)
                parser.advance();
            break;
        } case 'x': {
            long count = parser.expect_marker('x');
            if (!have_x0) {
                model.variables_initial_values.assign(n_var, 0.0);
                have_x0 = true;
            }
            for (long i = 0; i < count; ++i) {
                int variable_id = (int)parser.advance_int();
                model.variables_initial_values[variable_id] = parser.advance_double();
            }
            break;
        } case 'd': {
            long count = parser.expect_marker('d');
            for (long i = 0; i < count; ++i) {
                parser.advance();
                parser.advance();
            }
            break;
        } case 'S': {
            long count = parser.expect_marker('S');
            (void)count;
            long n = parser.advance_int();
            parser.advance(); // suffix name.
            for (long i = 0; i < n; ++i) {
                parser.advance();
                parser.advance();
            }
            break;
        } default:
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "unsupported .nl segment \"" + token + "\".");
        }
    }

    // Flatten the per-constraint linear-term buffers into CSR.
    for (int constraint_id = 0; constraint_id < n_con; ++constraint_id) {
        model.constraints_starts[constraint_id] = model.number_of_elements();
        for (const auto& term: constraint_linear_terms[constraint_id]) {
            model.elements_variables.push_back(term.first);
            model.elements_coefficients.push_back(term.second);
        }
    }

    // Classify variable types from their position in the .nl partition.
    auto mark_integer_range = [&](int start, int end) {
        for (int variable_id = start; variable_id < end; ++variable_id) {
            if (model.variables_lower_bounds[variable_id] == 0.0 && model.variables_upper_bounds[variable_id] == 1.0)
                model.variables_types[variable_id] = VariableType::Binary;
            else
                model.variables_types[variable_id] = VariableType::Integer;
        }
    };
    mark_integer_range(header.nlvb - header.nlvbi, header.nlvb);
    mark_integer_range(header.nlvc - header.nlvci, header.nlvc);
    mark_integer_range(nlv - header.nlvoi, nlv);
    for (int variable_id = nlv; variable_id < nlv + header.nbv; ++variable_id)
        model.variables_types[variable_id] = VariableType::Binary;
    for (int variable_id = nlv + header.nbv; variable_id < nlv + header.nbv + header.niv; ++variable_id)
        model.variables_types[variable_id] = VariableType::Integer;

    return model;
}

void mathoptsolverscmake::write_nl(
        const MathOptModel& model,
        const std::string& file_path)
{
    std::ofstream file(file_path);
    if (!file.good()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "unable to open file \"" + file_path + "\" for writing.");
    }
    write_nl(model, file);
}

void mathoptsolverscmake::write_nl(
        const MathOptModel& model,
        std::ostream& out)
{
    if (model.has_black_box()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "model has a black-box objective or constraint function, "
                "which has no .nl representation.");
    }

    int n_var = model.number_of_variables();
    int n_con = model.number_of_constraints();
    // has_nonlinear()/has_quadratic() are true if EITHER the objective or
    // some constraint has that kind of content; the constraint-side
    // parallel arrays (nonlinear_elements_constraints_starts, ...) are
    // only populated (and safe to index) when a constraint itself does,
    // so that specific condition is tracked separately below.
    bool constraints_have_nonlinear = !model.nonlinear_elements_operators.empty();
    bool constraints_have_quadratic = !model.quadratic_elements_variables_1.empty();

    // Which variables appear nonlinearly (expression tree or quadratic
    // term) in the objective and/or in some constraint.
    std::vector<bool> in_objective(n_var, false), in_constraints(n_var, false);
    for (int variable_id: model.objective_nonlinear_elements_variables)
        if (variable_id >= 0)
            in_objective[variable_id] = true;
    for (size_t i = 0; i < model.objective_quadratic_elements_variables_1.size(); ++i) {
        in_objective[model.objective_quadratic_elements_variables_1[i]] = true;
        in_objective[model.objective_quadratic_elements_variables_2[i]] = true;
    }
    for (int constraint_id = 0; constraint_id < n_con; ++constraint_id) {
        if (constraints_have_nonlinear) {
            int start = model.nonlinear_elements_constraints_starts[constraint_id];
            int end = model.nonlinear_constraint_end(constraint_id);
            for (int e = start; e < end; ++e)
                if (model.nonlinear_elements_variables[e] >= 0)
                    in_constraints[model.nonlinear_elements_variables[e]] = true;
        }
        if (constraints_have_quadratic) {
            int start = model.quadratic_elements_constraints_starts[constraint_id];
            int end = model.quadratic_constraint_end(constraint_id);
            for (int e = start; e < end; ++e) {
                in_constraints[model.quadratic_elements_variables_1[e]] = true;
                in_constraints[model.quadratic_elements_variables_2[e]] = true;
            }
        }
    }

    // Partition variables into .nl's required order: both, constraints-only,
    // objective-only (continuous before integer/binary within each), then
    // linear (continuous, then binary, then integer).
    std::vector<int> group_both, group_con_only, group_obj_only;
    std::vector<int> group_lin_cont, group_lin_bin, group_lin_int;
    for (int pass = 0; pass < 2; ++pass) {
        bool want_continuous = (pass == 0);
        for (int v = 0; v < n_var; ++v) {
            bool is_continuous = model.variables_types[v] == VariableType::Continuous;
            if (is_continuous != want_continuous)
                continue;
            if (in_objective[v] && in_constraints[v])
                group_both.push_back(v);
            else if (in_constraints[v])
                group_con_only.push_back(v);
            else if (in_objective[v])
                group_obj_only.push_back(v);
        }
    }
    for (int v = 0; v < n_var; ++v) {
        if (in_objective[v] || in_constraints[v])
            continue;
        VariableType type = model.variables_types[v];
        if (type == VariableType::Continuous)
            group_lin_cont.push_back(v);
        else if (type == VariableType::Binary)
            group_lin_bin.push_back(v);
        else
            group_lin_int.push_back(v);
    }
    auto count_noncontinuous = [&](const std::vector<int>& group) {
        int count = 0;
        for (int v: group)
            if (model.variables_types[v] != VariableType::Continuous)
                ++count;
        return count;
    };
    int nlvbi = count_noncontinuous(group_both);
    int nlvci = count_noncontinuous(group_con_only);
    int nlvoi = count_noncontinuous(group_obj_only);
    int nlvb = (int)group_both.size();
    int nlvc = nlvb + (int)group_con_only.size();
    int nlvo = nlvc + (int)group_obj_only.size();
    int nbv = (int)group_lin_bin.size();
    int niv = (int)group_lin_int.size();

    std::vector<int> new_order;
    new_order.reserve(n_var);
    new_order.insert(new_order.end(), group_both.begin(), group_both.end());
    new_order.insert(new_order.end(), group_con_only.begin(), group_con_only.end());
    new_order.insert(new_order.end(), group_obj_only.begin(), group_obj_only.end());
    new_order.insert(new_order.end(), group_lin_cont.begin(), group_lin_cont.end());
    new_order.insert(new_order.end(), group_lin_bin.begin(), group_lin_bin.end());
    new_order.insert(new_order.end(), group_lin_int.begin(), group_lin_int.end());
    std::vector<int> old_to_new(n_var);
    for (int new_id = 0; new_id < n_var; ++new_id)
        old_to_new[new_order[new_id]] = new_id;

    // Nonlinear constraints (expression tree or quadratic term) first.
    std::vector<bool> constraint_is_nonlinear(n_con, false);
    for (int c = 0; c < n_con; ++c) {
        bool nl = constraints_have_nonlinear && model.nonlinear_constraint_end(c) > model.nonlinear_elements_constraints_starts[c];
        bool q = constraints_have_quadratic && model.quadratic_constraint_end(c) > model.quadratic_elements_constraints_starts[c];
        constraint_is_nonlinear[c] = nl || q;
    }
    std::vector<int> new_con_order;
    for (int c = 0; c < n_con; ++c)
        if (constraint_is_nonlinear[c])
            new_con_order.push_back(c);
    int nlc = (int)new_con_order.size();
    for (int c = 0; c < n_con; ++c)
        if (!constraint_is_nonlinear[c])
            new_con_order.push_back(c);

    bool objective_is_nonlinear = !model.objective_nonlinear_elements_operators.empty();
    bool objective_is_quadratic = !model.objective_quadratic_elements_variables_1.empty();
    int nlo = (objective_is_nonlinear || objective_is_quadratic)? 1: 0;

    // Linear nonzero counts and range/equality row counts, for the header.
    int nzc = 0;
    for (int c = 0; c < n_con; ++c)
        for (int e = model.constraints_starts[c]; e < model.constraint_end(c); ++e)
            if (model.elements_coefficients[e] != 0.0)
                ++nzc;
    int nzo = 0;
    for (int v = 0; v < n_var; ++v)
        if (model.objective_coefficients[v] != 0.0)
            ++nzo;
    int nranges = 0, n_eqn = 0;
    for (int c = 0; c < n_con; ++c) {
        double lower = model.constraints_lower_bounds[c];
        double upper = model.constraints_upper_bounds[c];
        if (lower == upper)
            ++n_eqn;
        else if (lower != -inf && upper != inf)
            ++nranges;
    }

    out << std::setprecision(17);

    out << "g0\n"
        << n_var << " " << n_con << " 1 " << nranges << " " << n_eqn << "\n"
        << nlc << " " << nlo << "\n"
        << "0 0\n"
        << nlvc << " " << nlvo << " " << nlvb << "\n"
        << "0 0\n"
        << nbv << " " << niv << " " << nlvbi << " " << nlvci << " " << nlvoi << "\n"
        << nzc << " " << nzo << "\n"
        << "0 0\n"
        << "0 0 0 0 0\n";

    out << "b\n";
    for (int new_id = 0; new_id < n_var; ++new_id) {
        int old_id = new_order[new_id];
        write_bound(out, model.variables_lower_bounds[old_id], model.variables_upper_bounds[old_id]);
    }

    out << "r\n";
    for (int new_c = 0; new_c < n_con; ++new_c) {
        int old_c = new_con_order[new_c];
        write_bound(out, model.constraints_lower_bounds[old_c], model.constraints_upper_bounds[old_c]);
    }

    // C segments: one per constraint, in the new (nonlinear-first) order.
    ExpressionArrays con_nl_arrays{
        model.nonlinear_elements_operators,
        model.nonlinear_elements_values,
        model.nonlinear_elements_variables,
        model.nonlinear_elements_number_of_children};
    for (int new_c = 0; new_c < n_con; ++new_c) {
        int old_c = new_con_order[new_c];
        out << "C" << new_c << "\n";
        size_t nl_start = 0, nl_end = 0;
        if (constraints_have_nonlinear) {
            nl_start = (size_t)model.nonlinear_elements_constraints_starts[old_c];
            nl_end = (size_t)model.nonlinear_constraint_end(old_c);
        }
        size_t q_start = 0, q_end = 0;
        if (constraints_have_quadratic) {
            q_start = (size_t)model.quadratic_elements_constraints_starts[old_c];
            q_end = (size_t)model.quadratic_constraint_end(old_c);
        }
        write_combined_expression(
                out, con_nl_arrays, nl_start, nl_end,
                model.quadratic_elements_variables_1, model.quadratic_elements_variables_2, model.quadratic_elements_coefficients,
                q_start, q_end, old_to_new);
    }

    // O0 segment: the objective.
    out << "O0 " << (model.objective_direction == ObjectiveDirection::Maximize? 1: 0) << "\n";
    ExpressionArrays obj_nl_arrays{
        model.objective_nonlinear_elements_operators,
        model.objective_nonlinear_elements_values,
        model.objective_nonlinear_elements_variables,
        model.objective_nonlinear_elements_number_of_children};
    write_combined_expression(
            out, obj_nl_arrays, 0, model.objective_nonlinear_elements_operators.size(),
            model.objective_quadratic_elements_variables_1, model.objective_quadratic_elements_variables_2,
            model.objective_quadratic_elements_coefficients,
            0, model.objective_quadratic_elements_variables_1.size(), old_to_new);

    if (!model.variables_initial_values.empty()) {
        out << "x" << n_var << "\n";
        for (int new_id = 0; new_id < n_var; ++new_id)
            out << new_id << " " << model.variables_initial_values[new_order[new_id]] << "\n";
    }

    // k segment: cumulative nonzero Jacobian column counts (columns 0..n_var-2).
    std::vector<int> column_nnz(n_var, 0);
    for (int c = 0; c < n_con; ++c)
        for (int e = model.constraints_starts[c]; e < model.constraint_end(c); ++e)
            if (model.elements_coefficients[e] != 0.0)
                ++column_nnz[old_to_new[model.elements_variables[e]]];
    if (n_var > 0) {
        out << "k" << (n_var - 1) << "\n";
        int cumulative = 0;
        for (int v = 0; v < n_var - 1; ++v) {
            cumulative += column_nnz[v];
            out << cumulative << "\n";
        }
    }

    // J segments: linear part of each constraint, in the new order, sorted by new variable id.
    for (int new_c = 0; new_c < n_con; ++new_c) {
        int old_c = new_con_order[new_c];
        std::vector<std::pair<int, double>> pairs;
        for (int e = model.constraints_starts[old_c]; e < model.constraint_end(old_c); ++e) {
            double coefficient = model.elements_coefficients[e];
            if (coefficient != 0.0)
                pairs.push_back(std::make_pair(old_to_new[model.elements_variables[e]], coefficient));
        }
        std::sort(pairs.begin(), pairs.end());
        out << "J" << new_c << " " << pairs.size() << "\n";
        for (const auto& p: pairs)
            out << p.first << " " << p.second << "\n";
    }

    // G0 segment: linear part of the objective, sorted by new variable id.
    {
        std::vector<std::pair<int, double>> pairs;
        for (int v = 0; v < n_var; ++v) {
            double coefficient = model.objective_coefficients[v];
            if (coefficient != 0.0)
                pairs.push_back(std::make_pair(old_to_new[v], coefficient));
        }
        std::sort(pairs.begin(), pairs.end());
        out << "G0 " << pairs.size() << "\n";
        for (const auto& p: pairs)
            out << p.first << " " << p.second << "\n";
    }
}

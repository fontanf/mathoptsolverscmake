#include "mathoptsolverscmake/mathopt_lp.hpp"

#include "number_parsing.hpp"
#include "text_output.hpp"

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>

using namespace mathoptsolverscmake;

namespace
{

const double inf = std::numeric_limits<double>::infinity();

std::string to_lower(const std::string& s)
{
    std::string result = s;
    for (char& c: result)
        c = (char)std::tolower((unsigned char)c);
    return result;
}

/**
 * Characters allowed in a name in addition to letters and digits. Doesn't
 * include anything the tokenizer treats specially elsewhere (whitespace,
 * '+-*:<>=', '[]^/' for quadratic terms, '\' for line comments
 * and the slash-star block comment marker), so there is no ambiguity in
 * allowing them here; real MPS
 * files sometimes carry names built from row/column indices that use
 * these, e.g. "x(1,2)", "A###13" or "z1&2.0".
 */
bool is_extra_identifier_char(char c)
{
    return c == '_' || c == '.' || c == '(' || c == ')'
        || c == '#' || c == '&' || c == '%' || c == '!'
        || c == '?' || c == '@' || c == '{' || c == '}'
        || c == '$' || c == ',' || c == ';' || c == '\'';
}

bool is_identifier_start(char c)
{
    // A leading '.' (or one of the other extra characters) is allowed too
    // (real MPS files may use names such as "..." or ".001"); the
    // tokenizer only reaches this check once the number rule (a '.'
    // immediately followed by a digit) has already failed to match, so
    // there is no ambiguity with numeric literals.
    return std::isalpha((unsigned char)c) || is_extra_identifier_char(c);
}

bool is_identifier_char(char c)
{
    return std::isalnum((unsigned char)c) || is_extra_identifier_char(c);
}

/*
 * Tokenizer: turns the whole file content into a flat token stream. .lp
 * expressions are free-form (a term list can span multiple lines), so
 * unlike the line-oriented MPS reader, parsing works over this single
 * stream rather than line by line.
 */

enum class TokenType { End, Identifier, Number, Plus, Minus, Star, Colon, Le, Ge, Eq, LBracket, RBracket, Caret, Slash };

struct Token
{
    TokenType type;
    std::string text;  // Identifier
    double value;       // Number
};

std::vector<Token> tokenize(const std::string& content)
{
    std::vector<Token> tokens;
    size_t pos = 0;
    size_t size = content.size();

    while (pos < size) {
        char c = content[pos];

        if (std::isspace((unsigned char)c)) {
            ++pos;
            continue;
        }
        if (c == '\\') {
            while (pos < size && content[pos] != '\n')
                ++pos;
            continue;
        }
        if (c == '/' && pos + 1 < size && content[pos + 1] == '*') {
            pos += 2;
            while (pos + 1 < size && !(content[pos] == '*' && content[pos + 1] == '/'))
                ++pos;
            pos = (pos + 1 < size) ? pos + 2 : size;
            continue;
        }
        if (c == '+') { tokens.push_back({TokenType::Plus, "", 0.0}); ++pos; continue; }
        if (c == '-') { tokens.push_back({TokenType::Minus, "", 0.0}); ++pos; continue; }
        if (c == '*') { tokens.push_back({TokenType::Star, "", 0.0}); ++pos; continue; }
        if (c == ':') { tokens.push_back({TokenType::Colon, "", 0.0}); ++pos; continue; }
        if (c == '<') {
            ++pos;
            if (pos < size && content[pos] == '=') ++pos;
            tokens.push_back({TokenType::Le, "", 0.0});
            continue;
        }
        if (c == '>') {
            ++pos;
            if (pos < size && content[pos] == '=') ++pos;
            tokens.push_back({TokenType::Ge, "", 0.0});
            continue;
        }
        if (c == '=') {
            ++pos;
            if (pos < size && content[pos] == '<') { ++pos; tokens.push_back({TokenType::Le, "", 0.0}); continue; }
            if (pos < size && content[pos] == '>') { ++pos; tokens.push_back({TokenType::Ge, "", 0.0}); continue; }
            tokens.push_back({TokenType::Eq, "", 0.0});
            continue;
        }
        if (c == '[') { tokens.push_back({TokenType::LBracket, "", 0.0}); ++pos; continue; }
        if (c == ']') { tokens.push_back({TokenType::RBracket, "", 0.0}); ++pos; continue; }
        if (c == '^') { tokens.push_back({TokenType::Caret, "", 0.0}); ++pos; continue; }
        if (c == '/') { tokens.push_back({TokenType::Slash, "", 0.0}); ++pos; continue; }
        if (std::isdigit((unsigned char)c) || (c == '.' && pos + 1 < size && std::isdigit((unsigned char)content[pos + 1]))) {
            size_t start = pos;
            while (pos < size && std::isdigit((unsigned char)content[pos]))
                ++pos;
            if (pos < size && content[pos] == '.') {
                ++pos;
                while (pos < size && std::isdigit((unsigned char)content[pos]))
                    ++pos;
            }
            if (pos < size && (content[pos] == 'e' || content[pos] == 'E')) {
                size_t exponent_start = pos + 1;
                size_t digits_start = exponent_start;
                if (digits_start < size && (content[digits_start] == '+' || content[digits_start] == '-'))
                    ++digits_start;
                if (digits_start < size && std::isdigit((unsigned char)content[digits_start])) {
                    pos = exponent_start;
                    if (content[pos] == '+' || content[pos] == '-')
                        ++pos;
                    while (pos < size && std::isdigit((unsigned char)content[pos]))
                        ++pos;
                }
            }
            std::string text = content.substr(start, pos - start);
            double value = number_parsing::c_strtod(text.c_str(), nullptr);
            tokens.push_back({TokenType::Number, text, value});
            continue;
        }
        if (is_identifier_start(c)) {
            size_t start = pos;
            while (pos < size && is_identifier_char(content[pos]))
                ++pos;
            tokens.push_back({TokenType::Identifier, content.substr(start, pos - start), 0.0});
            continue;
        }
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "unexpected character '" + std::string(1, c) + "' in .lp input.");
    }
    tokens.push_back({TokenType::End, "", 0.0});
    return tokens;
}

enum class RelOp { Le, Ge, Eq };

bool is_relop(TokenType type)
{
    return type == TokenType::Le || type == TokenType::Ge || type == TokenType::Eq;
}

/** result of "relop1 value1 [relop2 value2]" applied around an expression. */
void apply_single_relop(
        RelOp relop,
        double value,
        bool value_on_left,
        double& lower,
        double& upper)
{
    if (value_on_left) {
        // "value relop X".
        if (relop == RelOp::Le) lower = value;
        else if (relop == RelOp::Ge) upper = value;
        else { lower = value; upper = value; }
    } else {
        // "X relop value".
        if (relop == RelOp::Le) upper = value;
        else if (relop == RelOp::Ge) lower = value;
        else { lower = value; upper = value; }
    }
}

void apply_range_relop(
        RelOp relop1,
        double value1,
        RelOp relop2,
        double value2,
        double& lower,
        double& upper)
{
    if (relop1 == RelOp::Le && relop2 == RelOp::Le) {
        lower = value1;
        upper = value2;
    } else if (relop1 == RelOp::Ge && relop2 == RelOp::Ge) {
        upper = value1;
        lower = value2;
    } else {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "a two-sided bound must use the same direction on both sides "
                "(\"lo <= ... <= hi\" or \"hi >= ... >= lo\").");
    }
}

enum class LpSection { None, ObjMin, ObjMax, Constraints, Bounds, General, Binary, Semi, End, Unsupported };

struct QuadraticTerm
{
    int variable_1_id;
    int variable_2_id;
    double coefficient;
};

struct ParsedExpression
{
    std::vector<std::pair<int, double>> terms; // variable_id -> coefficient, deduplicated.
    std::vector<QuadraticTerm> quadratic_terms; // coefficient * x[variable_1_id] * x[variable_2_id].
    double constant = 0.0;
};

class LpParser
{

public:

    explicit LpParser(const std::string& content):
        tokens_(tokenize(content))
    { }

    MathOptModel parse()
    {
        parse_objective();
        parse_constraints_section();
        while (true) {
            LpSection section = peek_section_keyword();
            if (section == LpSection::Bounds) {
                consume_section_keyword(section);
                parse_bounds_section();
            } else if (section == LpSection::General) {
                consume_section_keyword(section);
                parse_typed_variable_list(VariableType::Integer);
            } else if (section == LpSection::Binary) {
                consume_section_keyword(section);
                parse_typed_variable_list(VariableType::Binary);
            } else if (section == LpSection::Semi) {
                consume_section_keyword(section);
                parse_typed_variable_list(VariableType::SemiContinuous);
            } else if (section == LpSection::End) {
                break;
            } else if (section == LpSection::Unsupported) {
                // SOS sections: only reject the section if it actually has
                // entries (some writers emit empty section headers).
                std::string keyword_text = peek().text;
                advance();
                if (peek_section_keyword() == LpSection::None && !at_end()) {
                    throw std::invalid_argument(
                            FUNC_SIGNATURE + ": "
                            "unsupported .lp section \"" + keyword_text + "\".");
                }
            } else if (at_end()) {
                break;
            } else {
                throw std::invalid_argument(
                        FUNC_SIGNATURE + ": "
                        "expected a section keyword (Bounds/General/Binary/Semi-continuous/End), found \"" +
                        (peek().type == TokenType::Identifier? peek().text: "<symbol>") + "\".");
            }
        }
        return finalize();
    }

private:

    std::vector<Token> tokens_;
    size_t pos_ = 0;

    MathOptModel model_;
    std::unordered_map<std::string, int> variable_ids_;
    std::vector<std::string> constraint_labels_; // "" means "no explicit label".

    const Token& peek(size_t offset = 0) const
    {
        size_t index = pos_ + offset;
        return tokens_[index < tokens_.size()? index: tokens_.size() - 1];
    }

    const Token& advance()
    {
        const Token& token = peek();
        if (pos_ + 1 < tokens_.size())
            ++pos_;
        return token;
    }

    bool at_end() const { return peek().type == TokenType::End; }

    bool at_identifier(size_t offset, const std::string& lowered_text) const
    {
        return peek(offset).type == TokenType::Identifier
            && to_lower(peek(offset).text) == lowered_text;
    }

    /**
     * Whether the current token could start a "label:" (an identifier, or
     * a bare number for the numeric-looking row names some MPS files use,
     * e.g. row "2" -- matching HiGHS's own .lp reader).
     */
    bool is_label_start() const
    {
        return peek().type == TokenType::Identifier || peek().type == TokenType::Number;
    }

    /**
     * Which section keyword starts here, if any (does not consume tokens).
     * A keyword followed by ':' is a label, as in "end: x >= 1" (a constraint
     * named "end"), not a section keyword.
     */
    LpSection peek_section_keyword() const
    {
        if (peek().type != TokenType::Identifier)
            return LpSection::None;
        if (peek(1).type == TokenType::Colon)
            return LpSection::None;
        std::string word = to_lower(peek().text);
        if (word == "subject" && at_identifier(1, "to")) return LpSection::Constraints;
        if (word == "such" && at_identifier(1, "that")) return LpSection::Constraints;
        if (word == "st" || word == "s.t.") return LpSection::Constraints;
        if (word == "minimize" || word == "min" || word == "minimum") return LpSection::ObjMin;
        if (word == "maximize" || word == "max" || word == "maximum") return LpSection::ObjMax;
        if (word == "bounds" || word == "bound") return LpSection::Bounds;
        if (word == "general" || word == "generals" || word == "gen"
                || word == "integer" || word == "integers") return LpSection::General;
        if (word == "binary" || word == "binaries" || word == "bin") return LpSection::Binary;
        if (word == "end") return LpSection::End;
        if (word == "semi" || word == "semis")
            return LpSection::Semi;
        if (word == "sos" || word == "sos1" || word == "sos2")
            return LpSection::Unsupported;
        return LpSection::None;
    }

    void consume_section_keyword(LpSection section)
    {
        std::string word = to_lower(peek().text);
        if (section == LpSection::Constraints && (word == "subject" || word == "such")) {
            advance();
            advance();
            return;
        }
        advance();
        // "Semi-continuous" tokenizes as "semi", '-', "continuous".
        if (section == LpSection::Semi && word == "semi"
                && peek().type == TokenType::Minus && at_identifier(1, "continuous")) {
            advance();
            advance();
        }
    }

    int get_or_create_variable(const std::string& name)
    {
        auto it = variable_ids_.find(name);
        if (it != variable_ids_.end())
            return it->second;

        int variable_id = model_.number_of_variables();
        variable_ids_[name] = variable_id;
        model_.variables_names.push_back(name);
        model_.variables_lower_bounds.push_back(0.0);
        model_.variables_upper_bounds.push_back(inf);
        model_.variables_types.push_back(VariableType::Continuous);
        model_.objective_coefficients.push_back(0.0);
        return variable_id;
    }

    double parse_signed_number()
    {
        double sign = 1.0;
        if (peek().type == TokenType::Plus) {
            advance();
        } else if (peek().type == TokenType::Minus) {
            sign = -1.0;
            advance();
        }
        if (peek().type == TokenType::Number)
            return sign * advance().value;
        if (peek().type == TokenType::Identifier) {
            std::string word = to_lower(peek().text);
            if (word == "infinity" || word == "inf") {
                advance();
                return sign * inf;
            }
        }
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "expected a number.");
    }

    RelOp parse_relop()
    {
        if (peek().type == TokenType::Le) { advance(); return RelOp::Le; }
        if (peek().type == TokenType::Ge) { advance(); return RelOp::Ge; }
        if (peek().type == TokenType::Eq) { advance(); return RelOp::Eq; }
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "expected a relational operator (<=, >= or =).");
    }

    /** Parses a sum of signed terms, stopping at a relational operator (in a
     * constraint/bound context) or at the next section keyword (in the
     * objective). */
    ParsedExpression parse_linear_expression(bool stop_at_relop)
    {
        ParsedExpression expression;
        std::unordered_map<int, size_t> term_index;
        bool first_term = true;

        while (true) {
            if (at_end())
                break;
            if (stop_at_relop) {
                if (is_relop(peek().type))
                    break;
            } else if (peek_section_keyword() != LpSection::None) {
                break;
            }

            double sign = 1.0;
            if (peek().type == TokenType::Plus) {
                advance();
            } else if (peek().type == TokenType::Minus) {
                sign = -1.0;
                advance();
            } else if (!first_term) {
                throw std::invalid_argument(
                        FUNC_SIGNATURE + ": "
                        "expected '+' or '-' between terms of a linear expression.");
            }
            first_term = false;

            if (peek().type == TokenType::LBracket) {
                parse_quadratic_block(sign, expression);
                continue;
            }

            double coefficient = sign;
            bool has_explicit_coefficient = false;
            if (peek().type == TokenType::Number) {
                coefficient = sign * advance().value;
                has_explicit_coefficient = true;
                if (peek().type == TokenType::Star)
                    advance();
            }

            // A variable name is normally an Identifier, but a variable
            // whose name looks like a number (e.g. column "0002") tokenizes
            // as a Number; the only way to tell it apart from a bare
            // constant term is that it immediately follows an explicit
            // coefficient with no operator in between (e.g. "+1 0002"),
            // which real .lp writers (e.g. HiGHS) use for exactly this case.
            bool next_is_variable =
                (peek().type == TokenType::Identifier && peek_section_keyword() == LpSection::None)
                || (has_explicit_coefficient && peek().type == TokenType::Number);
            if (next_is_variable) {
                int variable_id = get_or_create_variable(advance().text);
                auto it = term_index.find(variable_id);
                if (it == term_index.end()) {
                    term_index[variable_id] = expression.terms.size();
                    expression.terms.push_back(std::make_pair(variable_id, coefficient));
                } else {
                    expression.terms[it->second].second += coefficient;
                }
            } else if (has_explicit_coefficient) {
                expression.constant += coefficient;
            } else {
                throw std::invalid_argument(
                        FUNC_SIGNATURE + ": "
                        "expected a number or a variable in a linear expression.");
            }
        }

        return expression;
    }

    /**
     * Parses "[ terms ] [/ divisor]", the quadratic part of an expression,
     * each term being "[coefficient] x ^ 2" or "[coefficient] x * y", and
     * adds its terms, multiplied by 'sign' and divided by the divisor, to
     * 'expression'. CPLEX's format writes the objective's quadratic part as
     * "[ ... ] / 2" and the constraints' without a divisor.
     */
    void parse_quadratic_block(double sign, ParsedExpression& expression)
    {
        advance(); // '['
        std::vector<QuadraticTerm> terms;
        bool first_term = true;
        while (peek().type != TokenType::RBracket) {
            if (at_end()) {
                throw std::invalid_argument(
                        FUNC_SIGNATURE + ": "
                        "missing ']' at the end of quadratic terms.");
            }
            double term_sign = 1.0;
            if (peek().type == TokenType::Plus) {
                advance();
            } else if (peek().type == TokenType::Minus) {
                term_sign = -1.0;
                advance();
            } else if (!first_term) {
                throw std::invalid_argument(
                        FUNC_SIGNATURE + ": "
                        "expected '+' or '-' between quadratic terms.");
            }
            first_term = false;

            double coefficient = term_sign;
            if (peek().type == TokenType::Number) {
                coefficient = term_sign * advance().value;
                if (peek().type == TokenType::Star)
                    advance();
            }
            int variable_1_id = parse_quadratic_variable();
            int variable_2_id = variable_1_id;
            if (peek().type == TokenType::Caret) {
                advance();
                if (peek().type != TokenType::Number || peek().value != 2.0) {
                    throw std::invalid_argument(
                            FUNC_SIGNATURE + ": "
                            "only squares (\"^ 2\") are supported in quadratic terms.");
                }
                advance();
            } else if (peek().type == TokenType::Star) {
                advance();
                variable_2_id = parse_quadratic_variable();
            } else {
                throw std::invalid_argument(
                        FUNC_SIGNATURE + ": "
                        "expected \"^ 2\" or \"* variable\" in a quadratic term.");
            }
            terms.push_back({variable_1_id, variable_2_id, coefficient});
        }
        advance(); // ']'

        double divisor = 1.0;
        if (peek().type == TokenType::Slash) {
            advance();
            if (peek().type != TokenType::Number || peek().value == 0.0) {
                throw std::invalid_argument(
                        FUNC_SIGNATURE + ": "
                        "expected a non-zero number after '/'.");
            }
            divisor = advance().value;
        }
        for (QuadraticTerm& term: terms) {
            term.coefficient *= sign / divisor;
            expression.quadratic_terms.push_back(term);
        }
    }

    int parse_quadratic_variable()
    {
        if (peek().type != TokenType::Identifier) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "expected a variable name in a quadratic term.");
        }
        return get_or_create_variable(advance().text);
    }

    void parse_objective()
    {
        LpSection section = peek_section_keyword();
        if (section == LpSection::ObjMin) {
            model_.objective_direction = ObjectiveDirection::Minimize;
        } else if (section == LpSection::ObjMax) {
            model_.objective_direction = ObjectiveDirection::Maximize;
        } else {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "expected \"Minimize\" or \"Maximize\" at the start of the .lp file.");
        }
        advance();

        // Optional "label:".
        if (is_label_start() && peek(1).type == TokenType::Colon) {
            advance();
            advance();
        }

        ParsedExpression expression = parse_linear_expression(/* stop_at_relop = */ false);
        for (const auto& term: expression.terms)
            model_.objective_coefficients[term.first] += term.second;
        for (const QuadraticTerm& term: expression.quadratic_terms) {
            model_.objective_quadratic_elements_variables_1.push_back(term.variable_1_id);
            model_.objective_quadratic_elements_variables_2.push_back(term.variable_2_id);
            model_.objective_quadratic_elements_coefficients.push_back(term.coefficient);
        }
        // expression.constant (the objective offset) is discarded: MathOptModel has no such field.
    }

    /** The "Subject To" section, which a model without constraints may omit. */
    void parse_constraints_section()
    {
        LpSection section = peek_section_keyword();
        if (section != LpSection::Constraints)
            return;
        consume_section_keyword(section);

        while (peek_section_keyword() == LpSection::None && !at_end())
            parse_constraint();
    }

    void parse_constraint()
    {
        std::string label;
        if (is_label_start() && peek(1).type == TokenType::Colon) {
            label = advance().text;
            advance();
        }

        bool leading_bound_value =
            (peek().type == TokenType::Number && is_relop(peek(1).type))
            || ((peek().type == TokenType::Plus || peek().type == TokenType::Minus)
                    && peek(1).type == TokenType::Number
                    && is_relop(peek(2).type));

        double lower = -inf;
        double upper = inf;
        ParsedExpression expression;

        if (leading_bound_value) {
            double value1 = parse_signed_number();
            RelOp relop1 = parse_relop();
            expression = parse_linear_expression(/* stop_at_relop = */ true);
            RelOp relop2 = parse_relop();
            double value2 = parse_signed_number();
            apply_range_relop(relop1, value1, relop2, value2, lower, upper);
        } else {
            expression = parse_linear_expression(/* stop_at_relop = */ true);
            RelOp relop = parse_relop();
            double value = parse_signed_number();
            apply_single_relop(relop, value, /* value_on_left = */ false, lower, upper);
        }

        model_.constraints_lower_bounds.push_back(lower);
        model_.constraints_upper_bounds.push_back(upper);
        model_.constraints_starts.push_back((int)model_.elements_variables.size());
        for (const auto& term: expression.terms) {
            model_.elements_variables.push_back(term.first);
            model_.elements_coefficients.push_back(term.second);
        }
        model_.quadratic_elements_constraints_starts.push_back(
                (int)model_.quadratic_elements_variables_1.size());
        for (const QuadraticTerm& term: expression.quadratic_terms) {
            model_.quadratic_elements_variables_1.push_back(term.variable_1_id);
            model_.quadratic_elements_variables_2.push_back(term.variable_2_id);
            model_.quadratic_elements_coefficients.push_back(term.coefficient);
        }
        constraint_labels_.push_back(label);
    }

    void parse_bounds_section()
    {
        while (peek_section_keyword() == LpSection::None && !at_end())
            parse_bound_line();
    }

    void parse_bound_line()
    {
        if (peek().type == TokenType::Identifier && peek_section_keyword() == LpSection::None) {
            if (at_identifier(1, "free")) {
                int variable_id = get_or_create_variable(advance().text);
                advance(); // "free"
                model_.variables_lower_bounds[variable_id] = -inf;
                model_.variables_upper_bounds[variable_id] = inf;
                return;
            }

            int variable_id = get_or_create_variable(advance().text);
            RelOp relop = parse_relop();
            double value = parse_signed_number();
            double lower = model_.variables_lower_bounds[variable_id];
            double upper = model_.variables_upper_bounds[variable_id];
            apply_single_relop(relop, value, /* value_on_left = */ false, lower, upper);
            model_.variables_lower_bounds[variable_id] = lower;
            model_.variables_upper_bounds[variable_id] = upper;
            return;
        }

        double value1 = parse_signed_number();
        RelOp relop1 = parse_relop();
        if (peek().type != TokenType::Identifier) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "expected a variable name in a Bounds entry.");
        }
        int variable_id = get_or_create_variable(advance().text);
        double lower = model_.variables_lower_bounds[variable_id];
        double upper = model_.variables_upper_bounds[variable_id];
        if (is_relop(peek().type)) {
            RelOp relop2 = parse_relop();
            double value2 = parse_signed_number();
            apply_range_relop(relop1, value1, relop2, value2, lower, upper);
        } else {
            apply_single_relop(relop1, value1, /* value_on_left = */ true, lower, upper);
        }
        model_.variables_lower_bounds[variable_id] = lower;
        model_.variables_upper_bounds[variable_id] = upper;
    }

    void parse_typed_variable_list(VariableType type)
    {
        while (peek_section_keyword() == LpSection::None && !at_end()) {
            if (peek().type != TokenType::Identifier) {
                throw std::invalid_argument(
                        FUNC_SIGNATURE + ": "
                        "expected a variable name in a General/Binary/Semi-continuous section.");
            }
            int variable_id = get_or_create_variable(advance().text);
            // A variable listed in both the General and the Semi-continuous
            // sections is semi-integer.
            VariableType current = model_.variables_types[variable_id];
            if ((type == VariableType::SemiContinuous && current == VariableType::Integer)
                    || (type == VariableType::Integer && current == VariableType::SemiContinuous)) {
                model_.variables_types[variable_id] = VariableType::SemiInteger;
                continue;
            }
            model_.variables_types[variable_id] = type;
            if (type == VariableType::Binary
                    && model_.variables_upper_bounds[variable_id] == inf) {
                model_.variables_upper_bounds[variable_id] = 1.0;
            }
        }
    }

    MathOptModel finalize()
    {
        // The constraints' quadratic starts are only kept if some constraint
        // has quadratic terms.
        if (model_.quadratic_elements_variables_1.empty())
            model_.quadratic_elements_constraints_starts.clear();

        bool any_label = false;
        for (const std::string& label: constraint_labels_)
            if (!label.empty())
                any_label = true;
        if (any_label) {
            model_.constraints_names.resize(constraint_labels_.size());
            for (size_t constraint_id = 0; constraint_id < constraint_labels_.size(); ++constraint_id) {
                model_.constraints_names[constraint_id] = constraint_labels_[constraint_id].empty()?
                    ("c" + std::to_string(constraint_id)):
                    constraint_labels_[constraint_id];
            }
        }
        return std::move(model_);
    }

};

/*
 * Writer helpers.
 */

bool is_valid_number_literal(const std::string& s)
{
    if (s.empty())
        return false;
    try {
        std::vector<Token> number_tokens = tokenize(s);
        return number_tokens.size() == 2 && number_tokens[0].type == TokenType::Number;
    } catch (const std::invalid_argument&) {
        return false;
    }
}

/**
 * 'allow_numeric' accepts a name that reads back as a plain number (e.g.
 * "2"), which is safe for a constraint name (read_lp() only recognizes it
 * as a label immediately before ':', matching HiGHS's own .lp reader) but
 * not for a variable name (bare in an expression, indistinguishable from
 * a coefficient).
 */
void check_lp_name(const std::string& name, bool allow_numeric = false)
{
    bool ok = !name.empty() && is_identifier_start(name[0]);
    for (size_t i = 1; ok && i < name.size(); ++i)
        ok = is_identifier_char(name[i]);
    // Only run the tokenizer for the names that aren't identifiers.
    if (!ok && allow_numeric && is_valid_number_literal(name))
        return;
    if (!ok) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "name \"" + name + "\" is not representable in .lp format "
                "(must start with a letter, '_' or '.' and contain only "
                "letters, digits, '_' and '.'" + (allow_numeric? ", or be a plain number": "") + ").");
    }
}

using text_output::TextOutput;
using text_output::Names;

/** A bound: a number, or "infinity" / "-infinity". */
void write_value(TextOutput& out, double value)
{
    if (value == inf)
        out << "infinity";
    else if (value == -inf)
        out << "-infinity";
    else
        out.number(value);
}

/** The sign and coefficient of a term: " +", " -", " +2" or " -2". */
void write_coefficient(
        TextOutput& out,
        double coefficient)
{
    out << ' ';
    if (coefficient == 1.0) {
        out << '+';
    } else if (coefficient == -1.0) {
        out << '-';
    } else if (coefficient >= 0.0) {
        out << '+';
        out.number(coefficient);
    } else {
        out << '-';
        out.number(-coefficient);
    }
}

void write_term(
        TextOutput& out,
        double coefficient,
        const Names& names,
        int variable_id)
{
    write_coefficient(out, coefficient);
    out << ' ';
    names.write(out, variable_id);
}

/**
 * Quadratic terms, "[ c x ^ 2 + c x * y ... ]", with their coefficients
 * multiplied by 'factor' (2 for the objective, written as "[ ... ] / 2").
 */
void write_quadratic_terms(
        TextOutput& out,
        const Names& names,
        const std::vector<int>& variables_1,
        const std::vector<int>& variables_2,
        const std::vector<double>& coefficients,
        int start,
        int end,
        double factor)
{
    out << " + [";
    for (int element_id = start; element_id < end; ++element_id) {
        double coefficient = coefficients[element_id];
        if (coefficient == 0.0)
            continue;
        write_coefficient(out, coefficient * factor);
        out << ' ';
        names.write(out, variables_1[element_id]);
        if (variables_1[element_id] == variables_2[element_id]) {
            out << " ^ 2";
        } else {
            out << " * ";
            names.write(out, variables_2[element_id]);
        }
    }
    out << " ]";
}

void write_expression(
        TextOutput& out,
        const MathOptModel& model,
        const Names& variable_names,
        int start,
        int end)
{
    for (int element_id = start; element_id < end; ++element_id) {
        double coefficient = model.elements_coefficients[element_id];
        if (coefficient == 0.0)
            continue;
        write_term(out, coefficient, variable_names, model.elements_variables[element_id]);
    }
}

}

MathOptModel mathoptsolverscmake::read_lp(const std::string& file_path)
{
    std::ifstream file(file_path);
    if (!file.good()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "unable to open file \"" + file_path + "\".");
    }
    return read_lp(file);
}

MathOptModel mathoptsolverscmake::read_lp(std::istream& in)
{
    std::string content(
            (std::istreambuf_iterator<char>(in)),
            std::istreambuf_iterator<char>());
    LpParser parser(content);
    return parser.parse();
}

void mathoptsolverscmake::write_lp(
        const MathOptModel& model,
        const std::string& file_path)
{
    std::ofstream file(file_path);
    if (!file.good()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "unable to open file \"" + file_path + "\" for writing.");
    }
    write_lp(model, file);
}

void mathoptsolverscmake::write_lp(
        const MathOptModel& model,
        std::ostream& out)
{
    if (model.has_nonlinear() || model.has_black_box()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "model has nonlinear terms; .lp format only supports linear "
                "and quadratic models.");
    }

    int number_of_variables = model.number_of_variables();
    int number_of_constraints = model.number_of_constraints();

    // Check the model's names (the default ones, "x<id>" and "c<id>", are
    // always valid) before writing anything.
    if (!model.variables_names.empty()) {
        for (int variable_id = 0; variable_id < number_of_variables; ++variable_id)
            check_lp_name(model.variables_names[variable_id]);
    }
    if (!model.constraints_names.empty()) {
        for (int constraint_id = 0; constraint_id < number_of_constraints; ++constraint_id)
            check_lp_name(model.constraints_names[constraint_id], /* allow_numeric = */ true);
    }
    Names variable_names(model.variables_names, 'x');
    Names constraint_names(model.constraints_names, 'c');

    TextOutput output(out);

    output << (model.objective_direction == ObjectiveDirection::Maximize? "Maximize": "Minimize") << "\n";
    output << " obj:";
    for (int variable_id = 0; variable_id < number_of_variables; ++variable_id) {
        double coefficient = model.objective_coefficients[variable_id];
        if (coefficient != 0.0)
            write_term(output, coefficient, variable_names, variable_id);
    }
    if (!model.objective_quadratic_elements_coefficients.empty()) {
        write_quadratic_terms(
                output,
                variable_names,
                model.objective_quadratic_elements_variables_1,
                model.objective_quadratic_elements_variables_2,
                model.objective_quadratic_elements_coefficients,
                0,
                (int)model.objective_quadratic_elements_coefficients.size(),
                2.0);
        output << " / 2";
    }
    output << "\n";

    output << "Subject To\n";
    for (int constraint_id = 0; constraint_id < number_of_constraints; ++constraint_id) {
        ConstraintSense sense = model.constraint_sense(constraint_id);
        // A range needs its lower bound written between the label and the
        // expression, as "name: lo <= expr <= up" (read_lp() expects the
        // label first); every other sense has a single bound written after.
        output << ' ';
        constraint_names.write(output, constraint_id);
        output << ':';
        if (sense == ConstraintSense::Range) {
            output << ' ';
            write_value(output, model.constraints_lower_bounds[constraint_id]);
            output << " <=";
        }
        write_expression(
                output,
                model,
                variable_names,
                model.constraints_starts[constraint_id],
                model.constraint_end(constraint_id));
        if (!model.quadratic_elements_constraints_starts.empty()) {
            int start = model.quadratic_elements_constraints_starts[constraint_id];
            int end = model.quadratic_constraint_end(constraint_id);
            if (start < end) {
                write_quadratic_terms(
                        output,
                        variable_names,
                        model.quadratic_elements_variables_1,
                        model.quadratic_elements_variables_2,
                        model.quadratic_elements_coefficients,
                        start,
                        end,
                        1.0);
            }
        }

        switch (sense) {
        case ConstraintSense::LessThanOrEqualTo: {
            output << " <=";
            write_value(output, model.constraints_upper_bounds[constraint_id]);
            break;
        } case ConstraintSense::GreaterThanOrEqualTo: {
            output << " >=";
            write_value(output, model.constraints_lower_bounds[constraint_id]);
            break;
        } case ConstraintSense::Equality: {
            output << " =";
            write_value(output, model.constraints_lower_bounds[constraint_id]);
            break;
        } case ConstraintSense::Range: {
            output << " <=";
            write_value(output, model.constraints_upper_bounds[constraint_id]);
            break;
        } case ConstraintSense::Free: {
            // .lp has no direct syntax for an unconstrained row.
            output << " >= -infinity";
            break;
        }
        }
        output << '\n';
    }

    output << "Bounds\n";
    for (int variable_id = 0; variable_id < number_of_variables; ++variable_id) {
        double lower = model.variables_lower_bounds[variable_id];
        double upper = model.variables_upper_bounds[variable_id];
        bool is_binary = (model.variables_types[variable_id] == VariableType::Binary);
        double implied_upper = is_binary? 1.0: inf;

        if (lower == 0.0 && upper == implied_upper)
            continue;

        if (lower == -inf && upper == inf) {
            output << ' ';
            variable_names.write(output, variable_id);
            output << " free\n";
        } else if (lower == upper) {
            output << ' ';
            variable_names.write(output, variable_id);
            output << " =";
            write_value(output, lower);
            output << '\n';
        } else {
            if (lower != 0.0) {
                write_value(output, lower);
                output << " <=";
            }
            output << ' ';
            variable_names.write(output, variable_id);
            if (upper != inf) {
                output << " <=";
                write_value(output, upper);
            }
            output << '\n';
        }
    }

    // A semi-integer variable is listed in both the General and the
    // Semi-continuous sections.
    bool has_general = false;
    bool has_binary = false;
    bool has_semi = false;
    for (int variable_id = 0; variable_id < number_of_variables; ++variable_id) {
        VariableType type = model.variables_types[variable_id];
        if (type == VariableType::Integer || type == VariableType::SemiInteger)
            has_general = true;
        else if (type == VariableType::Binary)
            has_binary = true;
        if (is_semi(type))
            has_semi = true;
    }
    if (has_general) {
        output << "General\n";
        for (int variable_id = 0; variable_id < number_of_variables; ++variable_id) {
            VariableType type = model.variables_types[variable_id];
            if (type == VariableType::Integer || type == VariableType::SemiInteger) {
                output << ' ';
                variable_names.write(output, variable_id);
                output << '\n';
            }
        }
    }
    if (has_binary) {
        output << "Binary\n";
        for (int variable_id = 0; variable_id < number_of_variables; ++variable_id) {
            if (model.variables_types[variable_id] == VariableType::Binary) {
                output << ' ';
                variable_names.write(output, variable_id);
                output << '\n';
            }
        }
    }

    if (has_semi) {
        output << "Semi-continuous\n";
        for (int variable_id = 0; variable_id < number_of_variables; ++variable_id) {
            if (is_semi(model.variables_types[variable_id])) {
                output << ' ';
                variable_names.write(output, variable_id);
                output << '\n';
            }
        }
    }

    output << "End\n";
    output.flush();
}

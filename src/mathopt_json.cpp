#include "mathoptsolverscmake/mathopt_json.hpp"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

using namespace mathoptsolverscmake;

namespace
{

/*
 * A minimal JSON DOM and recursive-descent parser, just expressive enough
 * to read back what write_json() produces (objects, arrays, strings with
 * escapes, numbers, true/false/null, and the bare Infinity/-Infinity/NaN
 * tokens write_json() uses for non-finite doubles).
 */

struct JsonValue
{
    enum class Type { Null, Bool, Number, String, Array, Object };

    Type type = Type::Null;
    bool bool_value = false;
    double number_value = 0.0;
    std::string string_value;
    std::vector<JsonValue> array_value;
    std::vector<std::pair<std::string, JsonValue>> object_value;

    const JsonValue* find(const std::string& key) const
    {
        for (const auto& entry: object_value)
            if (entry.first == key)
                return &entry.second;
        return nullptr;
    }
};

class JsonParser
{

public:

    explicit JsonParser(const std::string& text):
        text_(text),
        pos_(0)
    { }

    JsonValue parse()
    {
        JsonValue value = parse_value();
        skip_whitespace();
        if (pos_ != text_.size()) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "trailing content after the top-level JSON value.");
        }
        return value;
    }

private:

    const std::string& text_;
    size_t pos_;

    char peek() const
    {
        if (pos_ >= text_.size()) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "unexpected end of JSON input.");
        }
        return text_[pos_];
    }

    void skip_whitespace()
    {
        while (pos_ < text_.size() && std::isspace((unsigned char)text_[pos_]))
            ++pos_;
    }

    void expect(char c)
    {
        if (peek() != c) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "expected '" + std::string(1, c) + "' at position " + std::to_string(pos_) + ".");
        }
        ++pos_;
    }

    void expect_literal(const std::string& literal)
    {
        if (text_.compare(pos_, literal.size(), literal) != 0) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "expected \"" + literal + "\" at position " + std::to_string(pos_) + ".");
        }
        pos_ += literal.size();
    }

    JsonValue parse_value()
    {
        skip_whitespace();
        char c = peek();
        if (c == '{')
            return parse_object();
        if (c == '[')
            return parse_array();
        if (c == '"')
            return parse_string_value();
        if (c == 't') {
            expect_literal("true");
            JsonValue value; value.type = JsonValue::Type::Bool; value.bool_value = true;
            return value;
        }
        if (c == 'f') {
            expect_literal("false");
            JsonValue value; value.type = JsonValue::Type::Bool; value.bool_value = false;
            return value;
        }
        if (c == 'n') {
            expect_literal("null");
            JsonValue value; value.type = JsonValue::Type::Null;
            return value;
        }
        if (c == 'N') {
            expect_literal("NaN");
            JsonValue value; value.type = JsonValue::Type::Number;
            value.number_value = std::numeric_limits<double>::quiet_NaN();
            return value;
        }
        if (c == 'I') {
            expect_literal("Infinity");
            JsonValue value; value.type = JsonValue::Type::Number;
            value.number_value = std::numeric_limits<double>::infinity();
            return value;
        }
        return parse_number();
    }

    JsonValue parse_number()
    {
        size_t start = pos_;
        if (peek() == '-') {
            ++pos_;
            if (pos_ < text_.size() && text_[pos_] == 'I') {
                expect_literal("Infinity");
                JsonValue value; value.type = JsonValue::Type::Number;
                value.number_value = -std::numeric_limits<double>::infinity();
                return value;
            }
        }
        while (pos_ < text_.size()
                && (std::isdigit((unsigned char)text_[pos_])
                    || text_[pos_] == '.'
                    || text_[pos_] == 'e'
                    || text_[pos_] == 'E'
                    || text_[pos_] == '+'
                    || text_[pos_] == '-')) {
            ++pos_;
        }
        std::string token = text_.substr(start, pos_ - start);
        if (token.empty() || token == "-") {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "invalid number at position " + std::to_string(start) + ".");
        }
        JsonValue value;
        value.type = JsonValue::Type::Number;
        value.number_value = std::strtod(token.c_str(), nullptr);
        return value;
    }

    std::string parse_raw_string()
    {
        expect('"');
        std::string result;
        while (true) {
            if (pos_ >= text_.size()) {
                throw std::invalid_argument(
                        FUNC_SIGNATURE + ": "
                        "unterminated JSON string.");
            }
            char c = text_[pos_++];
            if (c == '"')
                break;
            if (c != '\\') {
                result += c;
                continue;
            }
            if (pos_ >= text_.size()) {
                throw std::invalid_argument(
                        FUNC_SIGNATURE + ": "
                        "unterminated escape sequence in JSON string.");
            }
            char escaped = text_[pos_++];
            switch (escaped) {
            case '"': result += '"'; break;
            case '\\': result += '\\'; break;
            case '/': result += '/'; break;
            case 'b': result += '\b'; break;
            case 'f': result += '\f'; break;
            case 'n': result += '\n'; break;
            case 'r': result += '\r'; break;
            case 't': result += '\t'; break;
            case 'u': {
                if (pos_ + 4 > text_.size()) {
                    throw std::invalid_argument(
                            FUNC_SIGNATURE + ": "
                            "truncated \\u escape in JSON string.");
                }
                unsigned int code = (unsigned int)std::strtoul(text_.substr(pos_, 4).c_str(), nullptr, 16);
                pos_ += 4;
                // Encode as UTF-8; surrogate pairs (code points beyond the
                // basic multilingual plane) are not supported, as variable
                // and constraint names are expected to be plain ASCII.
                if (code < 0x80) {
                    result += (char)code;
                } else if (code < 0x800) {
                    result += (char)(0xC0 | (code >> 6));
                    result += (char)(0x80 | (code & 0x3F));
                } else {
                    result += (char)(0xE0 | (code >> 12));
                    result += (char)(0x80 | ((code >> 6) & 0x3F));
                    result += (char)(0x80 | (code & 0x3F));
                }
                break;
            }
            default:
                throw std::invalid_argument(
                        FUNC_SIGNATURE + ": "
                        "invalid escape sequence \"\\" + std::string(1, escaped) + "\" in JSON string.");
            }
        }
        return result;
    }

    JsonValue parse_string_value()
    {
        JsonValue value;
        value.type = JsonValue::Type::String;
        value.string_value = parse_raw_string();
        return value;
    }

    JsonValue parse_array()
    {
        JsonValue value;
        value.type = JsonValue::Type::Array;
        expect('[');
        skip_whitespace();
        if (peek() == ']') {
            ++pos_;
            return value;
        }
        while (true) {
            value.array_value.push_back(parse_value());
            skip_whitespace();
            char c = peek();
            ++pos_;
            if (c == ',')
                continue;
            if (c == ']')
                break;
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "expected ',' or ']' in JSON array.");
        }
        return value;
    }

    JsonValue parse_object()
    {
        JsonValue value;
        value.type = JsonValue::Type::Object;
        expect('{');
        skip_whitespace();
        if (peek() == '}') {
            ++pos_;
            return value;
        }
        while (true) {
            skip_whitespace();
            std::string key = parse_raw_string();
            skip_whitespace();
            expect(':');
            JsonValue field_value = parse_value();
            value.object_value.push_back(std::make_pair(key, std::move(field_value)));
            skip_whitespace();
            char c = peek();
            ++pos_;
            if (c == ',')
                continue;
            if (c == '}')
                break;
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "expected ',' or '}' in JSON object.");
        }
        return value;
    }

};

/*
 * Extraction helpers: a missing field defaults to an empty vector (or the
 * given default for a scalar); a present-but-wrong-typed field throws.
 */

void require_array(const JsonValue& value, const std::string& key)
{
    if (value.type != JsonValue::Type::Array) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "field \"" + key + "\" must be an array.");
    }
}

std::vector<double> get_double_array(const JsonValue& object, const std::string& key)
{
    const JsonValue* value = object.find(key);
    if (!value)
        return {};
    require_array(*value, key);
    std::vector<double> result(value->array_value.size());
    for (size_t element_id = 0; element_id < result.size(); ++element_id) {
        const JsonValue& element = value->array_value[element_id];
        if (element.type != JsonValue::Type::Number) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "element " + std::to_string(element_id) + " of \"" + key + "\" must be a number.");
        }
        result[element_id] = element.number_value;
    }
    return result;
}

std::vector<int> get_int_array(const JsonValue& object, const std::string& key)
{
    std::vector<double> doubles = get_double_array(object, key);
    std::vector<int> result(doubles.size());
    for (size_t element_id = 0; element_id < doubles.size(); ++element_id)
        result[element_id] = (int)doubles[element_id];
    return result;
}

std::vector<std::string> get_string_array(const JsonValue& object, const std::string& key)
{
    const JsonValue* value = object.find(key);
    if (!value)
        return {};
    require_array(*value, key);
    std::vector<std::string> result(value->array_value.size());
    for (size_t element_id = 0; element_id < result.size(); ++element_id) {
        const JsonValue& element = value->array_value[element_id];
        if (element.type != JsonValue::Type::String) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "element " + std::to_string(element_id) + " of \"" + key + "\" must be a string.");
        }
        result[element_id] = element.string_value;
    }
    return result;
}

std::vector<char> get_char_array(const JsonValue& object, const std::string& key)
{
    std::vector<std::string> strings = get_string_array(object, key);
    std::vector<char> result(strings.size());
    for (size_t element_id = 0; element_id < strings.size(); ++element_id) {
        if (strings[element_id].size() != 1) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "element " + std::to_string(element_id) + " of \"" + key + "\" must be a single-character string.");
        }
        result[element_id] = strings[element_id][0];
    }
    return result;
}

double get_double(const JsonValue& object, const std::string& key, double default_value)
{
    const JsonValue* value = object.find(key);
    if (!value)
        return default_value;
    if (value->type != JsonValue::Type::Number) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "field \"" + key + "\" must be a number.");
    }
    return value->number_value;
}

/*
 * Writer helpers.
 */

std::string escape_json_string(const std::string& s)
{
    std::string result;
    result.reserve(s.size() + 2);
    for (char c: s) {
        switch (c) {
        case '"': result += "\\\""; break;
        case '\\': result += "\\\\"; break;
        case '\b': result += "\\b"; break;
        case '\f': result += "\\f"; break;
        case '\n': result += "\\n"; break;
        case '\r': result += "\\r"; break;
        case '\t': result += "\\t"; break;
        default:
            if ((unsigned char)c < 0x20) {
                char buffer[8];
                std::snprintf(buffer, sizeof(buffer), "\\u%04x", (unsigned char)c);
                result += buffer;
            } else {
                result += c;
            }
        }
    }
    return result;
}

void write_double(std::ostream& out, double value)
{
    if (std::isnan(value))
        out << "NaN";
    else if (std::isinf(value))
        out << (value > 0? "Infinity": "-Infinity");
    else
        out << value;
}

void write_double_array(std::ostream& out, const std::vector<double>& values)
{
    out << "[";
    for (size_t element_id = 0; element_id < values.size(); ++element_id) {
        if (element_id > 0)
            out << ", ";
        write_double(out, values[element_id]);
    }
    out << "]";
}

void write_int_array(std::ostream& out, const std::vector<int>& values)
{
    out << "[";
    for (size_t element_id = 0; element_id < values.size(); ++element_id) {
        if (element_id > 0)
            out << ", ";
        out << values[element_id];
    }
    out << "]";
}

void write_string_array(std::ostream& out, const std::vector<std::string>& values)
{
    out << "[";
    for (size_t element_id = 0; element_id < values.size(); ++element_id) {
        if (element_id > 0)
            out << ", ";
        out << "\"" << escape_json_string(values[element_id]) << "\"";
    }
    out << "]";
}

void write_char_array(std::ostream& out, const std::vector<char>& values)
{
    out << "[";
    for (size_t element_id = 0; element_id < values.size(); ++element_id) {
        if (element_id > 0)
            out << ", ";
        out << "\"" << escape_json_string(std::string(1, values[element_id])) << "\"";
    }
    out << "]";
}

}

MathOptModel mathoptsolverscmake::read_json(const std::string& file_path)
{
    std::ifstream file(file_path);
    if (!file.good()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "unable to open file \"" + file_path + "\".");
    }
    return read_json(file);
}

MathOptModel mathoptsolverscmake::read_json(std::istream& in)
{
    std::string content(
            (std::istreambuf_iterator<char>(in)),
            std::istreambuf_iterator<char>());

    JsonParser parser(content);
    JsonValue root = parser.parse();
    if (root.type != JsonValue::Type::Object) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "the top-level JSON value must be an object.");
    }

    MathOptModel model;

    const JsonValue* direction_value = root.find("objective_direction");
    if (!direction_value || direction_value->type != JsonValue::Type::String) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "missing or invalid \"objective_direction\" field.");
    }
    {
        std::istringstream direction_stream(direction_value->string_value);
        direction_stream >> model.objective_direction;
        if (!direction_stream) {
            throw std::invalid_argument(
                    FUNC_SIGNATURE + ": "
                    "invalid \"objective_direction\" value \"" + direction_value->string_value + "\".");
        }
    }

    model.variables_lower_bounds = get_double_array(root, "variables_lower_bounds");
    model.variables_upper_bounds = get_double_array(root, "variables_upper_bounds");
    {
        std::vector<std::string> type_names = get_string_array(root, "variables_types");
        if (type_names.empty()) {
            model.variables_types.assign(model.variables_lower_bounds.size(), VariableType::Continuous);
        } else {
            model.variables_types.resize(type_names.size());
            for (size_t variable_id = 0; variable_id < type_names.size(); ++variable_id) {
                std::istringstream type_stream(type_names[variable_id]);
                type_stream >> model.variables_types[variable_id];
                if (!type_stream) {
                    throw std::invalid_argument(
                            FUNC_SIGNATURE + ": "
                            "invalid variable type \"" + type_names[variable_id] + "\".");
                }
            }
        }
    }
    model.variables_names = get_string_array(root, "variables_names");
    model.variables_initial_values = get_double_array(root, "variables_initial_values");

    model.constraints_lower_bounds = get_double_array(root, "constraints_lower_bounds");
    model.constraints_upper_bounds = get_double_array(root, "constraints_upper_bounds");
    model.constraints_names = get_string_array(root, "constraints_names");

    model.objective_coefficients = get_double_array(root, "objective_coefficients");
    model.constraints_starts = get_int_array(root, "constraints_starts");
    model.elements_variables = get_int_array(root, "elements_variables");
    model.elements_coefficients = get_double_array(root, "elements_coefficients");

    model.objective_quadratic_elements_variables_1 = get_int_array(root, "objective_quadratic_elements_variables_1");
    model.objective_quadratic_elements_variables_2 = get_int_array(root, "objective_quadratic_elements_variables_2");
    model.objective_quadratic_elements_coefficients = get_double_array(root, "objective_quadratic_elements_coefficients");
    model.quadratic_elements_constraints_starts = get_int_array(root, "quadratic_elements_constraints_starts");
    model.quadratic_elements_variables_1 = get_int_array(root, "quadratic_elements_variables_1");
    model.quadratic_elements_variables_2 = get_int_array(root, "quadratic_elements_variables_2");
    model.quadratic_elements_coefficients = get_double_array(root, "quadratic_elements_coefficients");

    model.objective_nonlinear_elements_operators = get_char_array(root, "objective_nonlinear_elements_operators");
    model.objective_nonlinear_elements_values = get_double_array(root, "objective_nonlinear_elements_values");
    model.objective_nonlinear_elements_variables = get_int_array(root, "objective_nonlinear_elements_variables");
    model.objective_nonlinear_elements_number_of_children = get_int_array(root, "objective_nonlinear_elements_number_of_children");

    model.nonlinear_elements_constraints_starts = get_int_array(root, "nonlinear_elements_constraints_starts");
    model.nonlinear_elements_operators = get_char_array(root, "nonlinear_elements_operators");
    model.nonlinear_elements_values = get_double_array(root, "nonlinear_elements_values");
    model.nonlinear_elements_variables = get_int_array(root, "nonlinear_elements_variables");
    model.nonlinear_elements_number_of_children = get_int_array(root, "nonlinear_elements_number_of_children");

    model.feasibility_tolerance = get_double(root, "feasibility_tolerance", 0.0);
    model.integrality_tolerance = get_double(root, "integrality_tolerance", 0.0);

    return model;
}

void mathoptsolverscmake::write_json(
        const MathOptModel& model,
        const std::string& file_path)
{
    std::ofstream file(file_path);
    if (!file.good()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "unable to open file \"" + file_path + "\" for writing.");
    }
    write_json(model, file);
}

void mathoptsolverscmake::write_json(
        const MathOptModel& model,
        std::ostream& out)
{
    if (model.has_black_box()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "model has a black-box objective or constraint function, "
                "which has no JSON representation.");
    }

    out << std::setprecision(17);

    bool first_field = true;
    auto field = [&](const std::string& name) -> std::ostream&
    {
        if (!first_field)
            out << ",\n";
        first_field = false;
        out << "  \"" << name << "\": ";
        return out;
    };

    out << "{\n";

    {
        std::ostringstream direction_stream;
        direction_stream << model.objective_direction;
        field("objective_direction") << "\"" << direction_stream.str() << "\"";
    }

    write_double_array(field("variables_lower_bounds"), model.variables_lower_bounds);
    write_double_array(field("variables_upper_bounds"), model.variables_upper_bounds);
    {
        std::vector<std::string> type_names(model.variables_types.size());
        for (size_t variable_id = 0; variable_id < type_names.size(); ++variable_id) {
            std::ostringstream type_stream;
            type_stream << model.variables_types[variable_id];
            type_names[variable_id] = type_stream.str();
        }
        write_string_array(field("variables_types"), type_names);
    }
    if (!model.variables_names.empty())
        write_string_array(field("variables_names"), model.variables_names);
    if (!model.variables_initial_values.empty())
        write_double_array(field("variables_initial_values"), model.variables_initial_values);

    write_double_array(field("constraints_lower_bounds"), model.constraints_lower_bounds);
    write_double_array(field("constraints_upper_bounds"), model.constraints_upper_bounds);
    if (!model.constraints_names.empty())
        write_string_array(field("constraints_names"), model.constraints_names);

    write_double_array(field("objective_coefficients"), model.objective_coefficients);
    write_int_array(field("constraints_starts"), model.constraints_starts);
    write_int_array(field("elements_variables"), model.elements_variables);
    write_double_array(field("elements_coefficients"), model.elements_coefficients);

    if (!model.objective_quadratic_elements_variables_1.empty()) {
        write_int_array(field("objective_quadratic_elements_variables_1"), model.objective_quadratic_elements_variables_1);
        write_int_array(field("objective_quadratic_elements_variables_2"), model.objective_quadratic_elements_variables_2);
        write_double_array(field("objective_quadratic_elements_coefficients"), model.objective_quadratic_elements_coefficients);
    }
    if (!model.quadratic_elements_variables_1.empty()) {
        write_int_array(field("quadratic_elements_constraints_starts"), model.quadratic_elements_constraints_starts);
        write_int_array(field("quadratic_elements_variables_1"), model.quadratic_elements_variables_1);
        write_int_array(field("quadratic_elements_variables_2"), model.quadratic_elements_variables_2);
        write_double_array(field("quadratic_elements_coefficients"), model.quadratic_elements_coefficients);
    }

    if (!model.objective_nonlinear_elements_operators.empty()) {
        write_char_array(field("objective_nonlinear_elements_operators"), model.objective_nonlinear_elements_operators);
        write_double_array(field("objective_nonlinear_elements_values"), model.objective_nonlinear_elements_values);
        write_int_array(field("objective_nonlinear_elements_variables"), model.objective_nonlinear_elements_variables);
        write_int_array(field("objective_nonlinear_elements_number_of_children"), model.objective_nonlinear_elements_number_of_children);
    }
    if (!model.nonlinear_elements_operators.empty()) {
        write_int_array(field("nonlinear_elements_constraints_starts"), model.nonlinear_elements_constraints_starts);
        write_char_array(field("nonlinear_elements_operators"), model.nonlinear_elements_operators);
        write_double_array(field("nonlinear_elements_values"), model.nonlinear_elements_values);
        write_int_array(field("nonlinear_elements_variables"), model.nonlinear_elements_variables);
        write_int_array(field("nonlinear_elements_number_of_children"), model.nonlinear_elements_number_of_children);
    }

    write_double(field("feasibility_tolerance"), model.feasibility_tolerance);
    write_double(field("integrality_tolerance"), model.integrality_tolerance);

    out << "\n}\n";
}

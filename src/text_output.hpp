#pragma once

/**
 * Buffered text output shared by the file writers (LP, MPS, .nl, JSON).
 *
 * The text is built in a string and written to the stream in large chunks,
 * and numbers are formatted without going through the stream. The output is
 * the same as with an std::ostream using the classic locale and
 * std::setprecision(17), which is what the writers used before.
 *
 * Internal header: not installed.
 */

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <ostream>
#include <string>
#include <vector>

namespace mathoptsolverscmake
{
namespace text_output
{

class TextOutput
{

public:

    explicit TextOutput(std::ostream& out): out_(out)
    {
        buffer_.reserve(CAPACITY + 1024);
    }

    /** Write the buffered text to the stream; to call once done writing. */
    void flush()
    {
        out_.write(buffer_.data(), buffer_.size());
        buffer_.clear();
    }

    TextOutput& operator<<(char c)
    {
        buffer_.push_back(c);
        return check_flush();
    }

    TextOutput& operator<<(const char* text)
    {
        buffer_.append(text);
        return check_flush();
    }

    TextOutput& operator<<(const std::string& text)
    {
        buffer_.append(text);
        return check_flush();
    }

    /** An integer, in decimal. */
    TextOutput& integer(long long value)
    {
        char digits[24];
        int size = 0;
        unsigned long long magnitude = (value < 0)?
            0ULL - (unsigned long long)value: (unsigned long long)value;
        do {
            digits[size++] = (char)('0' + magnitude % 10);
            magnitude /= 10;
        } while (magnitude != 0);
        if (value < 0)
            buffer_.push_back('-');
        while (size > 0)
            buffer_.push_back(digits[--size]);
        return check_flush();
    }

    /** A number, as "%.17g" writes it ("inf", "-inf" and "nan" included). */
    TextOutput& number(double value)
    {
        // Integer values (most coefficients and bounds) are written directly.
        if (value == std::floor(value) && std::fabs(value) < 9007199254740992.0) {
            if (value == 0.0 && std::signbit(value))
                return *this << "-0";
            return integer((long long)value);
        }
        char text[32];
        int size = std::snprintf(text, sizeof(text), "%.17g", value);
        // snprintf follows the C locale's decimal separator; files need '.'.
        for (int i = 0; i < size; ++i)
            if (text[i] == ',')
                text[i] = '.';
        buffer_.append(text, size);
        return check_flush();
    }

private:

    static const std::size_t CAPACITY = 1 << 20;

    TextOutput& check_flush()
    {
        if (buffer_.size() >= CAPACITY)
            flush();
        return *this;
    }

    std::ostream& out_;

    std::string buffer_;

};

/**
 * Names of the variables or constraints: the model's ones, or, if it has none,
 * the default ones, prefix followed by the index ("x0", "c0"...), as
 * MathOptModel::variable_name() and constraint_name() return them.
 */
class Names
{

public:

    Names(const std::vector<std::string>& names, char prefix):
        names_(names), prefix_(prefix) { }

    void write(TextOutput& out, int id) const
    {
        if (names_.empty()) {
            out << prefix_;
            out.integer(id);
        } else {
            out << names_[id];
        }
    }

    /** True if the model has its own names (rather than the default ones). */
    bool has_names() const { return !names_.empty(); }

private:

    const std::vector<std::string>& names_;

    char prefix_;

};

}
}

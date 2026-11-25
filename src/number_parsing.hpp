#pragma once

/**
 * Locale-independent number parsing, for the file readers.
 *
 * std::strtod follows the decimal separator of the C locale, which isn't
 * always "C": the JVM sets it from the environment, and so do programs that
 * call setlocale(LC_ALL, ""). With a French locale, for example, strtod reads
 * "1.5" as 1. Files always use '.', so the readers parse numbers with the
 * "C" locale explicitly.
 *
 * Internal header: not installed.
 */

#include <clocale>
#include <cstdlib>

#if defined(_WIN32)
#  include <locale.h>
#elif defined(__APPLE__)
#  include <xlocale.h>
#else
#  include <locale.h>
#endif

namespace mathoptsolverscmake
{
namespace number_parsing
{

/** std::strtod with the "C" locale, whatever the current locale. */
inline double c_strtod(const char* text, char** end)
{
#if defined(_WIN32)
    static const _locale_t c_locale = _create_locale(LC_NUMERIC, "C");
    return _strtod_l(text, end, c_locale);
#else
    static const locale_t c_locale = newlocale(LC_NUMERIC_MASK, "C", (locale_t)0);
    return strtod_l(text, end, c_locale);
#endif
}

}
}

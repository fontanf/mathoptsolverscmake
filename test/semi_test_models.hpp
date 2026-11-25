#pragma once

/**
 * HiGHS's check/instances/semi-continuous and semi-integer models (MIT
 * license), whose optimal values are 8.22333 and 8.13333 (see HiGHS's
 * check/TestSemiVariables.cpp): x3 is semi-continuous / semi-integer, either
 * 0 or in [1.1, 10].
 */

#include <string>

namespace mathoptsolverscmaketest
{

inline std::string semi_mps(bool semi_integer)
{
    return std::string(
            "NAME lp_solve_sc\n"
            "OBJSENSE\n"
            "  MAX\n"
            "ROWS\n"
            " N  r_0\n"
            " L  r_1\n"
            " G  r_2\n"
            " G  r_3\n"
            " G  r_4\n"
            "COLUMNS\n"
            "    x1        r_0                  1   r_1                  1\n"
            "    x1        r_2                  2   r_3                 -1\n"
            "    x2        r_0                  2   r_1                  1\n"
            "    x2        r_2                 -1   r_3                  3\n"
            "    MARK0000  'MARKER'                 'INTORG'\n"
            "    x3        r_0               -0.1   r_4                  1\n"
            "    MARK0001  'MARKER'                 'INTEND'\n"
            "    x4        r_0                 -3   r_4                  1\n"
            "RHS\n"
            "    RHS       r_1                  5   r_4                0.5\n"
            "BOUNDS\n")
        + (semi_integer? " SI BND       x3                  10\n": " SC BND       x3                  10\n")
        + " LO BND       x3                 1.1\n"
          "ENDATA\n";
}

inline std::string semi_lp(bool semi_integer)
{
    return std::string(
            "Maximize\n"
            "x1 + 2 x2 - 0.1 x3 -3 x4\n"
            "Subject To\n"
            "x1 + x2 <= 5\n"
            "2 x1 - x2 >= 0\n"
            "-x1 + 3 x2 >= 0\n"
            "x3 + x4 >= .5\n"
            "Bounds\n"
            "x3 >= 1.1\n"
            "x3 <= 10\n"
            "Semi-continuous\n"
            "x3\n")
        + (semi_integer? "General\nx3\n": "")
        + "End\n";
}

const double SEMI_CONTINUOUS_OPTIMUM = 8.22333;
const double SEMI_INTEGER_OPTIMUM = 8.13333;

}

/*  YAP - Yet another PWA toolkit
    Copyright 2015, Technische Universitaet Muenchen,
    Authors: Daniel Greenwald, Johannes Rauch

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef yap_MathUtilities_h
#define yap_MathUtilities_h

#include <cmath>
#include <math.h>

namespace yap {

/// \brief Math Utilities
/// \author Johannes Rauch, Daniel Greenwald
/// This code has been taken from rootpwa and modified

/// return n!
inline double factorial(int n)
{ return std::tgamma(n + 1); }

/// returns whether val is an odd number
inline bool isOdd(int val)
{ return val & 0x1; }

/// returns whether val is an even number
inline bool isEven(int val)
{ return not isOdd(val); }

/// extracts sign from value
template <typename T>
typename std::enable_if<std::is_signed<T>::value, int>::type
signum(const T& val)
{ return (T(0) < val) - (val < T(0)); }

/// optimized function for (-1)^n
inline int powMinusOne(const int exponent)
{ return isOdd(exponent) ? -1 : +1; }

}

#endif

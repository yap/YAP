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

#include <complex>
#include <type_traits>

namespace yap {

/// \brief Math Utilities
/// \author Johannes Rauch, Daniel Greenwald

/// \return whether val is an odd number
constexpr bool is_odd(int val)
{ return val & 0x1; }

/// \return whether val is an even number
constexpr bool is_even(int val)
{ return not is_odd(val); }

/// extracts sign from value
template <typename T>
typename std::enable_if<std::is_signed<T>::value, T>::type
constexpr signum(const T& val)
{ return (T(0) < val) - (val < T(0)); }

/// optimized function for (-1)^n
constexpr int pow_negative_one(int exponent)
{ return is_odd(exponent) ? -1 : +1; }

/// create imaginary number
constexpr std::complex<double> operator"" _i(unsigned long long d)
{ return std::complex<double>{0.0, static_cast<double>(d)}; }

/// create imaginary number
constexpr std::complex<double> operator"" _i(long double d)
{ return std::complex<double>{0.0, static_cast<double>(d)}; }

/// \return pi
template <typename T = double>
constexpr T pi()
{ return acos((T) - 1); }

/// convert degrees to radians
template <typename T>
constexpr T rad(const T& d)
{ return d * pi<T>() / T(180); }

/// convert radians to degrees
template <typename T>
constexpr T deg(const T& r)
{ return r * T(180) / pi<T>(); }

}

#endif

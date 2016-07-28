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

/// \file

#ifndef yap_CompensatedSum_h
#define yap_CompensatedSum_h

#include <type_traits>

namespace yap {

/// \struct CompensatedSum
/// \brief Struct to calculate a compensated sum using the Kahan summation algorithm
/// \author Johannes Rauch, Daniel Greenwald
template <typename T>
struct CompensatedSum
{
    /// sum
    T sum;

    /// correction to sum
    T correction;

    /// constructor
    CompensatedSum(T val = 0) : sum(val), correction(0) {}

    /// (implicit) cast operator
    /// \return sum
    operator T() const
    { return sum; }

    /// addition assignment operator
    CompensatedSum& operator+=(typename std::conditional<std::is_fundamental<T>::value, const T, const T&>::type value)
    {
        T y = value - correction;
        T t = sum + y;
        correction = (t - sum) - y;
        sum = t;
        return *this;
    }
};

}

#endif

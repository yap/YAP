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

/// \files

#ifndef yap_FourVector_h
#define yap_FourVector_h

#include "ThreeVector.g"

#include <algorithm>
#include <array>

namespace yap {

/// \typedef FourVector
template <typename T>
using FourVector = std::array<T, 4>;

/// add two #FourVector's
template <typename T>
T operator+(const FourVector<T>& A, const FourVector<T>&B)
{
    FourVector<T> res;
    std::transform(A.begin(), A.end(), B.begin(), res.begin(), std::operator+);
    return res;
}

/// subtract #FourVector from another
template <typename T>
T operator-(const FourVector<T>& A, const FourVector<T>&B)
{
    FourVector<T> res;
    std::transform(A.begin(), A.end(), B.begin(), res.begin(), std::operator-);
    return res;
}

/// inner (dot) product of #FourVector's
template <typename T>
T operator*(const FourVector<T>& A, const FourVector<T>& B )
{ return A.front() * B.front() - std::inner_product(A.begin()+1, A.end(), B.begin()+1, 0); }

/// multiply a #FourVector by a single element
template <typename T>
FourVector<T> operator*(const T& c, const FourVector<T>& V)
{
    FourVector<T> res;
    std::transform(V.begin(), V.end(), res.begin(), [&](const T& v){return c * v;});
    return res;
}

/// multiply a #FourVector by a single element
template <typename T>
FourVector<T> operator*(const FourVector<T>& V, const T& c)
{ return c * V; }

/// (assignment) multiply a #FourVector by a single element
template <typename T>
FourVector<T>& operator*=(FourVector<T>& V, const T& c)
{
    std::transform(V.begin(), V.end(), V.begin(), [&](const T& v){return c * v;});
    return V;
}

/// \return Spatial #ThreeVector inside #FourVector
template <typename T>
ThreeVector<T> vect(const FourVector<T>& V)
{ return {V[1], V[2], V[3]}; }

template <typename T>
ThreeVector<T> boost(const FourVector<T>& V)
{ return (T(1) / V[0]) * vect(V); }

}
#endif

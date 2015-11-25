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

#ifndef yap_NVector_h
#define yap_NVector_h

#include <algorithm>
#include <array>
#include <cassert>
#include <numeric>
#include <type_traits>

namespace yap {

/// \class NVector
/// \brief N-dimensional vector
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup VectorAlgebra
template <typename T, size_t N>
class NVector : public std::array<T, N>
{
public:
    /// Default Constructor
    NVector() : std::array<T, N>() {}

    /// Constructor
    NVector(std::initializer_list<T> list) : NVector<T, N>()
        {
            assert(list.size() == N);
            std::copy(list.begin(), list.begin() + N, this->begin());
        }
};

/// addition assignment
template <typename T, size_t N>
NVector<T, N>& operator+=(NVector<T, N>& A, const NVector<T, N>& B)
{ std::transform(A.begin(), A.end(), B.begin(), A.begin(), [](const T& a, const T& b){return a + b;}); return A; }

/// subtraction assignment
template <typename T, size_t N>
NVector<T, N>& operator-=(NVector<T, N>& A, const NVector<T, N>& B)
{ std::transform(A.begin(), A.end(), B.begin(), A.begin(), [](const T& a, const T& b){return a - b;}); return A; }

/// (assignment) multiplication by a single element
template <typename T, size_t N>
typename std::enable_if<std::is_arithmetic<T>::value, NVector<T, N>& >::type
operator*=(NVector<T, N>& A, const T& B)
{ std::transform(A.begin(), A.end(), A.begin(), [&](const T & v) {return B * v;}); return A; }

/// addition
template <typename T, size_t N>
NVector<T, N> operator+(const NVector<T, N>& A, const NVector<T, N>& B)
{ NVector<T, N> res = A; res += B; return res; }

/// subtraction
template <typename T, size_t N>
NVector<T, N> operator-(const NVector<T, N>& A, const NVector<T, N>& B)
{ NVector<T, N> res = A; A -= B; return res; }

/// inner (dot) product of #NVector's
template <typename T, size_t N>
typename std::enable_if<std::is_arithmetic<T>::value, T>::type
operator*(const NVector<T, N>& A, const NVector<T, N>& B)
{ return std::inner_product(A.begin(), A.end(), B.begin(), T(0)); }

/// multiplication: #NVector<T> * T
template <typename T, size_t N>
typename std::enable_if<std::is_arithmetic<T>::value, NVector<T, N> >::type
operator*(const NVector<T, N>& A, const T& c)
{ NVector<T, N> res = A; res *= c; return res; }

/// multiplication: T * #NVector<T>
template <typename T, size_t N>
typename std::enable_if<std::is_arithmetic<T>::value, NVector<T, N> >::type
operator*(const T& c, const NVector<T, N>& A)
{ return operator*<T, N>(A, c); }

/// unary minus
template <typename T, size_t N>
NVector<T, N> operator-(const NVector<T, N>& A)
{ NVector<T, N> res = A; std::transform(res.begin(), res.end(), res.begin(), [](const T& t){return -t;}); return res; }

}
#endif

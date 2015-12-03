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
#include <cmath>
#include <numeric>
#include <string>
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

/// \return string
template <typename T, size_t N>
typename std::enable_if<std::is_arithmetic<T>::value, std::string >::type
to_string(const NVector<T, N>& V)
{
    std::string s = "(";
    std::for_each(V.begin(), V.end(), [&](const T & t) {s += std::to_string(t) + ", ";});
    s.erase(s.size() - 2, 2);
    s += ")";
    return s;
}

/// addition assignment
template <typename T, size_t N>
NVector<T, N>& operator+=(NVector<T, N>& A, const NVector<T, N>& B)
{ std::transform(A.begin(), A.end(), B.begin(), A.begin(), [](const T & a, const T & b) {return a + b;}); return A; }

/// subtraction assignment
template <typename T, size_t N>
NVector<T, N>& operator-=(NVector<T, N>& A, const NVector<T, N>& B)
{ std::transform(A.begin(), A.end(), B.begin(), A.begin(), [](const T & a, const T & b) {return a - b;}); return A; }

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
{ NVector<T, N> res = A; res -= B; return res; }

/// \return inner (dot) product of #NVector's
template <typename T, size_t N>
typename std::enable_if<std::is_arithmetic<T>::value, T>::type
operator*(const NVector<T, N>& A, const NVector<T, N>& B)
{ return std::inner_product(A.begin(), A.end(), B.begin(), T(0)); }

/// \return squared magnitude of #NVector (using associated inner product)
template <typename T, size_t N>
typename std::enable_if<std::is_arithmetic<T>::value, T>::type
norm(const NVector<T, N>& A)
{ return A * A; }

/// \return magnitude of #NVector (using associated inner product)
template <typename T, size_t N>
typename std::enable_if<std::is_arithmetic<T>::value, T>::type
abs(const NVector<T, N>& A)
{ return sqrt(norm(A)); }

/// \return unit vector in direction of vector
/// \param V NVector to use for direction of unit vector
template <typename T, size_t N>
typename std::enable_if<std::is_arithmetic<T>::value, NVector<T, N> >::type
unit(const NVector<T, N>& V)
{ T a = abs(V); return (a == 0) ? V : (T(1) / abs(V)) * V; }

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
{ NVector<T, N> res = A; std::transform(res.begin(), res.end(), res.begin(), [](const T & t) {return -t;}); return res; }

/// \name specifically for 3D vectors
/// @{

/// \return cross product
template <typename T>
NVector<T, 3> cross(const NVector<T, 3>& A, const NVector<T, 3>& B)
{
    return {A[1]* B[2] - A[2]* B[1], A[2]* B[0] - A[0]* B[2], A[0]* B[1] - A[1]* B[0]};
}

/// \return angle between two 3D vectors
template <typename T>
typename std::enable_if<std::is_arithmetic<T>::value, T>::type
angle(const NVector<T, 3>& A, const NVector<T, 3>& B)
{ return acos(A * B / abs(A) / abs(B)); }

/// @}

}
#endif

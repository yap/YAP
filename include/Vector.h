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

#ifndef yap_Vector_h
#define yap_Vector_h

#include "Matrix.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <string>

namespace yap {

/// \class Vector
/// \brief N-dimensional column vector
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup VectorAlgebra
template <typename T, size_t N>
class Vector : public std::array<T, N>
{
public:
    /// Constructor
    constexpr Vector(const std::array<T, N>& v) noexcept : std::array<T, N>(v) {}

    /// Default constructor
    Vector() = default;

    /// Use std::array's assignment operators
    using std::array<T, N>::operator=;

    /// inner (dot) product of #Vector's
    virtual T operator*(const Vector<T, N>& B) const
    { return std::inner_product(this->begin(), this->end(), B.begin(), T(0)); }

    /// unary minus
    virtual Vector<T, N> operator-() const
    { return T(-1) * *(this); }
};

/// \return string
template <typename T, size_t N>
std::string to_string(const Vector<T, N>& V)
{
    if (N == 0)
        return "(empty)";
    std::string s = "(";
    std::for_each(V.begin(), V.end(), [&](const T & t) {s += std::to_string(t) + ", ";});
    s.erase(s.size() - 2, 2);
    s += ")";
    return s;
}

/// streamer
template <typename T, size_t N>
std::ostream& operator<<(std::ostream& os, const Vector<T, N>& V)
{ os << to_string(V); return os; }

/// addition assignment
template <typename T, size_t N>
Vector<T, N>& operator+=(Vector<T, N>& A, const Vector<T, N>& B)
{ std::transform(A.begin(), A.end(), B.begin(), A.begin(), [](const T & a, const T & b) {return a + b;}); return A; }

/// addition
template <typename T, size_t N>
Vector<T, N> operator+(const Vector<T, N>& A, const Vector<T, N>& B)
{ auto v = A; v += B; return v; }

/// subtraction assignment
template <typename T, size_t N>
Vector<T, N>& operator-=(Vector<T, N>& A, const Vector<T, N>& B)
{ std::transform(A.begin(), A.end(), B.begin(), A.begin(), [](const T & a, const T & b) {return a - b;}); return A; }

/// subtraction
template <typename T, size_t N>
Vector<T, N> operator-(const Vector<T, N>& A, const Vector<T, N>& B)
{ auto v = A; v -= B; return v; }

/// (assignment) multiplication by a single element
template <typename T, size_t N>
Vector<T, N>& operator*=(Vector<T, N>& A, const T& B)
{ std::transform(A.begin(), A.end(), A.begin(), [&](const T & v) {return B * v;}); return A; }

/// multiplication: #Vector<T> * T
template <typename T, size_t N>
Vector<T, N> operator*(const Vector<T, N>& A, const T& c)
{ auto v = A; v *= c; return v; }

/// multiplication: T * #Vector<T>
template <typename T, size_t N>
constexpr Vector<T, N> operator*(const T& c, const Vector<T, N>& A)
{ return A * c; }

/// \return squared magnitude of #Vector (using associated inner product)
template <typename T, size_t N>
constexpr T norm(const Vector<T, N>& A)
{ return A * A; }

/// \return magnitude of #Vector (using associated inner product)
template <typename T, size_t N>
constexpr T abs(const Vector<T, N>& A)
{ return sqrt(norm(A)); }

/// \return unit vector in direction of vector
/// \param V Vector to use for direction of unit vector
template <typename T, size_t N>
Vector<T, N> unit(const Vector<T, N>& V)
{ T a = abs(V); return (a == 0) ? V : (T(1) / abs(V)) * V; }

/// Matrix * Vector
template <typename T, size_t R, size_t C>
Vector<T, R> operator*(const Matrix<T, R, C>& M, const Vector<T, C>& V)
{
    Vector<T, R> v;
    for (size_t r = 0; r < R; ++r)
        for (size_t c = 0; c < C; ++c)
            v[r] += M[r][c] * V[c];
    return v;
}

/// outer product
template <typename T, size_t N>
SquareMatrix<T, N> outer(const Vector<T, N>& A, const Vector<T, N>& B)
{
    SquareMatrix<T, N> m;
    for (size_t r = 0; r < N; ++r)
        for (size_t c = 0; c < N; ++r)
            m[r][c] = A[r] * B[c];
    return m;
}

}
#endif

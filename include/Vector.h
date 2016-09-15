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

#include "fwd/Vector.h"

#include "Matrix.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <numeric>
#include <ostream>
#include <string>
#include <type_traits>

namespace yap {

/// \class VectorIterator
/// \brief Iterator for Vector class
/// \author Daniel Greenwald
/// \ingroup VectorAlgebra
template <typename T, size_t N>
class VectorIterator : public std::iterator<std::random_access_iterator_tag, T, typename std::array<T, N>::difference_type>
{
public:

    /// Constructor
    VectorIterator(typename std::array<T, N>::iterator it)
        : Iterator_(it)
    {}

    /// inequality operator
    bool operator!=(VectorIterator b) const
    { return this->Iterator_ != b.Iterator_; }

    /// dereference operator
    T& operator*()
    { return *(this->Iterator_); }

    /// pointer operator
    T* operator->()
    { return std::array<T, N>::operator->(this->Iterator_); }

    /// pre-increment operator
    VectorIterator& operator++()
    { ++(this->Iterator_); return *this; }

    /// post-increment operator
    VectorIterator operator++(int)
    { VectorIterator it(*this); (this->Iterator_)++; return it; }

    /// pre-decrement operator
    VectorIterator& operator--()
    { --(this->Iterator_); return *this; }

    /// post-decrement operator
    VectorIterator operator--(int)
    { VectorIterator it(*this); --(this->Iterator_); return it; }

    /// addition assignment operator
    VectorIterator& operator+=(typename VectorIterator::difference_type n)
    { this->Iterator_ += n; return *this; }

    /// subtraction assignment operator
    VectorIterator& operator-=(typename VectorIterator::difference_type n)
    { this->Iterator_ -= n; return *this; }

    /// difference operator
    friend const typename VectorIterator::difference_type operator-(VectorIterator lhs, VectorIterator rhs)
    { return lhs.Iterator_ - rhs.Iterator_; }

    /// access operator
    VectorIterator operator[](typename VectorIterator::difference_type n) const
    { VectorIterator it(Iterator_[n]); return it; }

    /// less-than operator
    friend const bool operator<(VectorIterator lhs, VectorIterator rhs)
    { return lhs.Iterator_ < rhs.Iterator_; }

    /// greater-than operator
    friend bool operator>(VectorIterator lhs, VectorIterator rhs)
    { return lhs.Iterator_ > rhs.Iterator_; }

private:

    /// internal iterator
    typename std::array<T, N>::iterator Iterator_;
};

/// addition operator
template <typename T, size_t N>
const VectorIterator<T, N> operator+(VectorIterator<T, N> a, typename VectorIterator<T, N>::difference_type n)
{ return (VectorIterator<T, N>(a) += n); }

/// addition operator
template <typename T, size_t N>
const VectorIterator<T, N> operator+(typename VectorIterator<T, N>::difference_type n, VectorIterator<T, N> a)
{ return a + n; }

/// subtraction operator
template <typename T, size_t N>
const VectorIterator<T, N> operator-(VectorIterator<T, N> a, typename VectorIterator<T, N>::difference_type n)
{ return (VectorIterator<T, N>(a) -= n); }

/// greater-than-or-equal operator
template <typename T, size_t N>
bool operator>=(VectorIterator<T, N> a, VectorIterator<T, N> b)
{ return !(a < b); }

/// less-than-or-equal operator
template <typename T, size_t N>
bool operator<=(VectorIterator<T, N> a, VectorIterator<T, N> b)
{ return !(b > a); }

/// \class Vector
/// \brief N-dimensional column vector
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup VectorAlgebra
template <typename T, size_t N>
class Vector
{
public:
    /// Constructor
    constexpr Vector(const std::array<T, N>& v) noexcept : Elements_(v) {}

    /// Default constructor
    Vector(T element = 0)
    { Elements_.fill(element); }

    /// \return size
    constexpr size_t size() const
    { return N; }

    /// element access operator
    T& operator[](size_t i)
    { return Elements_[i]; }

    /// element access operator
    constexpr T operator[](size_t i) const
    { return Elements_[i]; }

    /// access to front
    T& front()
    { return Elements_.front(); }

    /// access to front
    const T& front() const
    { return Elements_.front(); }

    /// access to begin
    VectorIterator<T, N> begin()
    { return VectorIterator<T, N>(Elements_.begin()); }

    /// access to begin
    VectorIterator<const T, N> begin() const
    { return VectorIterator<const T, N>(Elements_.begin()); }

    /// access to end
    VectorIterator<T, N> end()
    { return VectorIterator<T, N>(Elements_.end()); }

    /// access to end
    VectorIterator<const T, N> end() const
    { return VectorIterator<const T, N>(Elements_.end()); }

    /// inner (dot) product of #Vector's
    virtual auto operator*(const Vector<T, N>& B) const
    -> const decltype(T(0) * T(0))
    { return std::inner_product(Elements_.begin(), Elements_.end(), B.Elements_.begin(), T(0) * T(0)); }

    /// equality operator
    friend constexpr bool operator==(const Vector<T, N>& lhs, const Vector<T, N>& rhs)
    { return lhs.Elements_ == rhs.Elements_; }

    /// Calculate the angle between two vectors
    friend const T angle(const Vector<T, N>& A, const Vector<T, N>& B)
    {
        T arg = A * B / abs(A) / abs(B);

        // correct for arg just outside boundary (by numerical precision)
        if (std::isfinite(arg) and fabs(arg) > 1)
            arg = (arg > 0) ? (T)1 : (T) - 1;

        return acos(arg);
    }

    /// \return squared magnitude of #Vector (using associated inner product)
    friend constexpr T norm(const Vector<T, N>& A)
    { return A * A; }

private:

    /// internal storage
    std::array<T, N> Elements_;
};

/// \return string
template <typename T, size_t N>
std::string to_string(const Vector<T, N>& V)
{
    return "(" + (
        (N == 0) ?
        "empty"
        : 
        std::accumulate(V.begin(), V.end(), std::string(""),
                        [](std::string& s, typename std::conditional<std::is_fundamental<T>::value, T, const T&>::type t)
                        { using std::to_string; return s += ", " + to_string(t);}).erase(0, 2)
        ) + ")";
}

/// streamer
template <typename T, size_t N>
std::ostream& operator<<(std::ostream& os, const Vector<T, N>& V)
{ os << to_string(V); return os; }

/// perform unary op on vector
template <typename T, size_t N, class UnaryOp>
Vector<T, N>& unary_op(Vector<T, N>& V, UnaryOp op1)
{ std::transform(V.begin(), V.end(), V.begin(), op1); return V; }

/// perform binary op on two vectors, storing result in first vector
template <typename T, size_t N, class BinaryOp>
Vector<T, N>& binary_op(Vector<T, N>& lhs, const Vector<T, N>& rhs, BinaryOp op2)
{ std::transform(lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(), op2); return lhs; }

/// negation
template <typename T, size_t N>
const Vector<T, N> operator-(Vector<T, N> V)
{ return unary_op(V, std::negate<T>()); }

/// Vector addition assignment
template <typename T, size_t N>
Vector<T, N>& operator+=(Vector<T, N>& lhs, const Vector<T, N>& rhs)
{ return binary_op(lhs, rhs, std::plus<T>()); }

/// Vector addition
template <typename T, size_t N>
const Vector<T, N> operator+(Vector<T, N> A, const Vector<T, N>& B)
{ return A += B; }

/// Vector subtraction assignment
template <typename T, size_t N>
Vector<T, N>& operator-=(Vector<T, N>& lhs, const Vector<T, N>& rhs)
{ return binary_op(lhs, rhs, std::minus<T>()); }

/// Vector subtraction
template <typename T, size_t N>
const Vector<T, N> operator-(Vector<T, N> A, const Vector<T, N>& B)
{ return A -= B; }

/// Vector scalar multiplication assignment
template <typename T, size_t N>
Vector<T, N>& operator*=(Vector<T, N>& V, typename std::conditional<std::is_fundamental<T>::value, T, const T&>::type c)
{return unary_op(V, std::bind(std::multiplies<T>(), c, std::placeholders::_1)); }

/// Vector scalar multiplication
template <typename T, size_t N>
const Vector<T, N> operator*(Vector<T, N> V, typename std::conditional<std::is_fundamental<T>::value, T, const T&>::type c)
{ return V *= c; }

/// Vector scalar multiplication
template <typename T, size_t N>
const Vector<T, N> operator*(typename std::conditional<std::is_fundamental<T>::value, T, const T&>::type c, Vector<T, N> V)
{ return V *= c; }

/// Vector scalar division assignment
template <typename T, size_t N>
Vector<T, N>& operator/=(Vector<T, N>& V, typename std::conditional<std::is_fundamental<T>::value, T, const T&>::type c)
{return unary_op(V, std::bind(std::divides<T>(), std::placeholders::_1, c)); }

/// Vector scalar division
template <typename T, size_t N>
const Vector<T, N> operator/(Vector<T, N> V, typename std::conditional<std::is_fundamental<T>::value, T, const T&>::type c)
{ return V /= c; }

/// \return magnitude of #Vector (using associated inner product)
template <typename T, size_t N>
constexpr T abs(const Vector<T, N>& A)
{ return sqrt(norm(A)); }

/// \return unit vector in direction of vector
/// \param V Vector to use for direction of unit vector
template <typename T, size_t N>
const Vector<T, N> unit(const Vector<T, N>& V)
{ T a = abs(V); return (a == 0) ? V : V / a; }

/// Matrix * Vector
template <typename T, size_t R, size_t C>
const Vector<T, R> operator*(const Matrix<T, R, C>& M, const Vector<T, C>& V)
{
    Vector<T, R> v;
    for (size_t r = 0; r < R; ++r)
        for (size_t c = 0; c < C; ++c)
            v[r] += M[r][c] * V[c];
    return v;
}

/// Matrix * Vector with different template types
template <typename T1, typename T2, size_t R, size_t C>
auto operator*(const Matrix<T1, R, C>& M, const Vector<T2, C>& V)
-> const Vector<typename std::remove_cv<decltype(operator*(M[0][0], V[0]))>::type, R>
{
    Vector<typename std::remove_cv<decltype(operator*(M[0][0], V[0]))>::type, R> v;
    for (size_t r = 0; r < R; ++r)
        for (size_t c = 0; c < C; ++c)
            v[r] += M[r][c] * V[c];
    return v;
}

/// outer product
template <typename T, size_t N>
const SquareMatrix<T, N> outer(const Vector<T, N>& A, const Vector<T, N>& B)
{
    SquareMatrix<T, N> m;
    for (size_t r = 0; r < N; ++r)
        for (size_t c = 0; c < N; ++c)
            m[r][c] = A[r] * B[c];
    return m;
}

}
#endif

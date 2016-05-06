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
    { VectorIterator it(*this); (this->Iterator_)--; return it; }

    /// addition assignment operator
    VectorIterator& operator+=(typename VectorIterator::difference_type n)
    { this->Iterator_ += n; return *this; }

    /// subtraction assignment operator
    VectorIterator& operator-=(typename VectorIterator::difference_type n)
    { this->Iterator_ -= n; return *this; }

    /// difference operator
    typename VectorIterator::difference_type operator-(VectorIterator a) const
    { return this->Iterator_ - a.Iterator_; }

    /// access operator
    VectorIterator operator[](typename VectorIterator::difference_type n) const
    { VectorIterator it(Iterator_[n]); return it; }

    /// less-than operator
    bool operator<(VectorIterator b) const
    { return this->Iterator_ < b.Iterator_; }

    /// greater-than operator
    bool operator>(VectorIterator b) const
    { return this->Iterator_ > b.Iterator_; }

private:

    /// internal iterator
    typename std::array<T, N>::iterator Iterator_;
};

/// addition operator
template <typename T, size_t N>
VectorIterator<T, N> operator+(VectorIterator<T, N> a, typename VectorIterator<T, N>::difference_type n)
{ VectorIterator<T, N> it(a); it += n; return it; }

/// addition operator
template <typename T, size_t N>
VectorIterator<T, N> operator+(typename VectorIterator<T, N>::difference_type n, VectorIterator<T, N> a)
{ return a + n; }

/// subtraction operator
template <typename T, size_t N>
VectorIterator<T, N> operator-(VectorIterator<T, N> a, typename VectorIterator<T, N>::difference_type n)
{ VectorIterator<T, N> it(a); it -= n; return it; }

/// greater-than-or-equal operator
template <typename T, size_t N>
bool operator>=(VectorIterator<T, N> a, VectorIterator<T, N> b)
{ return !(a < b); }

/// less-than-or-equal operator
template <typename T, size_t N>
bool operator<=(VectorIterator<T, N> a, VectorIterator<T, N> b)
{ return !(b > a); }

/// \class VectorExpression
/// \brief Class for holding vector expressions
/// \author Daniel Greenwald
/// \ingroup VectorAlgebra
template <typename T, size_t N, typename E>
class VectorExpression
{
public:

    /// access operator
    constexpr T operator[](size_t i) const
    { return static_cast<const E&>(*this)[i]; }

    /// \return size
    constexpr size_t size() const
    { return N; }

    /// cast to E
    operator E& ()
    { return static_cast<E&>(*this); }

    /// cast to E
    constexpr operator const E& () const
    { return static_cast<const E&>(*this); }

};

/// \class Vector
/// \brief N-dimensional column vector
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup VectorAlgebra
template <typename T, size_t N>
class Vector : public VectorExpression<T, N, Vector<T, N> >
{
public:
    /// Constructor
    constexpr Vector(const std::array<T, N>& v) noexcept : Elements_(v) {}

    /// expression constructor
    template <typename E>
    Vector(const VectorExpression<T, N, E>& V)
    { for (size_t i = 0; i < V.size(); ++i) Elements_[i] = V[i]; }

    /// Default constructor
    Vector()
    { Elements_.fill(T(0)); }

    /// element access operator
    T& operator[](size_t i)
    { return Elements_[i]; }

    /// element access operator
    const T operator[](size_t i) const
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
    virtual T operator*(const Vector<T, N>& B) const
    { return std::inner_product(Elements_.begin(), Elements_.end(), B.Elements_.begin(), T(0)); }

    /// equality operator
    friend bool operator==(const Vector<T, N>& lhs, const Vector<T, N>& rhs)
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

/// \class VectorAddition
/// \brief Expression for addition of two VectorExpressions
/// \author Daniel Greenwald
/// \ingroup VectorAlgebra
template <typename T, size_t N, typename E1, typename E2>
class VectorAddition : public VectorExpression<T, N, VectorAddition<T, N, E1, E2> >
{
public:
    /// Constructor
    constexpr VectorAddition(const VectorExpression<T, N, E1>& a, const VectorExpression<T, N, E2>& b)
        : VectorExpression<T, N, VectorAddition<T, N, E1, E2> >(), A_(a), B_(b)
    {}

    /// access operator
    constexpr T operator[](size_t i) const
    { return A_[i] + B_[i]; }

private:

    /// lhs expression
    const E1& A_;

    /// rhs expression
    const E2& B_;

};

/// addition of two vectors
template <typename T, size_t N, typename E1, typename E2>
constexpr VectorAddition<T, N, E1, E2> operator+(const VectorExpression<T, N, E1>& a, const VectorExpression<T, N, E2>& b)
{ return VectorAddition<T, N, E1, E2>(a, b); }

/// \class VectorSubtraction
/// \brief Expression for subtraction of two VectorExpressions
/// \author Daniel Greenwald
/// \ingroup VectorAlgebra
template <typename T, size_t N, typename E1, typename E2>
class VectorSubtraction : public VectorExpression<T, N, VectorSubtraction<T, N, E1, E2> >
{
public:
    /// Constructor
    constexpr VectorSubtraction(const VectorExpression<T, N, E1>& a, const VectorExpression<T, N, E2>& b)
        : VectorExpression<T, N, VectorSubtraction<T, N, E1, E2> >(), A_(a), B_(b)
    {}

    /// access operator
    constexpr T operator[](size_t i) const
    { return A_[i] - B_[i]; }

private:

    /// lhs expression
    const E1& A_;

    /// rhs expression
    const E2& B_;

};

/// subtraction of two vectors
template <typename T, size_t N, typename E1, typename E2>
constexpr VectorSubtraction<T, N, E1, E2> operator-(const VectorExpression<T, N, E1>& a, const VectorExpression<T, N, E2>& b)
{ return VectorSubtraction<T, N, E1, E2>(a, b); }

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

/// unary minus
template <typename T, size_t N>
constexpr Vector<T, N> operator-(const Vector<T, N>& V)
{ return T(-1) * V; }

/// subtraction assignment
template <typename T, size_t N>
Vector<T, N>& operator-=(Vector<T, N>& A, const Vector<T, N>& B)
{ std::transform(A.begin(), A.end(), B.begin(), A.begin(), [](const T & a, const T & b) {return a - b;}); return A; }

/// (assignment) multiplication by a single element
template <typename T, size_t N>
Vector<T, N>& operator*=(Vector<T, N>& A, const T& B)
{ std::transform(A.begin(), A.end(), A.begin(), [&](const T & v) {return B * v;}); return A; }

/// multiplication: #Vector<T> * T
template <typename T, size_t N>
const Vector<T, N> operator*(const Vector<T, N>& A, const T& c)
{ auto v = A; v *= c; return v; }

/// multiplication: T * #Vector<T>
template <typename T, size_t N>
constexpr Vector<T, N> operator*(const T& c, const Vector<T, N>& A)
{ return A * c; }

/// (assignment) division by a single element
template <typename T, size_t N>
Vector<T, N>& operator/=(Vector<T, N>& A, const T& B)
{ std::transform(A.begin(), A.end(), A.begin(), [&](const T & v) {return v / B;}); return A; }

/// division: #Vector<T> / T
template <typename T, size_t N>
const Vector<T, N> operator/(const Vector<T, N>& A, const T& c)
{ auto v = A; v /= c; return v; }

/// \return magnitude of #Vector (using associated inner product)
template <typename T, size_t N>
constexpr T abs(const Vector<T, N>& A)
{ return sqrt(norm(A)); }

/// \return unit vector in direction of vector
/// \param V Vector to use for direction of unit vector
template <typename T, size_t N>
const Vector<T, N> unit(const Vector<T, N>& V)
{ T a = abs(V); return (a == 0) ? V : (T(1) / abs(V)) * V; }

/// Matrix * Vector
template <typename T, size_t R, size_t C>
const Vector<T, R> operator*(const Matrix<T, R, C>& M, const Vector<T, C>& V)
{
    Vector<T, R> v = {};
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

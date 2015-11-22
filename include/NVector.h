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
#include <numeric>

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
    NVector(std::initializer_list<T> list) : std::array<T, N>()
    { std::copy(list.begin(), list.begin() + N, this->begin()); }

    /// Destructor
    virtual ~NVector() {}

    /// addition assignment
    NVector<T, N>& operator+=(const NVector<T, N>& B)
    { std::transform(this->begin(), this->end(), B.begin(), this->begin(), std::operator+); return *this; }

    /// subtraction assignment
    NVector<T, N>& operator-=(const NVector<T, N>& B)
    { std::transform(this->begin(), this->end(), B.begin(), this->begin(), std::operator-); return *this; }

    /// (assignment) multiplication by a single element
    NVector<T, N>& operator*=(const T& B)
    { std::transform(this->begin(), this->end(), this->begin(), [&](const T & v) {return B * v;}); return *this; }

    /// addition
    NVector<T, N> operator+(const NVector<T, N>& B) const
    {
        NVector<T, N> res;
        std::transform(this->begin(), this-> end(), B.begin(), res.begin(), std::operator+);
        return res;
    }

    /// subtraction
    NVector<T, N> operator-(const NVector<T, N>& B) const
    {
        NVector<T, N> res;
        std::transform(this->begin(), this->end(), B.begin(), res.begin(), std::operator-);
        return res;
    }

    /// inner (dot) product of #NVector's
    virtual T operator*(const NVector<T, N>& B) const
    { return std::inner_product(this->begin(), this->end(), B.begin(), 0); }

    /// multiplication: #NVector<T> * T
    NVector<T, N> operator*(const T& c) const
    {
        NVector<T, N> res;
        std::transform(this->begin(), this->end(), res.begin(), [&](const T & t) {return t * c;});
        return res;
    }
};

/// multiplication: T * #NVector<T>
template <typename T, size_t N>
NVector<T, N> operator*(const T& c, const NVector<T, N>& V)
{ return V * c; }

}
#endif

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

#ifndef yap_IntegralElement_h
#define yap_IntegralElement_h

#include "fwd/IntegralElement.h"

namespace yap {

/// \class IntegralElement
/// \brief Holds the values of a component of an integral
/// \author Daniel Greenwald
/// \ingroup Integration
template <typename T>
struct IntegralElement {

    /// integral value
    T value;

    /// constructor
    /// \param val initial value of integral component
    IntegralElement(T val = 0) : value(val) {}

    /// addition assignment operator
    IntegralElement& operator+=(const IntegralElement& B)
    { value += B.value; return *this; }

};

/// \return addition of two IntegralElements
template <typename T>
inline const IntegralElement<T> operator+(IntegralElement<T> A, const IntegralElement<T>& B)
{ return A += B; }

/// multiplication operator
template <typename T, typename U>
inline const IntegralElement<T> operator*(const IntegralElement<T>& A, const U& B)
{ return IntegralElement<T>(A.value * B); }

/// multiplication operator
template <typename T, typename U>
inline const IntegralElement<T> operator*(const U& A, const IntegralElement<T>& B)
{ return operator*(B, A); }

/// \return string of IntegralElement
template <typename T>
inline const std::string to_string(const IntegralElement<T>& a)
{ using std::to_string; return to_string(a.value); }

}

#endif

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

#include "complex_to_string.h"

namespace yap {

/// \class IntegralElement
/// \brief Holds the values of a component of an integral
/// \author Daniel Greenwald
/// \ingroup Integration
template <typename T>
class IntegralElement
{
public:

    /// constructor
    /// \param val initial value of integral component
    explicit IntegralElement(T val = 0) : Value_(val) {}

    /// reset
    void reset()
    { Value_ = 0; }

    /// access value
    const T value() const
    { return Value_; }

    /// access value
    T& value()
    { return Value_; }

    /// multiplication (by T) assignment operator
    template <typename U>
    IntegralElement& operator*=(const U& rhs)
    { Value_ *= rhs; return *this; }

    /// addition assignment operator
    IntegralElement& operator+=(const IntegralElement& B)
    { Value_ += B.Value_; return *this; }

    /// unary minus operator
    IntegralElement& operator-() const
    { Value_ *= -1; return *this; }

    /// subtraction assignment operator
    IntegralElement& operator-=(const IntegralElement& B)
    { return *this += -B; }

    /// multiplication assignment operator
    IntegralElement& operator*=(const IntegralElement& B)
    { Value_ *= B.Value_; return *this; }

    /// division assignment operator
    IntegralElement& operator/=(const IntegralElement& B)
    { Value_ /= B.Value_; return *this; }

    /// multiplication operator
    template <typename U>
    friend const IntegralElement operator*(IntegralElement A, const U& B)
    { return A *= B; }

    /// multiplication operator
    template <typename U>
    friend const IntegralElement operator*(const U& A, IntegralElement B)
    { return operator*(B, A); }

    /// access value by cast
    explicit operator T() const
    { return value; }

    /// access real value as complex value by cast
    template <typename U>
    explicit operator std::complex<U>() const
    { return static_cast<std::complex<U> >(Value_); }

    /// cast into complex integral element, from real one
    template <typename U>
    explicit operator IntegralElement<std::complex<U> >() const
    { return IntegralElement<std::complex<U> >(static_cast<std::complex<U> >(Value_)); }

private:

    /// integral value
    T Value_;

};

/// \return addition of two IntegralElements
template <typename T>
inline const IntegralElement<T> operator+(IntegralElement<T> A, const IntegralElement<T>& B)
{ return A += B; }

/// \return subtraction of two IntegralElements
template <typename T>
inline const IntegralElement<T> operator-(IntegralElement<T> A, const IntegralElement<T>& B)
{ return A -= B; }

/// \return multiplication of two IntegralElements
template <typename T>
inline const IntegralElement<T> operator*(IntegralElement<T> A, const IntegralElement<T>& B)
{ return A *= B; }

/// \return division IntegralElements into another
template <typename T>
inline const IntegralElement<T> operator/(IntegralElement<T> A, const IntegralElement<T>& B)
{ return A /= B; }

/// \return string of IntegralElement
template <typename T>
inline std::string to_string(const IntegralElement<T>& a)
{ using std::to_string; return to_string(a.value()); }

/// \name operations on ComplexIntegralElement
/// @{

/// \return conjugate of element
inline ComplexIntegralElement conj(const ComplexIntegralElement& Z)
{ return ComplexIntegralElement(conj(Z.value())); }

/// \return RealIntegralElement for real component of ComplexIntegralElement
inline RealIntegralElement real(const ComplexIntegralElement& Z)
{ return RealIntegralElement(real(Z.value())); }

/// \return RealIntegralElement for imaginary component of ComplexIntegralElement
inline RealIntegralElement imag(const ComplexIntegralElement& Z)
{ return RealIntegralElement(imag(Z.value())); }

/// @}


}

#endif

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

#ifndef yap_ComplexBasis_h
#define yap_ComplexBasis_h

#include "fwd/ComplexBasis.h"

#include "Matrix.h"
#include "Vector.h"

#include <complex>
#include <math>

namespace yap {

namespace complexBasis {

/// \class cartesian
/// complex number representation in cartesian coordinates, with covariance
template <typename T>
class cartesian {
public:
    /// constructor
    cartesian(const& std::complex<T>& val, const SquareMatrix<T, 2>& cov = SquareMatrix<T, 2>()) :
        value_(val), covariance_(cov)
    {}

    /// conversion constructor
    cartesian(const polar& pol) :
        value_( pol.r()*cos(pol.phi()), pol.r()*sin(pol.phi()) ),
        covariance_(jacobian(pol)*pol.covariance()*jacobian(pol).transpose())
    {}

    /// \return real part
    const T real() const
    { return real(value_); }

    /// \return imaginary part
    const T imag() const
    { return imag(value_); }

    /// casting operator to std::complex
    operator std::complex<T>&()
    { return value_; }

    /// \return covariance
    const SquareMatrix<T, 2>& covariance() const
    { return covariance_; }

private:
    /// jacobian for transformation from polar to cartesian
    SquareMatrix<T, 2> jacobian(const polar& pol) const
    { return SquareMatrix<T, 2>({cos(pol.phi()), -pol.r()*sin(pol.phi()),
                                 sin(pol.phi()),  pol.r()*cos(pol.phi()) }); }

    std::complex<T> value_;
    SquareMatrix<T, 2> covariance_;
};

/// \class polar
/// complex number representation in polar coordinates, with covariance
template <typename T>
class polar {
public:
    /// constructor
    polar(T r, T phi, const SquareMatrix<T, 2>& cov = SquareMatrix<T, 2>()) :
        value_({r, phi}), covariance_(cov)
    {}

    /// conversion constructor
    polar(const cartesian& cart) :
        value_( sqrt(pow(cart.real(), 2) + pow(cart.imag(), 2)), atan2(cart.imag(), cart.real()) ),
        covariance_(jacobian(cart)*cart.covariance()*jacobian(cart).transpose())
    {}

    /// \return radius
    const T r() const
    { return value_[0]; }

    /// \return angle
    const T phi() const
    { return value_[1]; }

    /// \return covariance
    const SquareMatrix<T, 2>& covariance() const
    { return covariance_; }

private:
    /// jacobian for transformation from cartesian to polar
    SquareMatrix<T, 2> jacobian(const cartesian& cart) const
    { return SquareMatrix<T, 2>({ cart.real()/r(),         cart.imag()/r,
                                 -cart.imag()/pow(r(), 2), cart.real()/pow(r(), 2) }); }

    /// radius, polar angle [rad]
    Vector<T, 2> value_;
    SquareMatrix<T, 2> covariance_;
};

}

}

#endif

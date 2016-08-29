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
#include <math.h>

namespace yap {

namespace complex_basis {

/// \typedef covariance_matrix
template <typename T>
using covariance_matrix = SquareMatrix<T, 2>;


namespace jacobian {

/// jacobian from polar to cartesian
template <typename T>
SquareMatrix<T, 2> polar_to_cartesian(const polar<T>& pol)
{
    std::complex<T> v(pol);
    return SquareMatrix<T, 2>({cos(std::arg(v)), -std::abs(v)*sin(std::arg(v)),
                               sin(std::arg(v)),  std::abs(v)*cos(std::arg(v)) });
}

/// jacobian for transformation from cartesian to polar
template <typename T>
SquareMatrix<T, 2> cartesian_to_polar(const cartesian<T>& cart)
{
    std::complex<T> v(cart);
    return SquareMatrix<T, 2>({ std::real(v)/std::abs(v),  std::imag(v)/std::abs(v),
                               -std::imag(v)/std::norm(v), std::real(v)/std::norm(v) });
}

}

/// \class basis
/// complex number representation
template <typename T>
class basis {
protected:
    /// constructor
    basis(const Vector<T, 2>& val, const covariance_matrix<T>& cov) :
        Value_(val), Covariance_(cov)
    {}

public:
    /// casting operator to std::complex
    virtual operator std::complex<T>() const = 0;

    /// \return value
    const Vector<T, 2>& value() const
    { return Value_; }

    /// \return covariance
    const covariance_matrix<T>& covariance() const
    { return Covariance_; }

private:
    const Vector<T, 2> Value_;
    const covariance_matrix<T> Covariance_;
};

/// \class cartesian
/// complex number representation in cartesian coordinates, with covariance
template <typename T>
class cartesian : public basis<T> {
public:
    /// constructor
    /// \param re real
    /// \param im imaginary
    /// \param cov covariance
    explicit cartesian(T re, T im, const covariance_matrix<T>& cov = covariance_matrix<T>()) :
        basis<T>(Vector<T, 2>({re, im}), cov)
    {}

    /// constructor
    /// \param val real and imaginary
    /// \param cov covariance
    explicit cartesian(const Vector<T, 2>& val, const covariance_matrix<T>& cov = covariance_matrix<T>()) :
        basis<T>(val, cov)
    {}

    /// constructor
    /// \param val complex value
    /// \param cov covariance
    explicit cartesian(const std::complex<T>& val, const covariance_matrix<T>& cov = covariance_matrix<T>()) :
        cartesian<T>(std::real(val), std::imag(val), cov)
    {}

    /// constructor
    /// \param val complex value
    /// \param cov_elements list of covariance entries, starting with the diagonal, 1st off diagonal ...
    explicit cartesian(const std::complex<T>& val, std::initializer_list<T> cov_elements) :
        cartesian<T>(val, symmetricMatrix<T, 2>(std::forward<std::initializer_list<T>>(cov_elements)))
    {}

    /// conversion constructor
    /// \param pol complex number in polar representation
    cartesian(const polar<T>& pol) :
        cartesian<T>(std::complex<T>(pol),
                jacobian::polar_to_cartesian(pol) * pol.covariance() * transpose(jacobian::polar_to_cartesian(pol)))
    {}

    /// casting operator to std::complex
    virtual operator std::complex<T>() const override
    { return std::complex<T>(basis<T>::value()[0], basis<T>::value()[1]); }
};

/// \class polar
/// complex number representation in polar coordinates, with covariance
template <typename T>
class polar : public basis<T> {
public:
    /// constructor
    /// \param r radius
    /// \param phi angle
    /// \param cov covariance
    explicit polar(T r, T phi, const covariance_matrix<T>& cov = covariance_matrix<T>()) :
        basis<T>(Vector<T, 2>({r, phi}), cov)
    {}

    /// constructor
    /// \param val radius and angle
    /// \param cov covariance
    explicit polar(const Vector<T, 2>& val, const covariance_matrix<T>& cov = covariance_matrix<T>()) :
        basis<T>(val, cov)
    {}

    /// constructor
    /// \param r radius
    /// \param phi angle
    /// \param cov_elements list of covariance entries, starting with the diagonal, 1st off diagonal ...
    explicit polar(T r, T phi, std::initializer_list<T> cov_elements) :
        polar<T>(r, phi,
                 symmetricMatrix<T, 2>(std::forward<std::initializer_list<T>>(cov_elements)))
    {}

    /// conversion constructor
    /// \param pol complex number in cartesian representation
    polar(const cartesian<T>& cart) :
        polar<T>(std::abs(std::complex<T>(cart)), std::arg(std::complex<T>(cart)),
                jacobian::cartesian_to_polar(cart) * cart.covariance() * transpose(jacobian::cartesian_to_polar(cart)))
    {}

    /// casting operator to std::complex
    virtual operator std::complex<T>() const override
    { return std::polar<T>(basis<T>::value()[0], basis<T>::value()[1]); }
};


/// convert to string
template <typename T>
inline std::string to_string(const cartesian<T>& z)
{
    std::string s = to_string(z.value()[0]) + " + i * " + to_string(z.value()[1]);
    s += "covariance = " + to_string(z.covariance());
    return s;
}

/// convert to string
template <typename T>
inline std::string to_string(const polar<T>& z)
{
    std::string s = to_string(z.value()[0]) + " * exp(i * " + to_string(z.value()[1]) + ")";
    s += "covariance = " + to_string(z.covariance());
    return s;
}

}

}

#endif

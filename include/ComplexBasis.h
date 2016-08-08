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

/// \typedef covariance_type
template <typename T>
using covariance_type = SquareMatrix<T, 2>;

/// \class basis
/// complex number representation
template <typename T>
class basis {
public:
    /// casting operator to std::complex
    virtual operator std::complex<T>() const = 0;

    /// \return value
    const Vector<T, 2>& value() const
    { return value_; }

    /// \return covariance
    const SquareMatrix<T, 2>& covariance() const
    { return covariance_; }

    /// \return names of basis as string
    virtual std::string basisString() const = 0;

protected:
    /// constructor
    basis(const Vector<T, 2>& val, const SquareMatrix<T, 2>& cov = SquareMatrix<T, 2>()) :
        value_(val), covariance_(cov)
    {}

    /// constructor
    basis(T v1, T v2, const SquareMatrix<T, 2>& cov = SquareMatrix<T, 2>()) :
        value_({v1, v2}), covariance_(cov)
    {}

    const Vector<T, 2> value_;

private:
    const covariance_type<T> covariance_;
};

/// \class cartesian
/// complex number representation in cartesian coordinates, with covariance
template <typename T>
class cartesian : public basis<T> {
public:
    /// constructor
    /// \param v1 real
    /// \param v2 imaginary
    /// \param cov covariance
    explicit cartesian(T v1, T v2, const SquareMatrix<T, 2>& cov = SquareMatrix<T, 2>()) :
        basis<T>(v1, v2, cov)
    {}

    /// constructor
    /// \param val real and imaginary
    /// \param cov covariance
    explicit cartesian(const Vector<T, 2>& val, const SquareMatrix<T, 2>& cov = SquareMatrix<T, 2>()) :
        basis<T>(val, cov)
    {}

    /// constructor
    /// \param val complex value
    /// \param cov covariance
    explicit cartesian(const std::complex<T>& val, const SquareMatrix<T, 2>& cov = SquareMatrix<T, 2>()) :
        basis<T>(std::real(val), std::imag(val))
    {}

    /// constructor
    /// \param val complex value
    /// \param cov_elements list of covariance entries, starting with the diagonal, 1st off diagonal ...
    explicit cartesian(const std::complex<T>& val, std::initializer_list<T> cov_elements) :
        basis<T>(std::real(val), std::imag(val),
                 symmetricMatrix<T, 2>(std::forward<std::initializer_list<T>>(cov_elements)))
    {}

    /// conversion constructor
    /// \param pol complex number in polar representation
    cartesian(const polar<T>& pol) :
        basis<T>(std::real(std::complex<T>(pol)), std::imag(std::complex<T>(pol)),
                 jac_p_c(pol) * pol.covariance() * transpose(jac_p_c(pol)))
    {}

    /// casting operator to std::complex
    virtual operator std::complex<T>() const override
    { return std::complex<T>(basis<T>::value_[0], basis<T>::value_[1]); }

    virtual std::string basisString() const override
    { return "(real, imag)"; }
};

/// \class polar
/// complex number representation in polar coordinates, with covariance
template <typename T>
class polar : public basis<T> {
public:
    /// constructor
    /// \param v1 radius
    /// \param v2 angle
    /// \param cov covariance
    explicit polar(T v1, T v2, const SquareMatrix<T, 2>& cov = SquareMatrix<T, 2>()) :
        basis<T>(v1, v2, cov)
    {}

    /// constructor
    /// \param val radius and angle
    /// \param cov covariance
    explicit polar(const Vector<T, 2>& val, const SquareMatrix<T, 2>& cov = SquareMatrix<T, 2>()) :
        basis<T>(val, cov)
    {}

    /// constructor
    /// \param v1 radius
    /// \param v2 angle
    /// \param cov_elements list of covariance entries, starting with the diagonal, 1st off diagonal ...
    explicit polar(T v1, T v2, std::initializer_list<T> cov_elements) :
        basis<T>(v1, v2,
                 symmetricMatrix<T, 2>(std::forward<std::initializer_list<T>>(cov_elements)))
    {}

    /// conversion constructor
    /// \param pol complex number in cartesian representation
    polar(const cartesian<T>& cart) :
        basis<T>(std::abs(std::complex<T>(cart)), std::arg(std::complex<T>(cart)),
                 jac_c_p(cart) * cart.covariance() * transpose(jac_c_p(cart)))
    {}

    /// casting operator to std::complex
    virtual operator std::complex<T>() const override
    { return std::polar<T>(basis<T>::value_[0], basis<T>::value_[1]); }

    virtual std::string basisString() const override
    { return "(abs, arg)"; }
};

/// jacobian from polar to cartesian
template <typename T>
SquareMatrix<T, 2> jac_p_c(const polar<T>& pol)
{
    std::complex<T> v(pol);
    return SquareMatrix<T, 2>({cos(std::arg(v)), -std::abs(v)*sin(std::arg(v)),
                               sin(std::arg(v)),  std::abs(v)*cos(std::arg(v)) });
}

/// jacobian for transformation from cartesian to polar
template <typename T>
SquareMatrix<T, 2> jac_c_p(const cartesian<T>& cart)
{
    std::complex<T> v(cart);
    return SquareMatrix<T, 2>({ std::real(v)/std::abs(v),         std::imag(v)/std::abs(v),
                               -std::imag(v)/pow(std::abs(v), 2), std::real(v)/pow(std::abs(v), 2) });
}

/// convert to string
template <typename T>
inline std::string to_string(const basis<T>& b)
{
    std::string s = b.basisString() + " = ";
    s += "(" + std::to_string(b.value()[0]) + ", " + std::to_string(b.value()[1]) + "); ";
    auto c = std::complex<double>(b);
    s += "as complex: (" + std::to_string(real(c)) + ", " + std::to_string(imag(c)) + ")\n";
    s += "covariance = " + to_string(b.covariance());
    return s;
}

}

}

#endif

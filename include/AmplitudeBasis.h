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

#ifndef yap_AmplitudeBasis_h
#define yap_AmplitudeBasis_h

#include "fwd/AmplitudeBasis.h"

#include "Exceptions.h"
#include "Matrix.h"
#include "Vector.h"

#include <complex>

namespace yap {

namespace basis {

/// \class basis
/// base class for amplitude bases
template <typename T>
class basis {

public:
    /// constructor
    /// \param c1 1st amplitude
    /// \param c2 2nd amplitude
    /// \param c3 3rd amplitude
    /// \param cov 3x3 covariance matrix
    explicit basis(std::complex<T> c1, std::complex<T> c2, std::complex<T> c3, SquareMatrix<T, 3> cov = zeroMatrix<T, 3, 3>()) :
        coordinates_({c1, c2, c3}), covariance_(cov)
    {}

    /// constructor
    /// \param c vector of amplitudes
    /// \param cov 3x3 covariance matrix
    explicit basis(Vector<std::complex<T>, 3> c, SquareMatrix<T, 3> cov = SquareMatrix<T, 3>::zeroMatrix()) :
        coordinates_(c), covariance_(cov)
    {}

    /// constructor
    /// \param cov 3x3 covariance matrix
    explicit basis(SquareMatrix<T, 3> cov) :
        coordinates_({0, 0, 0}), covariance_(cov)
    {}

    /// \return coordinates
    const Vector<std::complex<T>, 3>& coordinates() const
    { return coordinates_; }

    /// \return covariance matrix
    const SquareMatrix<T, 3>& covariance() const
    { return covariance_; }

protected:
    // conversion constructor
    explicit basis(const basis<T>& other, const SquareMatrix<T, 3>& jacobian) :
        coordinates_(jacobian * other.coordinates_),
        covariance_(jacobian * other.covariance_ * jacobian.transpose())
    {}

private:
    /// vector of coordinates
    Vector<std::complex<T>, 3> coordinates_;
    /// covariance matrix
    SquareMatrix<T, 3> covariance_;
};

/// \struct canonical
/// stores amplitudes in canonical basis (S, P, D)
/// \tparam T type stored in amplitudes
/// \defgroup AmplitudeBasis struct for converting among amplitude bases
template <typename T>
class canonical : public basis<T> {

public:
    /// forward all base class constructors
    // order of the amplitudes is: S, P, D
    template<typename... Args>
    canonical(Args&&... args) :
        basis<T>(std::forward<Args>(args)...)
    {}

    /// casting constructor
    explicit constexpr canonical(const transversity<T>& t) :
        basis<T>(t, t_to_c)
    {}

    /// casting constructor
    explicit constexpr canonical(const helicity<T>& h) :
        basis<T>(h, h_to_c)
    {}

    /// \return S amplitude (l=0)
    const std::complex<T> s() const
    { return basis<T>::coordinates()[0]; }

    /// \return P amplitude (l=1)
    const std::complex<T> p() const
    { return basis<T>::coordinates()[1]; }

    /// \return D amplitude (l=2)
    const std::complex<T> d() const
    { return basis<T>::coordinates()[2]; }

    /// access canonical amplitude via angular momentum
    /// \param l angular momentum to retrieve amplitude for
    const std::complex<T> operator[] (size_t l)
    {
        switch (l) {
            case 0:
                return s();
            case 1:
                return p();
            case 2:
                return d();
            default :
                throw exceptions::Exception("l must be 0, 1, or 2", "canonical::operator[]");
        }
    }

private:
    // transformation matrix (jacobian) from transversity to canonical
    static const SquareMatrix<T, 3> t_to_c;

    // transformation matrix (jacobian) from helicity to canonical
    static const SquareMatrix<T, 3> h_to_c;
};

/// \class transversity
/// Stores amplitudes in transversity basis (longitudinal, parallel, perpendicular)
/// \ingroup AmplitudeBasis
template <typename T>
class transversity : public basis<T> {

public:
    /// forward all base class constructors
    // order of the amplitudes is: longitudinal, parallel, perpendicular
    template<typename... Args>
    transversity(Args&&... args) :
        basis<T>(std::forward<Args>(args)...)
    {}

    /// casting constructor
    explicit constexpr transversity(const canonical<T>& c) :
        basis<T>(c, c_to_t)
    {}

    /// casting constructor
    explicit constexpr transversity(const helicity<T>& h) :
        basis<T>(h, h_to_t)
    {}

    /// \return longitudinal amplitude
    const std::complex<T> longitudinal() const
    { return basis<T>::coordinates()[0]; }

    /// \return parallel amplitude
    const std::complex<T> parallel() const
    { return basis<T>::coordinates()[1]; }

    /// \return perpendicular amplitude
    const std::complex<T> perpendicular() const
    { return basis<T>::coordinates()[2]; }

private:
    // transformation matrix (jacobian) from canonical to transversity
    static const SquareMatrix<T, 3> c_to_t;

    // transformation matrix (jacobian) from helicity to transversity
    static const SquareMatrix<T, 3> h_to_t;
};


/// \class helicity
/// Stores amplitudes in helicity basis (zero, plus, minus)
/// \ingroup AmplitudeBasis
template <typename T>
class helicity : public basis<T> {

public:
    /// forward all base class constructors
    // order of the amplitudes is: zero, plus, minus
    template<typename... Args>
    helicity(Args&&... args) :
        basis<T>(std::forward<Args>(args)...)
    {}

    /// casting constructor
    explicit constexpr helicity(const canonical<T>& c) :
        basis<T>(c, c_to_h)
    {}

    /// casting constructor
    explicit constexpr helicity(const transversity<T>& t) :
        basis<T>(t, t_to_h)
    {}

    /// \return zero amplitude A_0
    const std::complex<T> zero() const
    { return basis<T>::coordinates()[0]; }

    /// \return plus amplitude A_+
    const std::complex<T> plus() const
    { return basis<T>::coordinates()[1]; }

    /// \return minus amplitude A_-
    const std::complex<T> minus() const
    { return basis<T>::coordinates()[2]; }

private:
    // transformation matrix (jacobian) from canonical to helicity
    static const SquareMatrix<T, 3> c_to_h;

    // transformation matrix (jacobian) from transversity to helicity
    static const SquareMatrix<T, 3> t_to_h;
};


/// transformation matrix definitions

template <typename T>                        //  longitudinal parallel     perpendicular
const SquareMatrix<T, 3> canonical<T>::t_to_c({{-sqrt(1./3.), sqrt(2./3.), 0},     // S
                                               { 0,           0,           1},     // P
                                               { sqrt(2./3.), sqrt(1./3.), 0}});   // D

template <typename T>                        //  zero         plus          minus
const SquareMatrix<T, 3> canonical<T>::h_to_c({{-sqrt(1./3.), sqrt(1./3.),  sqrt(1./3.)},    // S
                                               { 0,           sqrt(0.5),   -sqrt(0.5)  },    // P
                                               { sqrt(2./3.), sqrt(1./6.),  sqrt(1./6.)}});  // D

template <typename T>                           //  S              P  D
const SquareMatrix<T, 3> transversity<T>::c_to_t({{-sqrt(1. / 3.), 0, sqrt(2. / 3.)},   // longitudinal
                                                  { sqrt(2. / 3.), 0, sqrt(1. / 3.)},   // parallel
                                                  { 0,             1, 0            }}); // perpendicular

template <typename T>                           // zero plus        minus
const SquareMatrix<T, 3> transversity<T>::h_to_t({{1,   0,          0        },    // longitudinal
                                                  {0,   sqrt(0.5),  sqrt(0.5)},    // parallel
                                                  {0,   sqrt(0.5), -sqrt(0.5)}});  // perpendicular

template <typename T>                       //  S             P          D
const SquareMatrix<T, 3> helicity<T>::c_to_h({{-sqrt(1./3.),  0,         sqrt(2./3.)},    // zero
                                              { sqrt(1./3.),  sqrt(0.5), sqrt(1./6.)},    // plus
                                              { sqrt(1./3.), -sqrt(0.5), sqrt(1./6.)}});  // minus

template <typename T>                       // longitudinal parallel    perpendicular
const SquareMatrix<T, 3> helicity<T>::t_to_h({{1,           0,          0        },    // zero
                                              {0,           sqrt(0.5),  sqrt(0.5)},    // plus
                                              {0,           sqrt(0.5), -sqrt(0.5)}});  // minus
}

}

#endif

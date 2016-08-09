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

#include "ComplexBasis.h"

#include "Exceptions.h"
#include "Matrix.h"
#include "Vector.h"

#include <complex>

namespace yap {

namespace amplitude_basis {

/// \typedef covariance_type
template <typename T>
using covariance_type = SquareMatrix<complex_basis::covariance_type<T>, 3>;

/// \class basis
/// base class for spin amplitude bases
/// 3x3 covariance matrix (between coordinates) of 2x2 covariances (between real and imaginary parts of amplitudes)
template <typename T>
class basis {

public:
    /// constructor
    /// \param c1 1st amplitude
    /// \param c2 2nd amplitude
    /// \param c3 3rd amplitude
    /// \param cov 3x3 covariance matrix of 2x2 covariances (between real and imaginary parts of amplitudes), defaults to 0
    explicit basis(const Vector<T, 2>& c1, const Vector<T, 2>& c2, const Vector<T, 2>& c3,
            covariance_type<T> cov = covariance_type<T>()) :
        coordinates_({c1, c2, c3}), covariance_(cov)
    {}

    /// constructor
    /// \param c vector of amplitudes
    /// \param cov 3x3 covariance matrix of 2x2 covariances (between real and imaginary parts of amplitudes), defaults to 0
    explicit basis(Vector<Vector<T, 2>, 3> c,
            covariance_type<T> cov = covariance_type<T>()) :
        coordinates_(c), covariance_(cov)
    {}

    /// constructor
    /// \param cov 3x3 covariance matrix of 2x2 covariances (between real and imaginary parts of amplitudes)
    explicit basis(const covariance_type<T>& cov) :
        coordinates_({0, 0, 0}), covariance_(cov)
    {}

    /// constructor
    /// \param c1 1st amplitude with covariance
    /// \param c2 2nd amplitude with covariance
    /// \param c3 3rd amplitude with covariance
    basis(const complex_basis::basis<T>& c1, const complex_basis::basis<T>& c2, const complex_basis::basis<T>& c3) :
        coordinates_({c1.value(), c2.value(), c3.value()}),
        covariance_(diagonalMatrix<complex_basis::covariance_type<T>, 3>({c1.covariance(), c2.covariance(), c3.covariance()}))
    {}

    /// \return coordinates
    const Vector<Vector<T, 2>, 3>& coordinates() const
    { return coordinates_; }

    /// \return covariance matrix
    const covariance_type<T>& covariance() const
    { return covariance_; }

protected:
    // conversion constructor
    explicit basis(const basis<T>& other, const SquareMatrix<T, 3>& jacobian) :
        coordinates_(jacobian * other.coordinates_),
        covariance_(jacobian * other.covariance_ * transpose(jacobian))
    {}

private:
    /// vector of coordinates
    Vector<Vector<T, 2>, 3> coordinates_;
    /// covariance matrix
    covariance_type<T> covariance_;
};

/// \struct canonical
/// stores amplitudes in canonical basis (S, P, D)
/// \tparam T type stored in amplitudes
/// \defgroup AmplitudeBasis struct for converting among amplitude bases
template <typename T>
class canonical : public basis<T> {

public:
    /// todo get rid of forwarding, or at least make sure casting constructor DOES NOT forward to base copy constructor
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
    const Vector<T, 2> s() const
    { return basis<T>::coordinates()[0]; }

    /// \return P amplitude (l=1)
    const Vector<T, 2> p() const
    { return basis<T>::coordinates()[1]; }

    /// \return D amplitude (l=2)
    const Vector<T, 2> d() const
    { return basis<T>::coordinates()[2]; }

    /// access canonical amplitude via angular momentum
    /// \param l angular momentum to retrieve amplitude for
    const Vector<T, 2> operator[] (size_t l) const
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
    const Vector<T, 2> longitudinal() const
    { return basis<T>::coordinates()[0]; }

    /// \return parallel amplitude
    const Vector<T, 2> parallel() const
    { return basis<T>::coordinates()[1]; }

    /// \return perpendicular amplitude
    const Vector<T, 2> perpendicular() const
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
    const Vector<T, 2> zero() const
    { return basis<T>::coordinates()[0]; }

    /// \return plus amplitude A_+
    const Vector<T, 2> plus() const
    { return basis<T>::coordinates()[1]; }

    /// \return minus amplitude A_-
    const Vector<T, 2> minus() const
    { return basis<T>::coordinates()[2]; }

private:
    // transformation matrix (jacobian) from canonical to helicity
    static const SquareMatrix<T, 3> c_to_h;

    // transformation matrix (jacobian) from transversity to helicity
    static const SquareMatrix<T, 3> t_to_h;
};


/// transformation matrix definitions

template <typename T>                        //  longitudinal parallel     perpendicular
const SquareMatrix<T, 3> canonical<T>::t_to_c({-sqrt(1./3.), sqrt(2./3.), 0,     // S
                                                0,           0,           1,     // P
                                                sqrt(2./3.), sqrt(1./3.), 0});   // D

template <typename T>                        //  zero         plus          minus
const SquareMatrix<T, 3> canonical<T>::h_to_c({-sqrt(1./3.), sqrt(1./3.),  sqrt(1./3.),    // S
                                                0,           sqrt(0.5),   -sqrt(0.5)  ,    // P
                                                sqrt(2./3.), sqrt(1./6.),  sqrt(1./6.)});  // D

template <typename T>                           //  S              P  D
const SquareMatrix<T, 3> transversity<T>::c_to_t({-sqrt(1. / 3.), 0, sqrt(2. / 3.),   // longitudinal
                                                   sqrt(2. / 3.), 0, sqrt(1. / 3.),   // parallel
                                                   0,             1, 0            }); // perpendicular

template <typename T>                           // zero plus        minus
const SquareMatrix<T, 3> transversity<T>::h_to_t({1,   0,          0        ,    // longitudinal
                                                  0,   sqrt(0.5),  sqrt(0.5),    // parallel
                                                  0,   sqrt(0.5), -sqrt(0.5)});  // perpendicular

template <typename T>                       //  S             P          D
const SquareMatrix<T, 3> helicity<T>::c_to_h({-sqrt(1./3.),  0,         sqrt(2./3.),    // zero
                                               sqrt(1./3.),  sqrt(0.5), sqrt(1./6.),    // plus
                                               sqrt(1./3.), -sqrt(0.5), sqrt(1./6.)});  // minus

template <typename T>                       // longitudinal parallel    perpendicular
const SquareMatrix<T, 3> helicity<T>::t_to_h({1,           0,          0        ,    // zero
                                              0,           sqrt(0.5),  sqrt(0.5),    // plus
                                              0,           sqrt(0.5), -sqrt(0.5)});  // minus
}

}

#endif

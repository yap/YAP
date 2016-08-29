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

/// \typedef covariance_matrix
template <typename T>
using covariance_matrix = SquareMatrix<complex_basis::covariance_matrix<T>, 3>;

/// \class basis
/// base class for spin amplitude bases
/// 3x3 covariance matrix (between coordinates) of 2x2 covariances (between real and imaginary parts of amplitudes)
template <typename T>
class basis {

protected:
    /// constructor
    /// \param c vector of amplitudes
    /// \param cov 3x3 covariance matrix of 2x2 covariances (between real and imaginary parts of amplitudes), defaults to 0
    basis(Vector<std::complex<T>, 3> c,
            covariance_matrix<T> cov = covariance_matrix<T>()) :
        Amplitudes_(c), Covariance_(cov)
    {}

    // conversion constructor
    explicit basis(const basis<T>& other, const SquareMatrix<T, 3>& jacobian) :
        Amplitudes_(jacobian * other.Amplitudes_),
        Covariance_(jacobian * other.Covariance_ * transpose(jacobian))
    {}

public:
    /// \return coordinates
    const Vector<std::complex<T>, 3>& coordinates() const
    { return Amplitudes_; }

    /// \return covariance matrix
    const covariance_matrix<T>& covariance() const
    { return Covariance_; }

private:
    /// vector of amplitudes
    Vector<std::complex<T>, 3> Amplitudes_;
    /// covariance matrix
    covariance_matrix<T> Covariance_;
};

/// \class canonical
/// stores amplitudes in canonical basis (S, P, D)
/// \tparam T type stored in amplitudes
/// \defgroup AmplitudeBasis Amplitude bases
template <typename T>
class canonical : public basis<T> {

public:
    /// constructor
    /// \param c1 S amplitude
    /// \param c2 P amplitude
    /// \param c3 D amplitude
    /// \param cov 3x3 covariance matrix of 2x2 covariances (between real and imaginary parts of amplitudes), defaults to 0
    explicit canonical(const std::complex<T>& c1, const std::complex<T>& c2, const std::complex<T>& c3,
            covariance_matrix<T> cov = covariance_matrix<T>()) :
        basis<T>(Vector<std::complex<T>, 3>({c1, c2, c3}), cov)
    {}

    /// constructor
    /// \param c vector of amplitudes
    /// \param cov 3x3 covariance matrix of 2x2 covariances (between real and imaginary parts of amplitudes), defaults to 0
    explicit canonical(Vector<std::complex<T>, 3> c,
            covariance_matrix<T> cov = covariance_matrix<T>()) :
        basis<T>(c, cov)
    {}

    /// constructor
    /// \param cov 3x3 covariance matrix of 2x2 covariances (between real and imaginary parts of amplitudes)
    explicit canonical(const covariance_matrix<T>& cov) :
        basis<T>(Vector<std::complex<T>, 3>({0, 0, 0}), cov)
    {}

    /// constructor
    /// \param c1 S amplitude with covariance
    /// \param c2 P amplitude with covariance
    /// \param c3 D amplitude with covariance
    canonical(const complex_basis::cartesian<T>& c1, const complex_basis::cartesian<T>& c2, const complex_basis::cartesian<T>& c3) :
        basis<T>(Vector<std::complex<T>, 3>({c1, c2, c3}),
                 diagonalMatrix<complex_basis::covariance_matrix<T>, 3>({c1.covariance(), c2.covariance(), c3.covariance()}))
    {}

    /// casting constructor
    explicit constexpr canonical(const transversity<T>& t) :
        basis<T>(t, jacobian_from_transversity)
    {}

    /// casting constructor
    explicit constexpr canonical(const helicity<T>& h) :
        basis<T>(h, jacobian_from_helicity)
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
    const std::complex<T> operator[] (size_t l) const
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
    static const SquareMatrix<T, 3> jacobian_from_transversity;

    // transformation matrix (jacobian) from helicity to canonical
    static const SquareMatrix<T, 3> jacobian_from_helicity;
};

/// \class transversity
/// Stores amplitudes in transversity basis (longitudinal, parallel, perpendicular)
/// \ingroup AmplitudeBasis
template <typename T>
class transversity : public basis<T> {

public:
    /// constructor
    /// \param c1 longitudinal amplitude
    /// \param c2 parallel amplitude
    /// \param c3 perpendicular amplitude
    /// \param cov 3x3 covariance matrix of 2x2 covariances (between real and imaginary parts of amplitudes), defaults to 0
    explicit transversity(const std::complex<T>& c1, const std::complex<T>& c2, const std::complex<T>& c3,
            covariance_matrix<T> cov = covariance_matrix<T>()) :
        basis<T>(Vector<std::complex<T>, 3>({c1, c2, c3}), cov)
    {}

    /// constructor
    /// \param c vector of amplitudes
    /// \param cov 3x3 covariance matrix of 2x2 covariances (between real and imaginary parts of amplitudes), defaults to 0
    explicit transversity(Vector<std::complex<T>, 3> c,
            covariance_matrix<T> cov = covariance_matrix<T>()) :
        basis<T>(c, cov)
    {}

    /// constructor
    /// \param cov 3x3 covariance matrix of 2x2 covariances (between real and imaginary parts of amplitudes)
    explicit transversity(const covariance_matrix<T>& cov) :
        basis<T>(Vector<std::complex<T>, 3>({0, 0, 0}), cov)
    {}

    /// constructor
    /// \param c1 longitudinal amplitude with covariance
    /// \param c2 parallel amplitude with covariance
    /// \param c3 perpendicular amplitude with covariance
    transversity(const complex_basis::cartesian<T>& c1, const complex_basis::cartesian<T>& c2, const complex_basis::cartesian<T>& c3) :
        basis<T>(Vector<std::complex<T>, 3>({c1, c2, c3}),
                 diagonalMatrix<complex_basis::covariance_matrix<T>, 3>({c1.covariance(), c2.covariance(), c3.covariance()}))
    {}

    /// casting constructor
    explicit constexpr transversity(const canonical<T>& c) :
        basis<T>(c, jacobian_from_canonical)
    {}

    /// casting constructor
    explicit constexpr transversity(const helicity<T>& h) :
        basis<T>(h, jacobian_from_helicity)
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
    static const SquareMatrix<T, 3> jacobian_from_canonical;

    // transformation matrix (jacobian) from helicity to transversity
    static const SquareMatrix<T, 3> jacobian_from_helicity;
};


/// \class helicity
/// Stores amplitudes in helicity basis (zero, plus, minus)
/// \ingroup AmplitudeBasis
template <typename T>
class helicity : public basis<T> {

public:
    /// constructor
    /// \param c1 zero amplitude
    /// \param c2 plus amplitude
    /// \param c3 minus amplitude
    /// \param cov 3x3 covariance matrix of 2x2 covariances (between real and imaginary parts of amplitudes), defaults to 0
    explicit helicity(const std::complex<T>& c1, const std::complex<T>& c2, const std::complex<T>& c3,
            covariance_matrix<T> cov = covariance_matrix<T>()) :
        basis<T>(Vector<std::complex<T>, 3>({c1, c2, c3}), cov)
    {}

    /// constructor
    /// \param c vector of amplitudes
    /// \param cov 3x3 covariance matrix of 2x2 covariances (between real and imaginary parts of amplitudes), defaults to 0
    explicit helicity(Vector<std::complex<T>, 3> c,
            covariance_matrix<T> cov = covariance_matrix<T>()) :
        basis<T>(c, cov)
    {}

    /// constructor
    /// \param cov 3x3 covariance matrix of 2x2 covariances (between real and imaginary parts of amplitudes)
    explicit helicity(const covariance_matrix<T>& cov) :
        basis<T>(Vector<std::complex<T>, 3>({0, 0, 0}), cov)
    {}

    /// constructor
    /// \param c1 zero amplitude with covariance
    /// \param c2 plus amplitude with covariance
    /// \param c3 minus amplitude with covariance
    helicity(const complex_basis::cartesian<T>& c1, const complex_basis::cartesian<T>& c2, const complex_basis::cartesian<T>& c3) :
        basis<T>(Vector<std::complex<T>, 3>({c1, c2, c3}),
                 diagonalMatrix<complex_basis::covariance_matrix<T>, 3>({c1.covariance(), c2.covariance(), c3.covariance()}))
    {}

    /// casting constructor
    explicit constexpr helicity(const canonical<T>& c) :
        basis<T>(c, jacobian_from_canonical)
    {}

    /// casting constructor
    explicit constexpr helicity(const transversity<T>& t) :
        basis<T>(t, jacobian_from_transversity)
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
    static const SquareMatrix<T, 3> jacobian_from_canonical;

    // transformation matrix (jacobian) from transversity to helicity
    static const SquareMatrix<T, 3> jacobian_from_transversity;
};


/// transformation matrix definitions

template <typename T>
const SquareMatrix<T, 3> canonical<T>::jacobian_from_transversity(
      //  longitudinal parallel     perpendicular
        {-sqrt(1./3.), sqrt(2./3.), 0,     // S
          0,           0,           1,     // P
          sqrt(2./3.), sqrt(1./3.), 0});   // D

template <typename T>
const SquareMatrix<T, 3> canonical<T>::jacobian_from_helicity(
      //  zero         plus          minus
        {-sqrt(1./3.), sqrt(1./3.),  sqrt(1./3.),    // S
          0,           sqrt(0.5),   -sqrt(0.5)  ,    // P
          sqrt(2./3.), sqrt(1./6.),  sqrt(1./6.)});  // D

template <typename T>
const SquareMatrix<T, 3> transversity<T>::jacobian_from_canonical(
      //  S              P  D
        {-sqrt(1. / 3.), 0, sqrt(2. / 3.),   // longitudinal
          sqrt(2. / 3.), 0, sqrt(1. / 3.),   // parallel
          0,             1, 0            }); // perpendicular

template <typename T>
const SquareMatrix<T, 3> transversity<T>::jacobian_from_helicity(
      // zero plus        minus
        {1,   0,          0        ,    // longitudinal
         0,   sqrt(0.5),  sqrt(0.5),    // parallel
         0,   sqrt(0.5), -sqrt(0.5)});  // perpendicular

template <typename T>
const SquareMatrix<T, 3> helicity<T>::jacobian_from_canonical(
      //  S             P          D
        {-sqrt(1./3.),  0,         sqrt(2./3.),    // zero
          sqrt(1./3.),  sqrt(0.5), sqrt(1./6.),    // plus
          sqrt(1./3.), -sqrt(0.5), sqrt(1./6.)});  // minus

template <typename T>
const SquareMatrix<T, 3> helicity<T>::jacobian_from_transversity(
      // longitudinal parallel    perpendicular
        {1,           0,          0        ,    // zero
         0,           sqrt(0.5),  sqrt(0.5),    // plus
         0,           sqrt(0.5), -sqrt(0.5)});  // minus
}

}

#endif

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
    basis(Vector<Vector<T, 2>, 3> c,
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
    const Vector<Vector<T, 2>, 3>& coordinates() const
    { return Amplitudes_; }

    /// \return covariance matrix
    const covariance_matrix<T>& covariance() const
    { return Covariance_; }

private:
    /// vector of amplitudes
    Vector<Vector<T, 2>, 3> Amplitudes_;
    /// covariance matrix
    covariance_matrix<T> Covariance_;
};

/// \class canonical
/// stores amplitudes in canonical basis (S, P, D)
/// \tparam T type stored in amplitudes
/// \defgroup AmplitudeBasis struct for converting among amplitude bases
template <typename T>
class canonical : public basis<T> {

public:
    /// constructor
    /// \param c1 S amplitude
    /// \param c2 P amplitude
    /// \param c3 D amplitude
    /// \param cov 3x3 covariance matrix of 2x2 covariances (between real and imaginary parts of amplitudes), defaults to 0
    explicit canonical(const Vector<T, 2>& c1, const Vector<T, 2>& c2, const Vector<T, 2>& c3,
            covariance_matrix<T> cov = covariance_matrix<T>()) :
        basis<T>(Vector<Vector<T, 2>, 3>({c1, c2, c3}), cov)
    {}

    /// constructor
    /// \param c vector of amplitudes
    /// \param cov 3x3 covariance matrix of 2x2 covariances (between real and imaginary parts of amplitudes), defaults to 0
    explicit canonical(Vector<Vector<T, 2>, 3> c,
            covariance_matrix<T> cov = covariance_matrix<T>()) :
        basis<T>(c, cov)
    {}

    /// constructor
    /// \param cov 3x3 covariance matrix of 2x2 covariances (between real and imaginary parts of amplitudes)
    explicit canonical(const covariance_matrix<T>& cov) :
        basis<T>(Vector<Vector<T, 2>, 3>({0, 0, 0}), cov)
    {}

    /// constructor
    /// \param c1 S amplitude with covariance
    /// \param c2 P amplitude with covariance
    /// \param c3 D amplitude with covariance
    canonical(const complex_basis::basis<T>& c1, const complex_basis::basis<T>& c2, const complex_basis::basis<T>& c3) :
        basis<T>(Vector<Vector<T, 2>, 3>({c1.value(), c2.value(), c3.value()}),
                 diagonalMatrix<complex_basis::covariance_matrix<T>, 3>({c1.covariance(), c2.covariance(), c3.covariance()}))
    {}

    /// casting constructor
    explicit constexpr canonical(const transversity<T>& t) :
        basis<T>(t, jacobian_t_to_c)
    {}

    /// casting constructor
    explicit constexpr canonical(const helicity<T>& h) :
        basis<T>(h, jacobian_h_to_c)
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
    static const SquareMatrix<T, 3> jacobian_t_to_c;

    // transformation matrix (jacobian) from helicity to canonical
    static const SquareMatrix<T, 3> jacobian_h_to_c;
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
    explicit transversity(const Vector<T, 2>& c1, const Vector<T, 2>& c2, const Vector<T, 2>& c3,
            covariance_matrix<T> cov = covariance_matrix<T>()) :
        basis<T>(Vector<Vector<T, 2>, 3>({c1, c2, c3}), cov)
    {}

    /// constructor
    /// \param c vector of amplitudes
    /// \param cov 3x3 covariance matrix of 2x2 covariances (between real and imaginary parts of amplitudes), defaults to 0
    explicit transversity(Vector<Vector<T, 2>, 3> c,
            covariance_matrix<T> cov = covariance_matrix<T>()) :
        basis<T>(c, cov)
    {}

    /// constructor
    /// \param cov 3x3 covariance matrix of 2x2 covariances (between real and imaginary parts of amplitudes)
    explicit transversity(const covariance_matrix<T>& cov) :
        basis<T>(Vector<Vector<T, 2>, 3>({0, 0, 0}), cov)
    {}

    /// constructor
    /// \param c1 longitudinal amplitude with covariance
    /// \param c2 parallel amplitude with covariance
    /// \param c3 perpendicular amplitude with covariance
    transversity(const complex_basis::basis<T>& c1, const complex_basis::basis<T>& c2, const complex_basis::basis<T>& c3) :
        basis<T>(Vector<Vector<T, 2>, 3>({c1.value(), c2.value(), c3.value()}),
                 diagonalMatrix<complex_basis::covariance_matrix<T>, 3>({c1.covariance(), c2.covariance(), c3.covariance()}))
    {}

    /// casting constructor
    explicit constexpr transversity(const canonical<T>& c) :
        basis<T>(c, jacobian_c_to_t)
    {}

    /// casting constructor
    explicit constexpr transversity(const helicity<T>& h) :
        basis<T>(h, jacobian_h_to_t)
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
    static const SquareMatrix<T, 3> jacobian_c_to_t;

    // transformation matrix (jacobian) from helicity to transversity
    static const SquareMatrix<T, 3> jacobian_h_to_t;
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
    explicit helicity(const Vector<T, 2>& c1, const Vector<T, 2>& c2, const Vector<T, 2>& c3,
            covariance_matrix<T> cov = covariance_matrix<T>()) :
        basis<T>(Vector<Vector<T, 2>, 3>({c1, c2, c3}), cov)
    {}

    /// constructor
    /// \param c vector of amplitudes
    /// \param cov 3x3 covariance matrix of 2x2 covariances (between real and imaginary parts of amplitudes), defaults to 0
    explicit helicity(Vector<Vector<T, 2>, 3> c,
            covariance_matrix<T> cov = covariance_matrix<T>()) :
        basis<T>(c, cov)
    {}

    /// constructor
    /// \param cov 3x3 covariance matrix of 2x2 covariances (between real and imaginary parts of amplitudes)
    explicit helicity(const covariance_matrix<T>& cov) :
        basis<T>(Vector<Vector<T, 2>, 3>({0, 0, 0}), cov)
    {}

    /// constructor
    /// \param c1 zero amplitude with covariance
    /// \param c2 plus amplitude with covariance
    /// \param c3 minus amplitude with covariance
    helicity(const complex_basis::basis<T>& c1, const complex_basis::basis<T>& c2, const complex_basis::basis<T>& c3) :
        basis<T>(Vector<Vector<T, 2>, 3>({c1.value(), c2.value(), c3.value()}),
                 diagonalMatrix<complex_basis::covariance_matrix<T>, 3>({c1.covariance(), c2.covariance(), c3.covariance()}))
    {}

    /// casting constructor
    explicit constexpr helicity(const canonical<T>& c) :
        basis<T>(c, jacobian_c_to_h)
    {}

    /// casting constructor
    explicit constexpr helicity(const transversity<T>& t) :
        basis<T>(t, jacobian_t_to_h)
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
    static const SquareMatrix<T, 3> jacobian_c_to_h;

    // transformation matrix (jacobian) from transversity to helicity
    static const SquareMatrix<T, 3> jacobian_t_to_h;
};


/// transformation matrix definitions

template <typename T>
const SquareMatrix<T, 3> canonical<T>::jacobian_t_to_c(
      //  longitudinal parallel     perpendicular
        {-sqrt(1./3.), sqrt(2./3.), 0,     // S
          0,           0,           1,     // P
          sqrt(2./3.), sqrt(1./3.), 0});   // D

template <typename T>
const SquareMatrix<T, 3> canonical<T>::jacobian_h_to_c(
      //  zero         plus          minus
        {-sqrt(1./3.), sqrt(1./3.),  sqrt(1./3.),    // S
          0,           sqrt(0.5),   -sqrt(0.5)  ,    // P
          sqrt(2./3.), sqrt(1./6.),  sqrt(1./6.)});  // D

template <typename T>
const SquareMatrix<T, 3> transversity<T>::jacobian_c_to_t(
      //  S              P  D
        {-sqrt(1. / 3.), 0, sqrt(2. / 3.),   // longitudinal
          sqrt(2. / 3.), 0, sqrt(1. / 3.),   // parallel
          0,             1, 0            }); // perpendicular

template <typename T>
const SquareMatrix<T, 3> transversity<T>::jacobian_h_to_t(
      // zero plus        minus
        {1,   0,          0        ,    // longitudinal
         0,   sqrt(0.5),  sqrt(0.5),    // parallel
         0,   sqrt(0.5), -sqrt(0.5)});  // perpendicular

template <typename T>
const SquareMatrix<T, 3> helicity<T>::jacobian_c_to_h(
      //  S             P          D
        {-sqrt(1./3.),  0,         sqrt(2./3.),    // zero
          sqrt(1./3.),  sqrt(0.5), sqrt(1./6.),    // plus
          sqrt(1./3.), -sqrt(0.5), sqrt(1./6.)});  // minus

template <typename T>
const SquareMatrix<T, 3> helicity<T>::jacobian_t_to_h(
      // longitudinal parallel    perpendicular
        {1,           0,          0        ,    // zero
         0,           sqrt(0.5),  sqrt(0.5),    // plus
         0,           sqrt(0.5), -sqrt(0.5)});  // minus
}

}

#endif

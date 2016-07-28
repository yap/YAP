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

#include <complex>

namespace yap {

namespace basis {

/// \struct canonical
/// stores amplitudes in canonical basis
/// \tparam T type stored in amplitudes
/// \defgroup AmplitudeBasis struct for converting among amplitude bases
template <typename T>
struct canonical {
    /// S amplitude (l=0)
    std::complex<T> S;
    /// P amplitude (l=1)
    std::complex<T> P;
    /// D amplitude (l=2)
    std::complex<T> D;

    /// constructor
    /// \param s S amplitude
    /// \param p P amplitude
    /// \param d D amplitude
    explicit constexpr canonical(std::complex<T> s, std::complex<T> p, std::complex<T> d) :
        S(s), P(p), D(d) {};

    /// casting constructor
    explicit constexpr canonical(const transversity<T>& t) :
        S(sqrt(2. / 3.) * t.parallel - sqrt(1. / 3.) * t.longitudinal),
        P(t.perpendicular),
        D(sqrt(1. / 3.) * t.parallel + sqrt(2. / 3.) * t.longitudinal)
    {}

    /// casting constructor
    explicit constexpr canonical(const helicity<T>& h) :
        S(sqrt(1. / 3.) * (h.plus + h.minus - h.zero)),
        P(sqrt(1. / 2.) * (h.plus - h.minus)),
        D(sqrt(1. / 6.) * (h.plus + h.minus) + sqrt(2. / 3.) * h.zero)
    {}

    /// access canonical amplitude via angular momentum
    /// \param L angular momentum to retrieve amplitude for
    const std::complex<T>& operator[] (size_t L)
    {
        switch (L) {
            case 0:
                return S;
            case 1:
                return P;
            case 2:
                return D;
            default :
                throw exceptions::Exception("L must be 0, 1, or 2", "canonical::operator[]");
        }
    }
};

/// \struct transversity
/// Stores amplitudes in transversity basis
/// \ingroup AmplitudeBasis
template <typename T>
struct transversity {
    /// longitudinal amplitude
    std::complex<T> longitudinal;
    /// parallel amplitude
    std::complex<T> parallel;
    /// perpendicular amplitude
    std::complex<T> perpendicular;

    /// constructor
    /// \param l    longitudinal amplitude
    /// \param par  parallel amplitude
    /// \param perp perpendicular amplitude
    explicit constexpr transversity(std::complex<T> l, std::complex<T> par, std::complex<T> perp) :
        longitudinal(l), parallel(par), perpendicular(perp) {};

    /// casting constructor
    explicit constexpr transversity(const canonical<T>& c) :
        longitudinal (-sqrt(1. / 3.) * c.S + sqrt(2. / 3.) * c.D),
        parallel     ( sqrt(2. / 3.) * c.S + sqrt(1. / 3.) * c.D),
        perpendicular( c.P)
    {}

    /// casting constructor
    explicit constexpr transversity(const helicity<T>& h) :
        longitudinal (h.zero),
        parallel     (sqrt(0.5) * (h.plus + h.minus)),
        perpendicular(sqrt(0.5) * (h.plus - h.minus))
    {}
};

/// \struct helicity
/// Stores amplitudes in helicity basis
/// \ingroup AmplitudeBasis
template <typename T>
struct helicity {
    /// zero amplitude A_0
    std::complex<T> zero;
    /// plus amplitude A_+
    std::complex<T> plus;
    /// minus amplitude A_-
    std::complex<T> minus;

    /// constructor
    /// \param z zero amplitude A_0
    /// \param p plus amplitude A_+
    /// \param m minus amplitude A_-
    explicit constexpr helicity(std::complex<T> z, std::complex<T> p, std::complex<T> m) :
        zero(z), plus(p), minus(m) {};

    /// casting constructor
    explicit constexpr helicity(const canonical<T>& c) :
        zero (-sqrt(1. / 3.) * c.S + sqrt(2. / 3.) * c.D),
        plus ( sqrt(1. / 3.) * c.S + sqrt(1. / 2.) * c.P + sqrt(1. / 6.) * c.D),
        minus( sqrt(1. / 3.) * c.S - sqrt(1. / 2.) * c.P + sqrt(1. / 6.) * c.D)
    {}

    /// casting constructor
    explicit constexpr helicity(const transversity<T>& t) :
        zero (t.longitudinal),
        plus (sqrt(0.5) * (t.parallel + t.perpendicular)),
        minus(sqrt(0.5) * (t.parallel - t.perpendicular))
    {}
};

}

}

#endif

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

#ifndef yap_ClebschGordan_h
#define yap_ClebschGordan_h

#include "MathUtilities.h"

namespace yap {

/// \namespace ClebschGordan
/// \brief Clebsch Gordan coefficients and related spin-checking functions
/// \author Johannes Rauch, Daniel Greenwald
namespace ClebschGordan {

/// \return Clebsch-Gordan coefficient (j1 m1 j2 m2 | J M)
/// Implemented from Eq. (16) from G. Racah, "Theory of Complex Spectra. II", Phys. Rev. 62, 438 (1942)
/// \param two_j1 2*spin of first particle
/// \param two_m1 2*spin-projection of first particle
/// \param two_j2 2*spin of second particle
/// \param two_m2 2*spin-projection of second particle
/// \param two_J  2*spin of composite system
/// \param two_M  2*spin-projection of composite system
double coefficient(unsigned two_j1, int two_m1, unsigned two_j2, int two_m2, unsigned two_J, int two_M);

/// \return Clebsch-Gordan coefficient (j1 m1 j2 m2 | J M), with M = m1 + m2
/// \param two_j1 2*spin of first particle
/// \param two_m1 2*spin-projection of first particle
/// \param two_j2 2*spin of second particle
/// \param two_m2 2*spin-projection of second particle
/// \param two_J  2*spin of composite system
double coefficient(unsigned two_j1, int two_m1, unsigned two_j2, int two_m2, unsigned two_J)
{ return coefficient(two_j1, two_m1, two_j2, two_m2, two_J, two_m1 + two_m2); }

/// \return consistency of spin and spin projection
/// \param two_J 2*spin
/// \param two_M 2*spin-projection
constexpr bool consistent(unsigned two_J, int two_M)
{ return (std::abs(two_M) <= twoJ) and isEven(two_J - two_M); }

/// \returns whether j1 and j2 can couple to J
/// \param two_j1 2*spin of first particle
/// \param two_j2 2*spin of second particle
/// \param two_J  2*spin of composite system
constexpr bool couple(unsigned two_j1, unsigned two_j2, unsigned two_J)
{ return (two_J <= (two_j1 + two_j2)) and (two_J >= std::abs(two_j1 - two_j2)) and is_even(two_J - two_j1 - two_j2); }
        
} // ClebschGordon namespace

} // yap namespace

#endif

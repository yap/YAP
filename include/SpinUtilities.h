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

#ifndef yap_SpinUtilities_h
#define yap_SpinUtilities_h

#include <string>

namespace yap {

/// \brief Clebsch Gordan coefficients and related functions
/// \author Johannes Rauch, Daniel Greenwald
/// This code has been taken from rootpwa and modified

/// \return Clebsch-Gordan coefficient (j1 m1 j2 m2 | J M)
/// \param two_j1 2*spin of first particle
/// \param two_m1 2*spin-projection of first particle
/// \param two_j2 2*spin of second particle
/// \param two_m2 2*spin-projection of second particle
/// \param two_J  2*spin of composite system
/// \param two_M  2*spin-projection of composite system
double clebschGordan(int two_j1, int two_m1, int two_j2, int two_m2, int two_J, int two_M);

/// \return consistency of spin and spin projection
/// \param two_J 2*spin
/// \param two_M 2*spin-projection
bool spinAndProjAreCompatible(int two_J, int two_M);

/// \returns whether j1 and j2 can couple to J
/// \param two_j1 2*spin of first particle
/// \param two_j2 2*spin of second particle
/// \param two_J  2*spin of composite system
bool spinStatesCanCouple(int two_j1, int two_j2, int two_J);

/// convert 2*J to string (e.g. 1/2, 1, 3/2, etc.)
std::string spinToString(int twoJ);

}

#endif

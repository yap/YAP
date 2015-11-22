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

#ifndef yap_Constants_h
#define yap_Constants_h

#include "SquareMatrix.h"
#include "ThreeVector.h"

#include <complex>

namespace yap {

/* template <typename T> */
/*     typedef FourVector; */

/// complex zero
extern const std::complex<double> Complex_0;

/// complex one
extern const std::complex<double> Complex_1;

/// complex i
extern const std::complex<double> Complex_i;

/// pi (11 digits)
extern const double PI;

/// convert deg to rad by multiplying by; rad to deg by dividing by
extern const double DEGTORAD;

/// X axis
extern const ThreeVector<double> Axis_X;

/// Y axis
extern const ThreeVector<double> Axis_Y;

/// Z axis
extern const ThreeVector<double> Axis_Z;

/// 0 as ThreeVector;
extern const ThreeVector<double> Vect3_0;

/* /// 0 as FourVector; */
/* extern const FourVector<double> Vect4_0; */

/// 3x3 Unit Matrix
extern const SquareMatrix<double, 3> Unit3x3;

}

#endif

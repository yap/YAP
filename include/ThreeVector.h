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

#ifndef yap_ThreeVector_h
#define yap_ThreeVector_h

#include "NVector.h"

namespace yap {

/// \typedef ThreeVector
/// \ingroup VectorAlgebra
template <typename T>
using ThreeVector = NVector<T, 3>;

/// cross product of #ThreeVector's
template <typename T>
ThreeVector<T> cross(const ThreeVector<T>& A, const ThreeVector<T>& B)
{ return {A[2]* B[3] - A[3]* B[2], A[3]* B[1] - A[1]* B[3], A[1]* B[2] - A[2]* B[1]}; }

}
#endif

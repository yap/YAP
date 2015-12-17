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

#ifndef yap_Rotation_h
#define yap_Rotation_h

#include "Constants.h"
#include "Matrix.h"
#include "ThreeVector.h"

namespace yap {

/// \return a 3D rotation matrix
/// \param V vector to rotate about
/// \param theta angle to rotate through
template <typename T>
ThreeMatrix<T> rotation(const ThreeVector<T>& V, const T& theta)
{
    T v = abs(V);

    if (v == 0 || theta == 0)
        return unitMatrix<T, 3>();
    else {
        // normalize direction vector
        ThreeVector<T> U = V * (T(1) / v);

        T cosTheta = cos(theta);

        return (1 - cosTheta) * outer(U, U) + cosTheta * unitMatrix<T, 3>() + (T)sin(theta) * skewSymmetric(U);
    }
}

}
#endif

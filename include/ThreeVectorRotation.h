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

#ifndef yap_ThreeVectorRotation_h
#define yap_ThreeVectorRotation_h

#include "Constants.h"
#include "SquareMatrix.h"
#include "ThreeVector.h"

namespace yap {

/// \class ThreeVectorRotation
/// \brief For rotating three-vectors
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup VectorAlgebra
template <typename T>
class ThreeVectorRotation : public SquareMatrix<T, 3>
{
public:

    /// constructor
    /// \param V ThreeVector defining axis of rotation, with magnitude = cos(theta)
    ThreeVectorRotation(const ThreeVector<T> V, T theta)
    {
        double v = sqrt(V * V);

        if (v == 0 || theta == 0)
            *this = Unit3x3;

        else {

            // normalize direction vector
            ThreeVector<T> U = V * (T(1) / v);

            double cosTheta = cos(theta);

            /* *this = (1 - cosTheta) * outer(V, V) + cosTheta * Unit3x3 + sin(theta) * skewSymmetric(V); */
        }
    }

    /// constructor
    /// \param V ThreeVector defining axis of rotation, with magnitude = theta
    ThreeVectorRotation(const ThreeVector<T> V)
        : ThreeVectorRotation(V, sqrt(V * V)) {}
};

}
#endif

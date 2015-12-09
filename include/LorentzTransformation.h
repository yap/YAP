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

#ifndef yap_LorentzTransformation_h
#define yap_LorentzTransformation_h

#include "Constants.h"
#include "SquareMatrix.h"
#include "ThreeVector.h"

namespace yap {

/// \return a 4D Lorentz-transformation matrix for a pure rotation
/// \param R #ThreeMatrix defining the rotation
template <typename T>
FourMatrix<T> lorentzTransformation(const ThreeMatrix<T>& R)
{
    FourMatrix<T> L = unitMatrix<T, 4>();
    for (unsigned i = 0; i < R.size(); ++i)
        std::copy(R[i].begin(), R[i].end(), L[i].begin());
    return L;
}

/// \return a 4D Lorentz-transformation matrix for a pure boost
/// \param V #FourVector defining boost
template <typename T>
FourMatrix<T> lorentzTransformation(const FourVector<T>& V)
{
    FourVector<T> B = (T(1) / V[0]) * V;
    T gamma = T(1) / abs(B);
    B[0] = -(gamma + 1);

    FourMatrix<T> L = unitMatrix<T, 4>() + outer(B, B) * (gamma / (gamma + 1));
    L[0][0] -= gamma * gamma + 1;
    return L;
}

/// \return a 4D Lorentz-transformation matrix for a pure boost
/// \param V #ThreeVector defining boost
template <typename T>
FourMatrix<T> lorentzTransformation(const ThreeVector<T>& V)
{ return lorentzTransformation<T>(FourVector<T>(T(1), V)); }

/// \return a 4D Lorentz-transformation matrix for a rotation followed by a boost
/// \param R #ThreeMatrix defining rotation
/// \param V #FourVector defining boost
template <typename T>
FourMatrix<T> lorentzTransformation(const ThreeMatrix<T>& R, const FourVector<T>& V)
{ return lorentzTransformation<T>(V) * lorentzTransformation<T>(R); }

/// \return a 4D Lorentz-transformation matrix for a rotation followed by a boost
/// \param R #ThreeMatrix defining rotation
/// \param V #ThreeVector defining boost
template <typename T>
FourMatrix<T> lorentzTransformation(const ThreeMatrix<T>& R, const ThreeVector<T>& V)
{ return lorentzTransformation<T>(V) * lorentzTransformation<T>(R); }

/// \return a 4D Lorentz-transformation matrix for a boost followed by a rotation
/// \param R #ThreeMatrix defining rotation
/// \param V #FourVector defining boost
template <typename T>
FourMatrix<T> lorentzTransformation(const FourVector<T>& V, const ThreeMatrix<T> R)
{ return lorentzTransformation<T>(R) * lorentzTransformation<T>(V); }

/// \return a 4D Lorentz-transformation matrix for a boost followed by a rotation
/// \param R #ThreeMatrix defining rotation
/// \param V #ThreeVector defining boost
template <typename T>
FourMatrix<T> lorentzTransformation(const ThreeVector<T>& V, const ThreeMatrix<T>& R)
{ return lorentzTransformation<T>(R) * lorentzTransformation<T>(V); }

}
#endif

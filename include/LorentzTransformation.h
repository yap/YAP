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

#include "Matrix.h"
#include "FourVector.h"
#include "ThreeVector.h"

#include <algorithm>

namespace yap {

/// \return a 4D Lorentz-transformation matrix for a pure rotation
/// \param R #ThreeMatrix defining the rotation
template <typename T>
FourMatrix<T> lorentzTransformation(const ThreeMatrix<T>& R)
{
    FourMatrix<T> L = unitMatrix<T, 4>();
    for (unsigned i = 0; i < R.size(); ++i)
        std::copy(R[i].begin(), R[i].end(), L[i + 1].begin() + 1);
    return L;
}

/// \return a 4D Lorentz-transformation matrix for a pure boost
/// \param V #FourVector of four-momentum defining boost
template <typename T>
FourMatrix<T> lorentzTransformation(const FourVector<T>& V)
{
    auto gamma = V[0] / abs(V); // E / m
    auto b = FourVector<T>((gamma + 1), gamma * vect(V) / V[0]) / sqrt(gamma + 1); // vect(V)/V[0] = beta

    return diagonalMatrix<T, 4>({ -1, 1, 1, 1}) + outer(b, b);
}

/// \return a 4D Lorentz-transformation matrix for a pure boost
/// \param V #ThreeVector defining boost
template <typename T>
constexpr FourMatrix<T> lorentzTransformation(const ThreeVector<T>& V)
{ return lorentzTransformation<T>(FourVector<T>(T(1), V)); }

/// \return a 4D Lorentz-transformation matrix for a pure boost
/// \param fourVecs the sum of these define the boost
template <typename T>
constexpr FourMatrix<T> lorentzTransformation(const std::vector<FourVector<T> >& fourVecs)
{ return lorentzTransformation(std::accumulate(fourVecs.begin(), fourVecs.end(), FourVector<T>({0, 0, 0, 0}))); }

/// \return a 4D Lorentz-transformation matrix for a rotation followed by a boost
/// \param R #ThreeMatrix defining rotation
/// \param V #FourVector defining boost
template <typename T>
constexpr FourMatrix<T> lorentzTransformation(const ThreeMatrix<T>& R, const FourVector<T>& V)
{ return lorentzTransformation<T>(V) * lorentzTransformation<T>(R); }

/// \return a 4D Lorentz-transformation matrix for a rotation followed by a boost
/// \param R #ThreeMatrix defining rotation
/// \param V #ThreeVector defining boost
template <typename T>
constexpr FourMatrix<T> lorentzTransformation(const ThreeMatrix<T>& R, const ThreeVector<T>& V)
{ return lorentzTransformation<T>(V) * lorentzTransformation<T>(R); }

/// \return a 4D Lorentz-transformation matrix for a boost followed by a rotation
/// \param R #ThreeMatrix defining rotation
/// \param V #FourVector defining boost
template <typename T>
constexpr FourMatrix<T> lorentzTransformation(const FourVector<T>& V, const ThreeMatrix<T> R)
{ return lorentzTransformation<T>(R) * lorentzTransformation<T>(V); }

/// \return a 4D Lorentz-transformation matrix for a boost followed by a rotation
/// \param R #ThreeMatrix defining rotation
/// \param V #ThreeVector defining boost
template <typename T>
constexpr FourMatrix<T> lorentzTransformation(const ThreeVector<T>& V, const ThreeMatrix<T>& R)
{ return lorentzTransformation<T>(R) * lorentzTransformation<T>(V); }

}
#endif

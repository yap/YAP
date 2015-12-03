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

/// \files

#ifndef yap_FourVector_h
#define yap_FourVector_h

#include "CoordinateSystem.h"
#include "SquareMatrix.h"
#include "ThreeVector.h"

#include <algorithm>
#include <type_traits>

namespace yap {

/// \typedef FourVector
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup VectorAlgebra
template <typename T>
using FourVector = typename std::enable_if<std::is_arithmetic<T>::value, NVector<T, 4> >::type;

/// \return #FourVector
/// \param E 0th component
/// \param P #ThreeVector component
template <typename T>
FourVector<T> fourVector(const T& E, const ThreeVector<T>& P)
{ return FourVector<T> {E, P[0], P[1], P[2]}; }

/// \return inner (dot) product for 4-vectors
template <typename T>
T operator*(const FourVector<T>& A, const FourVector<T>& B)
{ return A.front() * B.front() - std::inner_product(A.begin() + 1, A.end(), B.begin() + 1, 0); }

/// unary minus for 4-vector
template <typename T>
FourVector<T> operator-(const FourVector<T>& A)
{ FourVector<T> res = A; std::transform(res.begin() + 1, res.end(), res.begin() + 1, [](const T & t) {return -t;}); return res; }

/// \return Spatial #ThreeVector inside #FourVector
template <typename T>
ThreeVector<T> vect(const FourVector<T>& V)
{ return ThreeVector<T> {V[1], V[2], V[3]}; }

/// \return boost vector of this #FourVector
template <typename T>
ThreeVector<T> boost(const FourVector<T>& V)
{ return (V[0] != 0) ? (T(1) / V[0]) * vect(V) : ThreeVector<T> {0, 0, 0}; }

/// apply a three-rotation to a FourVector (rotating only the spatial components)
template <typename T>
FourVector<T>& operator*(const ThreeMatrix<T>& R, const FourVector<T>& V)
{ return fourVector<T>(V[0], R * vect(V)); }

/// Calculate helicity frame of V transformed from C,
/// with z = unit(V), y = C.z X z, x = y X z
/// \param V Fourector defining new Z direction
/// \param C CoordinateSystem aiding in defining new Y direction
template <typename T>
typename std::enable_if<std::is_arithmetic<T>::value, CoordinateSystem<T, 3> >::type
helicityFrame(const FourVector<T>& V, const CoordinateSystem<T, 3>& C)
{ return helicityFrame(vect<T>(V), C); }

}
#endif

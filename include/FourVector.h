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

#include "ThreeVector.h"
#include "ThreeVectorRotation.h"

#include <algorithm>
#include <array>

namespace yap {


/// \class FourVector
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup VectorAlgebra
template <typename T>
class FourVector : public NVector<T, 4>
{
public:

    /// Constructor
    FourVector(std::initializer_list<T> list) : NVector<T, 4>(list) {}

    /// constructor
    FourVector(const T& E = 0, const NVector<T, 3> P = Vect3_0) : FourVector( {E, P[0], P[1], P[2]}) {}

    /// inner (dot) product of FourVectors
    T operator*(const NVector<T, 4>& B) const override
    { return (*this)[0] * B[0] - std::inner_product(this->begin() + 1, this->end(), B.begin() + 1, 0); }

};


/// \return Spatial #ThreeVector inside #FourVector
template <typename T>
ThreeVector<T> vect(const FourVector<T>& V)
{ return {V[1], V[2], V[3]}; }

/// \return boost vector of this #FourVector
template <typename T>
ThreeVector<T> boost(const FourVector<T>& V)
{ return (V[0] != 0) ? (T(1) / V[0]) * vect(V) : Vect3_0; }

/// inner (dot) product of #FourVector's
template <typename T>
T operator*(const FourVector<T>& A, const FourVector<T>& B )
{ return A.front() * B.front() - std::inner_product(A.begin() + 1, A.end(), B.begin() + 1, 0); }

/// apply a three-rotation to a FourVector (rotating only the spatial components)
template <typename T>
FourVector<T>& operator*(const ThreeVectorRotation<T>& R, const FourVector<T>& V)
{ return FourVector<T>(V[0], R * vect(V)); }

}
#endif

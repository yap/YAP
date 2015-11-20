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

#ifndef yap_ThreeVector_h
#define yap_ThreeVector_h

#include <algorithm>
#include <array>

namespace yap {

/// \typedef ThreeVector
template <typename T>
using ThreeVector = std::array<T, 3>;

/// add two #ThreeVector's
template <typename T>
T operator+(const ThreeVector<T>& A, const ThreeVector<T>& B)
{
    ThreeVector<T> res;
    std::transform(A.begin(), A.end(), B.begin(), res.begin(), std::operator+);
    return res;
}

/// subtract #ThreeVector from another
template <typename T>
T operator-(const ThreeVector<T>& A, const ThreeVector<T>& B)
{
    ThreeVector<T> res;
    std::transform(A.begin(), A.end(), B.begin(), res.begin(), std::operator-);
    return res;
}

/// inner (dot) product of #ThreeVector's
template <typename T>
T operator*(const ThreeVector<T>& A, const ThreeVector<T>& B )
{ return std::inner_product(A.begin(), A.end(), B.begin(), 0); }

/// outer (cross) product of #ThreeVector's
template <typename T>
ThreeVector<T> cross(const ThreeVector<T>& A, const ThreeVector<T>& B)
{ return {A[2] * B[3] - A[3] * B[2], A[3] * B[1] - A[1] * B[3], A[1] * B[2] - A[2] * B[1]}; }

/// multiply a #ThreeVector by a single element
template <typename T>
ThreeVector<T> operator*(const T& c, const ThreeVector<T>& V)
{
    ThreeVector<T> res;
    std::transform(V.begin(), V.end(), res.begin(), [&](const T& v){return c * v;});
    return res;
}

/// multiply a #ThreeVector by a single element
template <typename T>
ThreeVector<T> operator*(const ThreeVector<T>& V, const T& c)
{ return c * V; }

/// (assignment) multiply a #ThreeVector by a single element
template <typename T>
ThreeVector<T>& operator*=(ThreeVector<T>& V, const T& c)
{
    std::transform(V.begin(), V.end(), V.begin(), [&](const T& v){return c * v;});
    return V;
}

/// \return components of first vector parallel and perpendicular to second
/// \param A Vector to be split
/// \param B Vector defining parallel direction
template <typename T>
std::pair<ThreeVector<T>, ThreeVector<T> > ParallelPerpendicular(const ThreeVector<T>& A, const ThreeVector<T>& B)
{
    ThreeVector<T> parallel = (A * B) / (B * B) * B;
    ThreeVector<T> perpendicular = A - parallel;
    return std::make_pair(parallel, perpendicular);
}

/// \typedef ThreeVectorRotation
/// first element is rotation axis, second is rotation angle
template <typename T>
using ThreeVectorRotation = std::pair<ThreeVector<T>, double>;

/// apply a rotation to a ThreeVector
template <typename T>
ThreeVector<T>& operator*(const ThreeVectorRotation<T>& R, const ThreeVector<T>& V)
{
    auto P = ParallelPerpendicular(V, R.first);
    ThreeVector<T> normal = cross(R.first, P.second);
    return cos(R.second) * P.second
           + sin(R.second) * sqrt(P.second * P.second / normal * normal) * normal
           + P.first;
}

/// (assignment) apply a rotation to a ThreeVector
template <typename T>
ThreeVector<T>& operator*=(ThreeVector<T>& V, const ThreeVectorRotation<T>& R)
{
    V = R * V;
    return V;
}

}
#endif

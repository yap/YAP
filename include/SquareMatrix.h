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

#ifndef yap_SquareMatrix_h
#define yap_SquareMatrix_h

#include "NVector.h"

#include <algorithm>

namespace yap {

/// \typedef SquareMatrix
/// \ingroup VectorAlgebra
template <typename T, size_t N>
using SquareMatrix = NVector<NVector<T, N>, N>;

/// \typedef ThreeMatrix
/// \ingroup VectorAlgebra
template <typename T>
using ThreeMatrix = SquareMatrix<T, 3>;

/// \typedef FourMatrix
/// \ingroup VectorAlgebra
template <typename T>
using FourMatrix = SquareMatrix<T, 4>;

/// unit matrix
template <typename T, size_t N>
SquareMatrix<T, N> unitMatrix()
{
    SquareMatrix<T, N> u;
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
            u[i][j] = (T)(i == j);
    return u;
}

/// matrix multiplication
template <typename T, size_t N>
SquareMatrix<T, N> operator*(const SquareMatrix<T, N> A, const SquareMatrix<T, N> B)
{
    SquareMatrix<T, N> res;
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
            for (size_t n = 0; n < N; ++n)
                res[i][j] += A[i][n] * B[n][j];
    return res;
}

/// multiplication by a single element
template <typename T, size_t N>
SquareMatrix<T, N> operator*(const T& c, const SquareMatrix<T, N>& M)
{
    SquareMatrix<T, N> res = M;
    std::transform(res.begin(), res.end(), res.begin(), [&](const NVector<T, N>& r) {return c * r;});
    return res;
}

/// multiplication by a single element
template <typename T, size_t N>
SquareMatrix<T, N> operator*(const SquareMatrix<T, N>& M, const T& c)
{ return c * M; }

/// multiplication: SquareMatrix * NVector
template <typename T, size_t N>
NVector<T, N> operator*(const SquareMatrix<T, N>& M, const NVector<T, N>& V)
{
    NVector<T, N> res;
    std::transform(M.begin(), M.end(), res.begin(), [&](const NVector<T, N>& m) {return m * V;});
    return res;
}

/// outer product of two #NVector's
template<typename T, size_t N>
SquareMatrix<T, N> outer(const NVector<T, N>& A, const NVector<T, N>& B)
{
    SquareMatrix<T, N> res;
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
            res[i][j] = A[i] * B[j];
    return res;
}

/// skew symmetric matrix formed from a 3 vector
template<typename T>
SquareMatrix<T, 3> skewSymmetric(const NVector<T, 3> V)
{
    return {
        {    0., -V[2],  V[1] },
        {  V[2],    0., -V[0] },
        { -V[1],  V[0],    0. }
    };
}

}
#endif

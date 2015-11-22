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
template <typename T, unsigned N>
using SquareMatrix = std::array<NVector<T, N>, N>;

/// addition
template <typename T, unsigned N>
SquareMatrix<T, N> operator+(const SquareMatrix<T, N>& A, const SquareMatrix<T, N>& B)
{
    SquareMatrix<T, N> res;
    std::transform(A.begin(), A.end(), B.begin(), res.begin(), operator+);
    return res;
}

/// subtraction
template <typename T, unsigned N>
SquareMatrix<T, N> operator-(const SquareMatrix<T, N>& A, const SquareMatrix<T, N>& B)
{
    SquareMatrix<T, N> res;
    std::transform(A.begin(), A.end(), B.begin(), res.begin(), operator-);
    return res;
}

/// multiplication
template <typename T, unsigned N>
SquareMatrix<T, N> operator*(const SquareMatrix<T, N> A, const SquareMatrix<T, N> B)
{
    SquareMatrix<T, N> res;
    for (unsigned i = 0; i < N; ++i)
        for (unsigned j = 0; j < N; ++j)
            for (unsigned n = 0; n < N; ++n)
                res[i][j] += A[i][n] * B[n][j];
    return res;
}

/// multiplication by a single element
template <typename T, unsigned N>
SquareMatrix<T, N> operator*(const T& c, const SquareMatrix<T, N>& M)
{
    SquareMatrix<T, N> res;
    std::transform(M.begin(), M.end(), res.begin(), [&](const NVector<T, N>& r) {return c * r;});
    return res;
}

/// multiplication by a single element
template <typename T, unsigned N>
SquareMatrix<T, N> operator*(const SquareMatrix<T, N>& M, const T& c)
{ return c * M; }

/// multiplication: SquareMatrix * NVector
template <typename T, unsigned N>
NVector<T, N> operator*(const SquareMatrix<T, N>& M, const NVector<T, N>& V)
{
    NVector<T, N> res;
    std::transform(M.begin(), M.end(), res.begin(), [&](const NVector<T, N>& m) {return m * V;});
    return res;
}

/// outer product of two #NVector's
template<typename T, unsigned N>
SquareMatrix<T, N> outer(const NVector<T, N>& A, const NVector<T, N>& B)
{
    SquareMatrix<T, N> res;
    for (unsigned i = 0; i < N; ++i)
        for (unsigned j = 0; j < N; ++j)
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

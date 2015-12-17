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

#ifndef yap_Matrix_h
#define yap_Matrix_h

#include <array>

namespace yap {

/// \typedef Matrix
/// \ingroup VectorAlgebra
/// \param R number of rows
/// \param C number of columns
template <typename T, size_t R, size_t C>
using Matrix = std::array<std::array<T, C>, R>;

/// \typedef SquareMatrix
/// \ingroup VectorAlgebra
template <typename T, size_t N>
using SquareMatrix = Matrix<T, N, N>;

/// \typedef ThreeMatrix
/// \ingroup VectorAlgebra
template <typename T>
using ThreeMatrix = SquareMatrix<T, 3>;

/// \typedef FourMatrix
/// \ingroup VectorAlgebra
template <typename T>
using FourMatrix = SquareMatrix<T, 4>;

/// \return string
template <typename T, size_t R, size_t C>
std::string to_string(const Matrix<T, R, C>& M)
{
    if (M.empty() or M.front().empty())
        return "(empty)";
    std::string s = "(";
    for (const auto& row : M) {
        s += "(";
        for (const auto& elt : row)
            s += std::to_string(elt) + ", ";
        s.erase(s.size() - 2, 2);
        s += "), ";
    }
    s.erase(s.size() - 2, 2);
    s += ")";
    return s;
}

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

/// transpose a matrix
template <typename T, size_t R, size_t C>
Matrix<T, C, R> transpose(const Matrix<T, R, C>& M)
{
    Matrix<T, C, R> res;
    for (size_t r = 0; r < R; ++r)
        for (size_t c = 0; c < C; ++c)
            res[c][r] = M[r][c];
}

/// uniary minus
template <typename T, size_t R, size_t C>
constexpr Matrix<T, C, R> operator-(const Matrix<T, R, C>& M)
{ return T(-1) * M; }

/// matrix multiplication
template <typename T, size_t R, size_t K, size_t C>
std::enable_if_t < (R != 1) or (C != 1), Matrix<T, R, C> >
operator*(const Matrix<T, R, K> A, const Matrix<T, K, C> B)
{
    Matrix<T, R, C> res;
    for (size_t r = 0; r < R; ++r)
        for (size_t c = 0; c < C; ++c)
            for (size_t k = 0; k < K; ++k)
                res[r][c] += A[r][k] * B[k][c];
    return res;
}

/// matrix multiplication yielding single value
template <typename T, size_t K>
T operator*(const Matrix<T, 1, K> A, const Matrix<T, K, 1> B)
{
    T res;
    for (size_t k = 0; k < K; ++k)
        res += A[0][k] * B[k][0];
    return res;
}

/// assignment by multiplication by a single element
template <typename T, size_t R, size_t C>
Matrix<T, R, C>& operator*=(Matrix<T, R, C>& M, const T& c)
{
    for (auto& row : M)
        for (auto& elt : row)
            elt *= c;
    return M;
}

/// multiplication by a single element
template <typename T, size_t R, size_t C>
Matrix<T, R, C> operator*(const T& c, const Matrix<T, R , C>& M)
{ auto m = M; m *= c; return m; }

/// multiplication by a single element
template <typename T, size_t R, size_t C>
Matrix<T, R, C> operator*(const Matrix<T, R, C>& M, const T& c)
{ return c * M; }

/// addition assignment
template <typename T, size_t R, size_t C>
Matrix<T, R, C>& operator+=(Matrix<T, R, C>& lhs, const Matrix<T, R, C>& rhs)
{
    for (size_t r = 0; r < R; ++r)
        for (size_t c = 0; c < C; ++c)
            lhs[r][c] += rhs[r][c];
    return lhs;
}

/// addition
template <typename T, size_t R, size_t C>
Matrix<T, R, C> operator+(const Matrix<T, R, C>& lhs, const Matrix<T, R, C>& rhs)
{ auto m = lhs; m += rhs; return m; }

/// subtraction assignment
template <typename T, size_t R, size_t C>
Matrix<T, R, C>& operator-=(Matrix<T, R, C>& lhs, const Matrix<T, R, C>& rhs)
{
    for (size_t r = 0; r < R; ++r)
        for (size_t c = 0; c < C; ++c)
            lhs[r][c] -= rhs[r][c];
    return lhs;
}

/// subtraction
template <typename T, size_t R, size_t C>
Matrix<T, R, C> operator-(const Matrix<T, R, C>& lhs, const Matrix<T, R, C>& rhs)
{ auto m = lhs; m -= rhs; return m; }

}
#endif

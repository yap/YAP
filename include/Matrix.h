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

#include "Exceptions.h"
#include "MathUtilities.h"

#include <array>
#include <string>
#include <vector>

namespace yap {

/// \class Matrix
/// \ingroup VectorAlgebra
/// \param R number of rows
/// \param C number of columns
template <typename T, size_t R, size_t C>
class Matrix : public std::array<std::array<T, C>, R>
{
public:
    /// Constructor
    constexpr Matrix(const std::array<std::array<T, C>, R>& m) noexcept : std::array<std::array<T, C>, R>(m) {}

    /// Default constructor;
    /// produces a zero matrix
    Matrix()
    {
        for (size_t i = 0; i < this->size(); ++i)
            this->at(i).fill(T(0));
    }

    /// Use std::array's assignment operators
    using std::array<std::array<T, C>, R>::operator=;
};

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

/// zero square matrix
template <typename T, size_t N>
SquareMatrix<T, N> zeroMatrix()
{
    SquareMatrix<T, N> u;
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
            u[i][j] = (T)(0);
    return u;
}

/// zero matrix
template <typename T, size_t R, size_t C>
Matrix<T, R, C> zeroMatrix()
{
    Matrix<T, R, C> u;
    for (size_t r = 0; r < R; ++r)
        for (size_t c = 0; c < C; ++c)
            u[r][c] = (T)(0);
    return u;
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

/// diagonal matrix
template <typename T, size_t N>
SquareMatrix<T, N> diagonalMatrix(std::array<T, N> d)
{
    SquareMatrix<T, N> D = zeroMatrix<T, N, N>();
    for (size_t i = 0; i < N; ++i)
        D[i][i] = d[i];
    return D;
}

/// transpose a matrix
template <typename T, size_t R, size_t C>
Matrix<T, C, R> transpose(const Matrix<T, R, C>& M)
{
    Matrix<T, C, R> res;
    for (size_t r = 0; r < R; ++r)
        for (size_t c = 0; c < C; ++c)
            res[c][r] = M[r][c];
    return res;
}

/// unary minus
template <typename T, size_t R, size_t C>
constexpr Matrix<T, C, R> operator-(const Matrix<T, R, C>& M)
{ return T(-1) * M; }

/// matrix multiplication
template <typename T, size_t R, size_t K, size_t C>
typename std::enable_if < (R != 1) or (C != 1), Matrix<T, R, C> >::type
const operator*(const Matrix<T, R, K>& A, const Matrix<T, K, C>& B)
{
    Matrix<T, R, C> res = zeroMatrix<T, R, C>();
    for (size_t r = 0; r < R; ++r)
        for (size_t c = 0; c < C; ++c)
            for (size_t k = 0; k < K; ++k)
                res[r][c] += A[r][k] * B[k][c];
    return res;
}

/// matrix multiplication yielding single value
template <typename T, size_t K>
const T operator*(const Matrix<T, 1, K>& A, const Matrix<T, K, 1>& B)
{
    T res(0);
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
const Matrix<T, R, C> operator*(const T& c, const Matrix<T, R, C>& M)
{ auto m = M; m *= c; return m; }

/// multiplication by a single element
template <typename T, size_t R, size_t C>
const Matrix<T, R, C> operator*(const Matrix<T, R, C>& M, const T& c)
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
const Matrix<T, R, C> operator+(const Matrix<T, R, C>& lhs, const Matrix<T, R, C>& rhs)
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
const Matrix<T, R, C> operator-(const Matrix<T, R, C>& lhs, const Matrix<T, R, C>& rhs)
{ auto m = lhs; m -= rhs; return m; }

/// Minor matrix
/// \param r row
/// \param c column
template <typename T, size_t N>
const SquareMatrix < T, N - 1 > minor_matrix(const SquareMatrix<T, N>& M, size_t r, size_t c)
{
    SquareMatrix < T, N - 1 > m;
    for (size_t i = 0; i < N; ++i)
        if (i != r)
            for (size_t j = 0; j < N; ++j)
                if (j != c)
                    m[i - (size_t)(i > r)][j - (size_t)(j > c)] = M[i][j];
    return m;
}

/// Minor
/// \param r row
/// \param c column
template <typename T, size_t N>
const T minor_det(const SquareMatrix<T, N>& M, size_t r, size_t c)
{ return det(minor_matrix(M, r, c)); }

/// cofactor
/// \param r row
/// \param c column
template <typename T, size_t N>
const T cofactor(const SquareMatrix<T, N>& M, size_t r, size_t c)
{ return pow_negative_one(r + c) * M[r][c] * minor_det(M, r, c); }

/// Determinant
/// This algorithm is simple, and should not be used for large matrices
template <typename T, size_t N>
const T det(SquareMatrix<T, N> M)
{
    switch (N) {
        case 0:
            throw exceptions::Exception("zero-size matric", "determinant");
        case 1:
            return M[0][0];
        case 2:
            return M[0][0] * M[1][1] - M[0][1] * M[1][0];
        case 3:
            return M[0][0] * M[1][1] * M[2][2]
                   +  M[0][1] * M[1][2] * M[2][0]
                   +  M[0][2] * M[1][0] * M[2][1]
                   -  M[0][2] * M[1][1] * M[2][0]
                   -  M[0][0] * M[1][2] * M[2][1]
                   -  M[0][1] * M[1][0] * M[2][2];
        default:
            T d(0);
            for (size_t i = 0; i < N; ++i)
                d += cofactor(M, i, 0);
            return d;
    }
}

/// diagonal minor matrix
template <typename T, size_t M, size_t N>
const SquareMatrix<T, M> diagonal_minor(const SquareMatrix<T, N> m, std::vector<size_t> indices)
{
    if (indices.size() != M)
        throw exceptions::Exception("wrong number of indices provided", "diagonal_minor");

    SquareMatrix<T, M> dm;
    for (size_t r = 0; r < indices.size(); ++r)
        for (size_t c = 0; c < indices.size(); ++c)
            dm[indices[r]][indices[c]] = m[r][c];
    return dm;
}

/// trace
template <typename T, size_t N>
const T trace(const SquareMatrix<T, N>& M)
{
    T t(0);
    for (size_t i = 0; i < N; ++i)
        t += M[i][i];
    return t;
}

}
#endif

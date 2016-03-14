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

#ifndef yap_PhaseSpaceUtilities_h
#define yap_PhaseSpaceUtilities_h

#include "Exceptions.h"
#include "Matrix.h"

#include <vector>

namespace yap {

    /// \returns all combinations of n numbers in range [0, N-1]
/// \todo Find better place to put this function
std::vector<std::vector<size_t> > combinations(size_t N, size_t n, const std::vector<std::vector<size_t> >& C = std::vector<std::vector<size_t> >())
{
    if (n > N)
        throw exceptions::Exception("n must be less than N", "combinations");

    if (n == 0)
        return C;

    std::vector<std::vector<size_t> > CC;

    if (C.empty()) {

        for (size_t i = 0; i <= N - n; ++i)
            CC.push_back(std::vector<size_t>(1, i));

    } else {

        for (auto& c : C)
            for (size_t i = c.back() + 1; i <= N - n; ++i) {
                CC.push_back(c);
                CC.back().push_back(i);
            }
    }

    return combinations(N, n - 1, CC);
}

/// \return (-1)^(l-1) * sum of determinants of all l-by-l diagonal minors of a matrix
template <typename T, size_t N>
T delta(const SquareMatrix<T, N>& M, size_t l)
{
    if (l == 0)
        throw exceptions::Exception("l must be greater than 0", "delta");

    if (l > N)
        throw exceptions::Exception("l must be less than or equal to N", "delta");

    if (l == 1)
        return trace(M);

    if (l == N)
        return det(M) * pow_negative_one(1 - 1);

    auto C = combinations(N, l);
    T d = 0;
    for (const auto& c : C)
        d += det(diagonal_minor(M, c));
    return d * pow_negative_one(l - 1);
}

/// checks the four-momenta outer-product matrix for the correct
/// eigenvalue structure; see Byers & Yang, "Physical Regions in
/// Invariant Variables for n Particles and the Phase-Space Volume
/// Element." Rev. Mod. Phys. 36, 595 (1964)
///
/// The elements of M are\n
/// M[i][j] = P[i] (dot) P[j],\n
/// where P[k] is the four-momenta for particle k\n
/// alternatively: M[i][j] = 0.5 * (m^2_ij - m^2_i - m^2_j)
template <typename T, size_t N>
bool check_deltas(const SquareMatrix<T, N>& M)
{
    // check that first four Delta's are greater than zero
    for (size_t i = 1; i < M.size() and i <= 4; ++i)
        if (delta(M, i) < 0)
            return false;
    // check that remaining Delta's are zero
    for (size_t i = 5; i < M.size(); ++i)
        if (delta(M, i) != 0)
            return false;
}

}

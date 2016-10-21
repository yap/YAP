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

#ifndef yap__Group_h
#define yap__Group_h

#include <algorithm>
#include <vector>

namespace yap {

/// Groups members of a vector into subgroups by equality under comparitor
/// \tparam T argument type
/// \tparam Compare attribute type
/// \param V vector to group
/// \param C Comparitor functor
template <typename container, typename Compare>
std::vector<std::vector<typename container::value_type> > group(const container& K, Compare C)
{
    /// create sortable container (vector) from generic container
    std::vector<typename container::value_type> V(K.begin(), K.end());
    /// sort by C
    std::stable_sort(V.begin(), V.end(), C);
    /// Group by C
    std::vector<std::vector<typename container::value_type> > W;
    for (auto it = V.begin(); it != V.end();) {
        // find range equal to *it according to C
        auto its = std::equal_range(it, V.end(), *it, C);
        // add to W new vector for this group
        W.push_back(std::vector<typename container::value_type>(its.first, its.second));
        // continue beyond this group
        it = its.second;
    }
    return W;
}

/// Groups members of a vector into subgroups by equality under comparitor
/// \tparam container Container type
/// \tparam Compare attribute type
/// \param K container 
/// \param c1 first Comparitor functor
/// \param c2 first Comparitor functor
/// \param C further Comparitor functors
template <typename container, typename Compare1, typename Compare2, typename ... Compares>
std::vector<std::vector<typename container::value_type> > group(const container& K, Compare1 c1, Compare2 c2, Compares ... C)
{
    // group by c1
    auto V = group(K, c1);
    // run recursively on new groups
    std::vector<std::vector<typename container::value_type> > W;
    for (const auto& v : V) {
        auto vg = group(v, c2, C...);
        W.insert(W.end(), vg.begin(), vg.end());
    }
    return W;
}

}

#endif

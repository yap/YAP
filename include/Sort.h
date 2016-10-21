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

#ifndef yap__Sort_h
#define yap__Sort_h

#include <algorithm>
#include <vector>

namespace yap {

/// sort elements in range [first, last) by less-than functor c1
/// \tparam ItType iterator type
/// \tparam Compare less-than functor type
/// \param first iterator to start of range
/// \param last iterator beyond end of range
/// \param c less-than functor object for comparison
template <typename ItType, typename Compare>
void lexicographically_sort(ItType first, ItType last, Compare c)
{ std::stable_sort(first, last, c); }

/// sort elements in range [first, last) lexicographically by c1, c2, C...
/// \tparam ItType iterator type
/// \tparam Compare1 first less-than functor type
/// \tparam Compare2 second less-than functor type
/// \tparam Compares further less-than functor type
/// \param first iterator to start of range
/// \param last iterator beyond end of range
/// \param c1 first less-than functor object for comparison
/// \param c2 second less-than functor object for comparison
/// \param C further less-than functor object for comparison
template <typename ItType, typename Compare1, typename Compare2, typename ... Compares>
void lexicographically_sort(ItType first, ItType last, Compare1 c1, Compare2 c2, Compares ... C)
{
    
    // sort by c1
    lexicographically_sort<ItType, Compare1>(first, last, c1);
    
    // then group by c1 and recurse on groupings
    for (auto it = first; it != last; ) {
        // find range equal to *it
        auto its = std::equal_range(it, last, *it, c1);
        // recurse on this range with the next comparison functor object
        lexicographically_sort(its.first, its.second, c2, C...);
        // continue beyond this group
        it = its.second;
    }
}

/// sort a container of objects by lexicographically applying the
/// provided sorters, which are binary predicates performing a less
/// than function
/// \tparam container Container type to sort
/// \tparam Compare1 first less-than functor type
/// \tparam Compares further less-than functor type
/// \param K container to be sorted
/// \param c1 first less-than functor object for comparison
/// \param C further less-than functor object for comparison
template <typename container, typename Compare1, typename ... Compares>
std::vector<typename container::value_type> sort(container K, Compare1 c1, Compares ... C)
{
    // create vector that can be sorted
    std::vector<typename container::value_type> V(K.begin(), K.end());
    // sort it
    lexicographically_sort(V.begin(), V.end(), c1, C...);
    // return it
    return V;
}



}

#endif

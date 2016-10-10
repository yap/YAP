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

#ifndef yap__Filter_h
#define yap__Filter_h

#include "Exceptions.h"

#include <algorithm>

namespace yap {

/// Throws if container's size is not 1
/// \return lone element in container
template <typename container>
typename container::value_type lone_elt(container& C)
{
    if (C.size() != 1)
        throw yap::exceptions::Exception("Container size not 1 (" + std::to_string(C.size()) + ")", "lone_elt");
    return *C.begin();
}

/// Throws if container's size is not 1
/// \return lone element in container
template <typename container>
typename container::value_type lone_elt(container&& C)
{
    if (C.size() != 1)
        throw yap::exceptions::Exception("Container size not 1 (" + std::to_string(C.size()) + ")", "lone_elt");
    return *C.begin();
}

/// filter through only the members of a container that evaluate to
/// true with all of the predicates given
template <typename container, typename First, typename Second, typename ... UnaryPredicates>
container filter(container C, First p1, Second p2, UnaryPredicates ... P)
{
    for (auto it = C.begin(); it != C.end(); ) {
        if (p1(**it))
            ++it;
        else
            it = C.erase(it);
    }
    return filter<container, Second, UnaryPredicates...>(C, p2, P...);
}

/// filter through only the members of a container that evaluate to
/// true with all of the predicates given
template <typename container, typename First>
container filter(container C, First p)
{
    for (auto it = C.begin(); it != C.end(); ) {
        if (p(**it))
            ++it;
        else
            it = C.erase(it);
    }
    return C;
}

}

#endif

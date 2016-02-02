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

#ifndef yap_container_utils_h
#define yap_container_utils_h

#include <algorithm>

/// check if two containers overlap (no sorting necessary).
/// BinaryPredicate must have signature bool(const InputIt1::type&, const InputIt2::type&)
template <class InputIt1, class InputIt2, class BinaryPredicate>
bool overlap(InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2, BinaryPredicate p)
{
    using T1 = typename std::iterator_traits<InputIt1>::value_type;
    using T2 = typename std::iterator_traits<InputIt2>::value_type;
    return std::any_of(first1, last1, [&](const T1 & a) { return std::any_of(first2, last2, [&](const T2 & b) {return p(a, b);}); });
}

/// check if two containers are disjoint (no sorting necessary).
/// BinaryPredicate must have signature bool(const InputIt1::type&, const InputIt2::type&)
template <class InputIt1, class InputIt2, class BinaryPredicate>
bool disjoint(InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2, BinaryPredicate p)
{
    using T1 = typename std::iterator_traits<InputIt1>::value_type;
    using T2 = typename std::iterator_traits<InputIt2>::value_type;
    return std::none_of(first1, last1, [&](const T1 & a) { return std::none_of(first2, last2, [&](const T2 & b) {return p(a, b);}); });
}

/// check if first container contains second (no sorting necessary).
/// BinaryPredicate must have signature bool(const InputIt1::type&, const InputIt2::type&)
template <class InputIt1, class InputIt2, class BinaryPredicate>
bool contains(InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2, BinaryPredicate p)
{
    using T1 = typename std::iterator_traits<InputIt1>::value_type;
    using T2 = typename std::iterator_traits<InputIt2>::value_type;
    return std::all_of(first2, last2, [&](const T2 & b) { return std::any_of(first1, last1, [&](const T1 & a) {return p(a, b);}); });
}

/// useful for removing duplicates from an unsorted vector (preserving order of first occurance)
/// \tparam InputIt iterator type
// \tparam Compare less-than comparitor
// \tparam BinaryPredicate equality comparitor
template <class InputIt/*, class Compare, class BinaryPredicate*/>
InputIt ordered_unique(InputIt first, InputIt last/*, Compare lt, BinaryPredicate eq*/)
{
    // make vector of iterators
    std::vector<InputIt> v;
    v.reserve(std::distance(first, last));
    for (InputIt i = first; i != last; ++i)
        v.push_back(i);
    // sort it
    std::sort(v.begin(), v.end(), [](const InputIt & A, const InputIt & B) {return /*lt(*A, *B)*/ *A < *B;});
    // apply unique to it
    v.erase(std::unique(v.begin(), v.end(), [](const InputIt & A, const InputIt & B) {return /*eq(*A, *B)*/ *A == *B;}), v.end());
    std::sort(v.begin(), v.end());

    size_t j = 0;
    for (InputIt i = first; i != last && j != v.size(); ++i)
        if (i == v[j]) {
            std::iter_swap(i, first);
            ++j;
            ++first;
        }

    return first;
}

/// \todo allow passing of binary predicate with defaulting to lambda

/// check if two vectors overlap
template <typename T>
bool overlap(const std::vector<T>& A, const std::vector<T>& B)
{ return overlap(A.begin(), A.end(), B.begin(), B.end(), [](const T & a, const T & b) {return a == b;}); }

/// check if two vectors are disjoint
template <typename T>
bool disjoint(const std::vector<T>& A, const std::vector<T>& B)
{ return disjoint(A.begin(), A.end(), B.begin(), B.end(), [](const T & a, const T & b) {return a == b;}); }

/// check if the first vector contains the second one
template <typename T>
bool contains(const std::vector<T>& A, const std::vector<T>& B)
{ return contains(A.begin(), A.end(), B.begin(), B.end(), [](const T & a, const T & b) {return a == b;}); }


#endif

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

#ifndef yap__FiltersFwd_h
#define yap__FiltersFwd_h

#include <memory>
#include <set>

namespace yap {

template <typename container>
typename container::value_type lone_elt(container& C);

template <typename container>
typename container::value_type lone_elt(container&& C);

template <typename T, typename Last, typename ... UnaryPredicates>
    std::set<std::shared_ptr<T> > filter(const std::set<std::shared_ptr<T> >& S, Last p, UnaryPredicates ... P);

template <typename T>
    const std::set<std::shared_ptr<T> >& filter(const std::set<std::shared_ptr<T> >& S);

template <class F, typename T>
const bool by_ptr(const F& f, const std::shared_ptr<T>& ptr);

template <class F, typename T>
const bool by_ptr(const F& f, const T* ptr);

struct filter_decay_tree;
struct filter_free_amplitude;
struct filter_decay_channel;
struct filter_spin_amplitude;
struct filter_particle;
struct filter_decaying_particle;

class has_decay_channel;
class has_spin_amplitude;
class has_decay_tree;
class has_free_amplitude;

class to;
class from;
class l_equals;
class m_equals;
class is_named;

}

#endif

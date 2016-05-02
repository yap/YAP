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
/// Contains forward declarations only

#ifndef yap_WeakPtrCacheFwd_h
#define yap_WeakPtrCacheFwd_h

namespace yap {

template <class T>
class WeakPtrCache;

// The following cannot be forward declared since they are nested
// inside WeakPtrCache
/* using WeakPtrCache::type = T; */
/* using WeakPtrCache::shared_ptr_type = std::shared_ptr<T>; */
/* using WeakPtrCache::weak_ptr_type = std::weak_ptr<T>; */
/* using WeakPtrCache::cache_type = std::set<weak_ptr_type, std::owner_less<weak_ptr_type> >; */

}

#endif

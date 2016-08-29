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

#ifndef yap_ParticleCombinationFwd_h
#define yap_ParticleCombinationFwd_h

#include <functional>
#include <map>
#include <memory>
#include <set>
#include <vector>

namespace yap {

class ParticleCombination;

/// \typedef ParticleCombinationEqualTo
/// std::function to compare two ParticleCombination's
using ParticleCombinationEqualTo = std::function<bool(const std::shared_ptr<const ParticleCombination>&, const std::shared_ptr<const ParticleCombination>&)>;

/// \typedef ParticleCombinationVector
using ParticleCombinationVector = std::vector<std::shared_ptr<ParticleCombination> >;

/// \typedef ParticleCombinationSet
using ParticleCombinationSet = std::set<std::shared_ptr<ParticleCombination> >;

/// \typedef ParticleCombinationMap
/// \tparam T Object to store in map, with shared_ptr to ParticleCombination as key
template<typename T>
using ParticleCombinationMap = std::map<std::shared_ptr<ParticleCombination>, T,
      std::owner_less<std::shared_ptr<ParticleCombination> > >;

}

#endif

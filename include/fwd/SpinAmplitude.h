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

#ifndef yap_SpinAmplitudeFwd_h
#define yap_SpinAmplitudeFwd_h

#include <map>
#include <memory>
#include <vector>

namespace yap {

class SpinAmplitude;

// The following cannot be forward declared since they are nested
// inside SpinAmplitude
/* using SpinAmplitude::SpinProjectionPair = std::array<int, 2>; */
/* using SpinAmplitude::AmplitudeSubmap = std::map<SpinProjectionPair, std::shared_ptr<ComplexCachedDataValue> >; */
/* using SpinAmplitude::AmplitudeMap = std::map<int, AmplitudeSubmap>; */

/// \typedef SpinAmplitudeVector
using SpinAmplitudeVector = std::vector<std::shared_ptr<SpinAmplitude> >;

/// \typedef SpinAmplitudeMap
/// \tparam T Object to store in map, with shared_ptr to SpinAmplitude as key
template<typename T>
using SpinAmplitudeMap = std::map<std::shared_ptr<SpinAmplitude>, T,
      std::owner_less<std::shared_ptr<SpinAmplitude> > >;

}

#endif

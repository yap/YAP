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
/// contains forward declaration only

#ifndef yap_PhaseSpaceFactorFwd_h
#define yap_PhaseSpaceFactorFwd_h

#include "fwd/SpinAmplitude.h"

#include <map>
#include <memory>

namespace yap {

class PhaseSpaceFactor;
class PhaseSpaceFactorFactory;

/// \typedef PhaseSpaceFactorMap
using PhaseSpaceFactorMap = std::map<std::shared_ptr<SpinAmplitude>, std::shared_ptr<PhaseSpaceFactor> >;

extern std::shared_ptr<PhaseSpaceFactorFactory> DefaultPHSPFactory;

}

#endif

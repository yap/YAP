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

#ifndef yap_ReportsParticleCombinations_h
#define yap_ReportsParticleCombinations_h

#include "ParticleCombination.h"

/// \name ReportsParticleCombinations
/// \brief Base class for all classes that report a list of particleCombinations
/// \author Daniel Greenwald

namespace yap {

class ReportsParticleCombinations
{

public:

    /// \return vector of ParticleCombinations
    virtual ParticleCombinationVector particleCombinations() const = 0;

    /// Default constructor
    ReportsParticleCombinations() {}

    /// virtual destructor
    virtual ~ReportsParticleCombinations() = default;

    /// default copy constructor
    ReportsParticleCombinations(const ReportsParticleCombinations& other) = default;

    /// default move constructor
    ReportsParticleCombinations(ReportsParticleCombinations&& other) = default;

    /// default copy assignment operator
    ReportsParticleCombinations& operator=(const ReportsParticleCombinations& rhs) = default;

    /// default move assignment operator
    ReportsParticleCombinations& operator=(ReportsParticleCombinations&& rhs) = default;

};

}

#endif

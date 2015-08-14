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

#ifndef yap_ParticleCombination_h
#define yap_ParticleCombination_h

#include "ParticleIndex.h"

#include <memory>
#include <vector>

namespace yap {

/// \class ParticleCombination
/// \brief Stores combinations of ParticleIndex types
/// \author Johannes Rauch, Daniel Greenwald

class ParticleCombination
{
public:

    ParticleCombination() = default;

    /// Final-state-particle constructor
    ParticleCombination(ParticleIndex index);

    /// \name Getters
    /// @{

    /// Get vector of indices
    const std::vector<ParticleIndex>& indices() const
    { return Indices_; }

    /// Get vector of daughters as weak_ptr's
    std::vector<std::weak_ptr<ParticleCombination> > daughters() const
    { return std::vector<std::weak_ptr<ParticleCombination> >(Daughters_.begin(), Daughters_.end()); }

    /// @}

    /// Add daughter ParticleCombination
    /// \param daughter Shared pointer to ParticleCombination object representing a daughter
    /// \return Success of action
    bool addDaughter(std::shared_ptr<ParticleCombination> daughter);

    /// Checks consistency of combination
    /// by checking for absence of duplicate entries
    bool consistent() const;

    /// equality operator
    friend bool operator==(const ParticleCombination& A, const ParticleCombination& B);

protected:
    std::vector<std::shared_ptr<ParticleCombination> > Daughters_;
    std::vector<ParticleIndex> Indices_;
};

}

#endif

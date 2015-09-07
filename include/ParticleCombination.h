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
#include <set>
#include <vector>

namespace yap {

/// \class ParticleCombination
/// \brief Stores combinations of ParticleIndex types
/// \author Johannes Rauch, Daniel Greenwald

class ParticleCombination
{
public:

    /// \todo Private constructor to force use of static creation functions?
    ParticleCombination();

    /// Final-state-particle constructor
    ParticleCombination(ParticleIndex index);

    /// Resonance particle constructor
    ParticleCombination(std::vector<std::shared_ptr<ParticleCombination> > c);

    /// \name Getters
    /// @{

    /// Get vector of indices
    const std::vector<ParticleIndex>& indices() const
    { return Indices_; }

    /// Get vector of daughters as weak_ptr's
    //std::vector<std::weak_ptr<ParticleCombination> > daughters() const
    //{ return std::vector<std::weak_ptr<ParticleCombination> >(Daughters_.begin(), Daughters_.end()); }

    /// Get vector of daughters
    const std::vector<std::shared_ptr<ParticleCombination> >& daughters() const
    { return Daughters_;}

    /// get parent
    const ParticleCombination* parent() const
    { return Parent_; }

    /// @}

    /// Add daughter ParticleCombination
    /// \param daughter Shared pointer to ParticleCombination object representing a daughter
    /// \return Success of action
    bool addDaughter(std::shared_ptr<ParticleCombination> daughter);

    /// Checks consistency of combination
    /// by checking for absence of duplicate entries
    bool consistent() const;

    /// cast into string
    operator std::string() const;

    /// check if this and B share one or more ParticleIndex's
    bool sharesIndices(std::shared_ptr<ParticleCombination> B);


    /// meh
    void setParents();

    /// equality operator
    friend bool operator==(const ParticleCombination& A, const ParticleCombination& B);

    /// inequality operator
    friend bool operator!=(const ParticleCombination& A, const ParticleCombination& B)
    { return !(A == B); }

protected:
    ParticleCombination* Parent_;
    std::vector<std::shared_ptr<ParticleCombination> > Daughters_;
    std::vector<ParticleIndex> Indices_;


/// \name Static methods for creating/retrieving ParticleCombination's
/// @{

// Following code is for managing unique shared pointers for particle
// combinations across all of YAP

public:

    /// return existing shared_ptr for final-state-particle ParticleCombination, if exists; otherwise creates and returns
    /// \param i ParticleIndex for FSP
    static std::shared_ptr<ParticleCombination> uniqueSharedPtr(std::shared_ptr<ParticleCombination> pc);

    /// return existing shared_ptr for final-state-particle ParticleCombination, if exists; otherwise creates and returns
    /// \param i ParticleIndex for FSP
    static std::shared_ptr<ParticleCombination> uniqueSharedPtr(ParticleIndex i);

    /// return existing shared_ptr for ParticleCombination, if exists; otherwise creates and returns
    /// \param c vector of shared_ptr's to ParticleCombination objects describing new ParticleCombination
    static std::shared_ptr<ParticleCombination> uniqueSharedPtr(std::vector<std::shared_ptr<ParticleCombination> > c);

    static const std::set<std::shared_ptr<ParticleCombination> >& particleCombinationSet()
    { return ParticleCombinationSet_; }

    static void printParticleCombinationSet();

private:

    /// Static set of all particle combinations created throughout code
    static std::set<std::shared_ptr<ParticleCombination> > ParticleCombinationSet_;

/// @}
};


}

#endif

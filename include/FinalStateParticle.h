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

#ifndef yap_FinalStateParticle_h
#define yap_FinalStateParticle_h

#include "Amp.h"
#include "Constants.h"
#include "Particle.h"
#include "ParticleIndex.h"

#include <memory>

namespace yap {

class ParticleCombination;

/// \class FinalStateParticle
/// \brief Class representing a final-state particle
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Particle

class FinalStateParticle : public Particle
{
public:

    /// Constructor
    /// \param indices index or indices (if there are identical final state particles) that this particle has in the DataPoint
    FinalStateParticle(const QuantumNumbers& q, double mass, std::string name, int pdg, std::vector<ParticleIndex>& indices);

    /// \return 1 + 0i
    virtual const Amp& amplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc) override
    { static Amp a = Complex_1; return a; }

    virtual bool consistent() const override;

    /// \return PDG code indicating particle type
    int pdgCode() const
    { return PDGCode_; }

    /// \return list of all ParticleCombinations
    std::vector<std::shared_ptr<ParticleCombination> > particleCombinations() const
    { return SymmetrizationIndices_; }

    // for internal use only
    virtual void setSymmetrizationIndexParents() override;

private:

    /// add symmetrizationIndex to SymmetrizationIndices_
    void addSymmetrizationIndex(std::shared_ptr<ParticleCombination> c);

    /// PDG code of the particle
    int PDGCode_;

    std::vector<std::shared_ptr<ParticleCombination> > SymmetrizationIndices_;

};

}

#endif

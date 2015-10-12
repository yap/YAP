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

#ifndef yap_ParticleFactory_h
#define yap_ParticleFactory_h

#include "ParticleIndex.h"
#include "QuantumNumbers.h"

#include <limits>
#include <map>
#include <memory>
#include <vector>

namespace yap {

class DecayingParticle;
class FinalStateParticle;
class InitialStateParticle;
class MassShape;
class Particle;
class Resonance;

/// \struct ParticleTableEntry
/// \brief Data container for storing particle information in database
/// \author Johannes Rauch, Daniel Greenwald
struct ParticleTableEntry : public QuantumNumbers {
    ParticleTableEntry(int pdg = 0, std::string name = "", QuantumNumbers q = QuantumNumbers(), double mass = -1, std::vector<double> parameters = {});
    bool consistent() const override;
    double mass() const
    { return MassShapeParameters_.at(0); }
    int PDG_;
    std::string Name_;
    std::vector<double> MassShapeParameters_; ///< first entry should always be nominal mass
};

/// \class ParticleFactory
/// \brief Factory class for easy creation of Particle objects from PDG codes.
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Particle

class ParticleFactory
{
public:

    /// Constructor
    /// \param pdlFile Path to a pdl file like used by EvtGen
    ParticleFactory(const std::string pdlFile);

    /// Create a FinalStateParticle from a PDG code
    /// \param PDG PDG code of particle to create
    /// \return HelicityStates object for new final state particle
    std::shared_ptr<FinalStateParticle> createFinalStateParticle(int PDG, std::vector<ParticleIndex> indices);

    /// Create an InitialStateParticle from a PDG code and a MassShape
    /// \param PDG PDG code of particle to create
    /// \param radialSize radial size of particle to create [GeV^-1]
    /// \return pointer to new InitialStateParticle object
    std::shared_ptr<InitialStateParticle> createInitialStateParticle(int PDG, double radialSize);

    /// Create a Resonance from a PDG code and a MassShape
    /// \param PDG PDG code of particle to create
    /// \param radialSize Radial size of particle to create [GeV^-1]
    /// \param massShape Pointer to MassShape object describing resonance
    /// \return pointer to new Resonance object
    std::shared_ptr<Resonance> createResonance(int PDG, double radialSize, std::shared_ptr<MassShape> massShape);

    /// \name Particle table access
    /// @{

    /// get ParticleTableEntry from #particleTable_ with safety checks
    /// \param PDG pdg code labeling particle table entry
    const ParticleTableEntry& particleTableEntry(int PDG) const;

    /// add ParticleTableEntry to #particleTable_
    /// \param entry ParticleTableEntry to add to #particleTable_
    /// \return Success of action
    bool addParticleTableEntry(ParticleTableEntry entry);

    /// @}

private:
    /// read pdl file and fill #particleTable_
    void readPDT(const std::string fname);

    /// maps PDGCodes to ParticleTableEntry's
    std::map<int, ParticleTableEntry> particleTable_;

};

}

#endif

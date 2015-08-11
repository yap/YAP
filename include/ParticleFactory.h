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

#include "FinalStateParticle.h"
#include "InitialStateParticle.h"
#include "MassShape.h"
#include "QuantumNumbers.h"
#include "Resonance.h"

namespace yap {

class Particle;

/// \class Particle
/// \brief Factory class for easy creation of Particle objects from PDG codes.
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Particle

class ParticleFactory
{
public:

    /// Constructor
    ParticleFactory()
    {}

    /// Create a FinalStateParticle from a PDG code
    /// \param PDG PDG code of particle to create
    /// \return pointer to new FinalStateParticle object
    FinalStateParticle* createFinalStateParticle(int PDG);

    /// Create an InitialStateParticle from a PDG code and a MassShape
    /// \param PDG PDG code of particle to create
    /// \param radialSize radial size of particle to create [GeV^-1]
    /// \return pointer to new InitialStateParticle object
    InitialStateParticle* createInitialStateParticle(int PDG, double radialSize);

    /// Create a Resonance from a PDG code and a MassShape
    /// \param PDG PDG code of particle to create
    /// \param radialSize Radial size of particle to create [GeV^-1]
    /// \param massShape Pointer to MassShape object describing resonance
    /// \return pointer to new Resonance object
    Resonance* createResonance(int PDG, double radialSize, MassShape* massShape);

    /// Create a Resonance from a PDG code. Use BreitWigner as MassShape
    /// \param PDG PDG code of particle to create
    /// \param radialSize Radial size of particle to create [GeV^-1]
    /// \return pointer to new Resonance object
    Resonance* createResonanceBreitWigner(int PDG, double radialSize);

    /// Create a Resonance from a PDG code. Use RelativisticBreitWigner as MassShape
    /// \param PDG PDG code of particle to create
    /// \param radialSize Radial size of particle to create [GeV^-1]
    /// \return pointer to new Resonance object
    Resonance* createResonanceRelativisticBreitWigner(int PDG, double radialSize);

    /// Create and add DecayChannel to parent
    /// \param parent Particle to add decay channel to
    /// \param daughterA first daughter particle to add to decay channel
    /// \param daughterB second daughter particle to add to decay channel
    /// \param L relative angular momentum between daughters in decay channel
    void createChannel(DecayingParticle* parent, Particle* daughterA, Particle* daughterB, unsigned L);

    /// Create quantum number object from PDG code
    /// \param PDG PDG code to look up
    /// \return Quantum numbers corresponding to particle
    QuantumNumbers createQuantumNumbers(int PDG);

private:

};

}

#endif

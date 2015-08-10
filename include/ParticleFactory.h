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

class ParticleFactory
{
public:
    /// Constructor
    ParticleFactory() {;}

    /// Create a FinalStateParticle from a PDG code
    FinalStateParticle* createFinalStateParticle(int PDG);

    /// Create an InitialStateParticle from a PDG code and a MassShape
    InitialStateParticle* createInitialStateParticle(int PDG, double radialSize);

    /// Create a FinalStateParticle from a PDG code and a MassShape
    Resonance* createResonance(int PDG, double radialSize, MassShape* massShape);

    /// Create a FinalStateParticle from a PDG code. Use BreitWigner as MassShape
    Resonance* createResonanceBreitWigner(int PDG, double radialSize);

    /// Create a FinalStateParticle from a PDG code. Use RelativisticBreitWigner as MassShape
    Resonance* createResonanceRelativisticBreitWigner(int PDG, double radialSize);

    /// Create and add DecaChannel to parent
    void createChannel(DecayingParticle* parent, Particle* daughterA, Particle* daughterB, unsigned L);

    QuantumNumbers createQuantumNumbers(int PDG);

private:

};

}

#endif

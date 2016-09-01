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

#ifndef yap_PhaseSpaceFactor_h
#define yap_PhaseSpaceFactor_h

#include "fwd/PhaseSpaceFactor.h"

#include "fwd/DataPoint.h"
#include "fwd/DecayChannel.h"
#include "fwd/MassShape.h"
#include "fwd/ParticleCombination.h"
#include "fwd/SpinAmplitude.h"

#include "AmplitudeComponent.h"

#include <complex>
#include <memory>

namespace yap {

/// Base class for calculating a phase-space factor
/// \defgroup PhaseSpaceFactor Classes related to phase-space factor calculation
/// \author Daniel Greenwald
class PhaseSpaceFactor : public RecalculableAmplitudeComponent
{
public:
    /// Constructor
    PhaseSpaceFactor(const ParticleCombinationEqualTo& equal)
        : RecalculableAmplitudeComponent(equal) {}
    
    /// grant friend status to DecayChannel to call addParticleCombination
    friend class DecayChannel;
};

/// Base class for factory creating PhaseSpaceFactor objects
/// \ingroup PhaseSpaceFactor
/// \author Daniel Greenwald
class PhaseSpaceFactorFactory
{
protected:

    /// (create and) return a PhaseSpaceFactor object applicaple for
    /// the final state of DecayChannel and the orbital angular
    /// momentum of SpinAmplitude
    /// \param dc DecayChannel
    /// \param sa SpinAmplitude
    /// \param ms MassShape
    virtual std::shared_ptr<PhaseSpaceFactor> phaseSpaceFactor(const DecayChannel& dc, const SpinAmplitude& sa, std::shared_ptr<MassShape> ms) = 0;

    /// grant friend status to DecayingParticle to call phaseSpaceFactor
    friend class DecayingParticle;

    /// grant friend status to Resonance to call phaseSpaceFactor
    friend class Resonance;
};

}

#endif

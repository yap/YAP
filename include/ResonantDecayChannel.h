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

#ifndef yap_ResonantDecayChannel_h
#define yap_ResonantDecayChannel_h

#include "DecayChannel.h"

#include <complex>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace yap {

/// \class ResonantDecayChannel
/// \brief Class implementing a decay channel of a resonance decaying to 2 daughter particles.
/// \author Johannes Rauch, Daniel Greenwald
class ResonantDecayChannel :
    public DecayChannel
{
public:

    /// \name Constructors
    /// @{

    /// N-particle Constructor (only valid for 2 particles).
    /// ResonantDecayChannel inherits ISP from daughters.
    /// \param daughters Vector of shared_ptr's to daughter Particle's
    ResonantDecayChannel(const ParticleVector& daughters);

    /// @}

    /// check consistency of object
    virtual bool consistent() const;

    /// \name Getters
    /// @{

    /// Get SpinAmplitude objects
    virtual const SpinAmplitudeVector& spinAmplitudes() const override
    { return SpinAmplitudes_; }

    /// @}

    /// add a spin amplitude
    void addSpinAmplitude(std::shared_ptr<SpinAmplitude> sa);

    /// Grant friend status to DecayingParticle to set itself as owner
    /// and to call fixSolitaryFreeAmplitudes()
    friend DecayingParticle;

protected:

    /// Add particle combination
    virtual void addParticleCombination(std::shared_ptr<ParticleCombination> c);

    /// set raw pointer to owning DecayingParticle
    virtual void setDecayingParticle(DecayingParticle* dp);

private:

    /// Vector of SpinAmplitudes
    SpinAmplitudeVector SpinAmplitudes_;

};


}

#endif

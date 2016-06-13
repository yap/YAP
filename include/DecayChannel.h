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

#ifndef yap_DecayChannel_h
#define yap_DecayChannel_h

#include "fwd/DecayChannel.h"

#include "fwd/CachedDataValue.h"
#include "fwd/DataPartition.h"
#include "fwd/DataPoint.h"
#include "fwd/DecayingParticle.h"
#include "fwd/FinalStateParticle.h"
#include "fwd/FreeAmplitude.h"
#include "fwd/Model.h"
#include "fwd/Parameter.h"
#include "fwd/Particle.h"
#include "fwd/SpinAmplitude.h"
#include "fwd/StatusManager.h"

#include "ParticleCombination.h"

#include <complex>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace yap {

/// \class DecayChannel
/// \brief Class implementing a decay channel.
/// \author Johannes Rauch, Daniel Greenwald
class DecayChannel
{
public:

    /// \name Constructors
    /// @{

    /// N-particle Constructor (at the moment only valid for 2 particles).
    /// DecayChannel inherits ISP from daughters.
    /// \param daughters Vector of shared_ptr's to daughter Particle's
    DecayChannel(const ParticleVector& daughters);

    /// @}

    /// check consistency of object
    virtual bool consistent() const;

    /// \return vector of shared_ptr's to final-state particles of channel (recursively checked)
    std::vector<std::shared_ptr<FinalStateParticle> > finalStateParticles() const;

    /// \name Getters
    /// @{

    /// Get Daughters
    const ParticleVector& daughters() const
    { return Daughters_; }

    /// Get SpinAmplitude objects
    const SpinAmplitudeVector& spinAmplitudes() const
    { return SpinAmplitudes_; }

    /// \return vector of ParticleCombinations
    const ParticleCombinationVector& particleCombinations() const
    { return ParticleCombinations_; }

    /// \return raw pointer to model through first Daughter
    const Model* model() const;

    /// \return raw pointer to owning DecayingParticle
    DecayingParticle* decayingParticle() const
    { return DecayingParticle_; }

    /// \return set of FreeAmplitude's for this decay channel
    FreeAmplitudeSet freeAmplitudes() const;

    /// @}

    /// add a spin amplitude
    void addSpinAmplitude(std::shared_ptr<SpinAmplitude> sa);

    /// Grant friend status to DecayingParticle to set itself as owner
    /// and to call fixSolitaryFreeAmplitudes()
    friend DecayingParticle;

protected:

    /// Add particle combination
    virtual void addParticleCombination(std::shared_ptr<ParticleCombination> c);

    /// prune ParticleCombinations_ to only contain ParticleCombination's tracing back up the ISP
    virtual void pruneParticleCombinations();

    /// call fixSolitaryFreeAmplitudes() on each daughter
    void fixSolitaryFreeAmplitudes();

    /// set raw pointer to owning DecayingParticle
    void setDecayingParticle(DecayingParticle* dp);

private:

    /// daughters of the decay
    ParticleVector Daughters_;

    /// Vector of SpinAmplitudes
    SpinAmplitudeVector SpinAmplitudes_;

    /// raw pointer owning DecayingParticle
    DecayingParticle* DecayingParticle_;

    /// vector of shared_ptr<ParticleCombination>
    ParticleCombinationVector ParticleCombinations_;

};

/// convert to string
std::string to_string(const DecayChannel& dc);

/// << operator
inline std::ostream& operator<< (std::ostream& os, const DecayChannel& dc)
{ os << to_string(dc); return os; }

}

#endif

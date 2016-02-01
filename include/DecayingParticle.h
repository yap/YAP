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

#ifndef yap_DecayingParticle_h
#define yap_DecayingParticle_h

#include "BlattWeisskopf.h"
#include "DataAccessor.h"
#include "DataPoint.h"
#include "DecayChannel.h"
#include "Exceptions.h"
#include "HelicitySpinAmplitude.h"
#include "make_unique.h"
#include "Particle.h"
#include "QuantumNumbers.h"
#include "SpinAmplitudeCache.h"

#include <complex>
#include <map>
#include <memory>
#include <vector>

namespace yap {

class FinalStateParticle;
class InitialStateParticle;
class ParticleCombination;

/// \class DecayingParticle
/// \brief Class for a particle that will decay
/// \authors Johannes Rauch, Daniel Greenwald
/// \ingroup Particle
///
/// The amplitude function returns a sum over the amplitudes of all
/// #DecayChannel's for the decay of the particle (denoted P, with
/// daughters in channel c denoted D1, D2; and amplitude A_c):\n
/// A_c = a_c * Blatt-Weisskopf(P->D1+D2) * SpinAmplitude(P->D1+D2) * A(D1->xx) * A(D2->xx)\n
/// with free amplitude a_c.

class DecayingParticle : public Particle, public DataAccessor
{
public:

    /// Constructor
    DecayingParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize);

    /// Calculate complex amplitude
    /// \param d DataPoint to calculate with
    /// \param pc (shared_ptr to) ParticleCombination to calculate for
    /// \param two_m 2 * the spin projection to calculate for
    /// \param dataPartitionIndex partition index for parallelization
    virtual std::complex<double> amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc,
                                           int two_m, unsigned dataPartitionIndex) const override;

    /// Check consistency of object
    virtual bool consistent() const override;

    /// Add a DecayChannel and set its parent to this DecayingParticle.
    /// \param c unique_ptr to DecayChannel, should be constructed in function call, or use std::move(c)
    virtual void addChannel(std::unique_ptr<DecayChannel> c);

    /// Add a DecayChannel and set its parent to this DecayingParticle.
    /// \param daughters ParticleVector of daughters to create DecayChannel object from
    virtual void addChannel(const ParticleVector& daughters)
    { addChannel(std::make_unique<DecayChannel>(daughters)); }


    /// Return final state particles of a channel (vector should be identical for all channels)
    /// \return vector of shared_ptr's to FinalStateParticles of this decaying particle (in channel i)
    /// \param i index of DecayChannel to return FinalStateParticles of.
    std::vector< std::shared_ptr<FinalStateParticle> > finalStateParticles(unsigned i = 0) const;

    /// \name Getters
    /// @{

    /// return channels
    const DecayChannelVector& channels() const
    { return Channels_;}

    /// \return Number of decay channels for this object
    unsigned int nChannels() const
    { return Channels_.size(); }

    /// Return Channel i
    const DecayChannel* channel(unsigned i) const
    { return Channels_.at(i).get(); }

    /// Return Channel i
    DecayChannel* channel(unsigned i)
    { return Channels_.at(i).get(); }

    /// \return Radial size [GeV^-1]
    std::shared_ptr<RealParameter> radialSize()
    { return RadialSize_; }

    /// @}

    /// Print complete decay chain
    void printDecayChain() const
    { printDecayChainLevel(0); }

    // virtual ParameterSet ParametersItDependsOn() override;

    virtual CachedDataValueSet CachedDataValuesItDependsOn() override;

    /// \return raw pointer to initial state particle through first DecayChannel
    InitialStateParticle* initialStateParticle() override
    { return Channels_.empty() ? nullptr : Channels_[0]->initialStateParticle(); }

    /// grant friend status to DecayChannel to see BlattWeiskopffs_, and get dataAccessors
    friend DecayChannel;

protected:

    /// add ParticleCombination to SymmetrizationIndices_ and BlattWeisskopfs_
    virtual void addParticleCombination(std::shared_ptr<ParticleCombination> c) override;

    void printDecayChainLevel(int level) const;

    /// \return set of shared_ptr's of DataAccessor's
    virtual DataAccessorSet dataAccessors() override;

    /// \return vector of shared_ptr's to all free amplitudes from this point in decay tree and down
    virtual ComplexParameterVector freeAmplitudes() const;

private:

    /// vector of decay channel objects
    DecayChannelVector Channels_;

    /// map of Blatt-Weisskopf barrier factors, key = angular momentum
    std::map<unsigned, std::shared_ptr<BlattWeisskopf> > BlattWeisskopfs_;

    /// Radial size parameter [GeV^-1]
    std::shared_ptr<RealParameter> RadialSize_;

    /// Cached amplitudes for each spin projection
    /// key = 2 * spin projection
    std::map<int, std::shared_ptr<ComplexCachedDataValue> > Amplitudes_;

};

}

#endif

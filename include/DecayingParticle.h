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

#include "fwd/BlattWeisskopf.h"
#include "fwd/DecayChannel.h"
#include "fwd/DecayTree.h"
#include "fwd/FinalStateParticle.h"
#include "fwd/FreeAmplitude.h"
#include "fwd/Model.h"
#include "fwd/ParticleCombination.h"
#include "fwd/QuantumNumbers.h"

#include "Particle.h"

#include <memory>

namespace yap {

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

class DecayingParticle : public Particle
{
protected:

    /// Constructor
    /// see #create
    DecayingParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize);

public:

    /// create
    /// \param q QuantumNumbers of decaying particle
    /// \param mass mass of decaying particle
    /// \param radialSize radial size of decaying particle
    static std::shared_ptr<DecayingParticle> create(const QuantumNumbers& q, double mass, std::string name, double radialSize)
    { return std::shared_ptr<DecayingParticle>(new DecayingParticle(q, mass, name, radialSize)); }

    /// \return DecayTrees
    /// map key is spin projection
    const DecayTreeVectorMap& decayTrees() const
    { return DecayTrees_; }

    /// Check consistency of object
    virtual bool consistent() const override;

    /// Check if a DecayChannel is valid for DecayingParticle;
    /// will throw if invalid
    virtual void checkDecayChannel(const DecayChannel& c) const;

    /// Add a DecayChannel and set its parent to this DecayingParticle.
    /// \param c unique_ptr to DecayChannel, should be constructed in function call, or use std::move(c)
    /// \return shared_ptr to DecayChannel that has been added
    virtual std::shared_ptr<DecayChannel> addChannel(std::shared_ptr<DecayChannel> c);

    /// Add a DecayChannel and set its parent to this DecayingParticle.
    /// \param daughters ParticleVector of daughters to create DecayChannel object from
    /// \return shared_ptr to DecayChannel that has been added
    std::shared_ptr<DecayChannel> addChannel(const ParticleVector& daughters);

    /// Add a DecayChannel and set its parent to this DecayingParticle.
    /// \param daughters... shared_ptr's to daughters to create DecayChannel object from
    /// \return shared_ptr to DecayChannel that has been added
    template <typename ... Types>
    std::shared_ptr<DecayChannel> addChannel(Types ... daughters)
    { ParticleVector V{daughters...}; return addChannel(V); }
    /* { ParticleVector V; fill_vector(daughters..., V); return addChannel(V); } */

    /// Return final state particles of a channel
    /// \return vector of shared_ptr's to FinalStateParticles of this decaying particle (in channel i)
    /// \param i index of DecayChannel to return FinalStateParticles of.
    FinalStateParticleVector finalStateParticles(unsigned i) const;

    /// \name Getters
    /// @{

    /// \return channel
    std::shared_ptr<DecayChannel> channel(const ParticleVector& daughters);

    /// \return channels
    const DecayChannelVector& channels() const
    { return Channels_;}

    /// \return Radial size [GeV^-1]
    std::shared_ptr<RealParameter> radialSize()
    { return RadialSize_; }

    /// \return Blatt-Weisskopf factors
    const BlattWeisskopfMap& blattWeisskopfs() const
    { return BlattWeisskopfs_; }

    FreeAmplitudeSet freeAmplitudes() const;

    /// @}

    /// Print complete decay chain
    void printDecayChain() const
    { printDecayChainLevel(0); }

    /// \return raw pointer to Model through first DecayChannel
    const Model* model() const override;

    /// grant friend status to DecayChannel to call fixSolitaryFreeAmplitudes()
    /// and storeBlattWeisskopf()
    friend DecayChannel;

    /// grant friend status to Model to call fixSolitaryFreeAmplitudes()
    friend Model;

protected:

    /// register any necessary DataAccessor's with model
    virtual void registerWithModel()
    {}

    /// tell DecayingParticle to store a BlattWeisskopf object for L
    /// \param l orbital angular momentum of breakup
    void storeBlattWeisskopf(unsigned l);

    /// add ParticleCombination to SymmetrizationIndices_ and BlattWeisskopfs_
    virtual void addParticleCombination(const std::shared_ptr<ParticleCombination>& c) override;

    /// prune ParticleCombinations_ to only contain ParticleCombination's tracing back up the ISP
    virtual void pruneParticleCombinations() override;

    /// if only one decay channel is available, fix its free amplitude to the current value
    void fixSolitaryFreeAmplitudes();

    void printDecayChainLevel(int level) const;

    /// modify a DecayTree
    /// \param dt DecayTree to modify
    virtual void modifyDecayTree(DecayTree& dt) const;

private:

    /// vector of decay channel objects
    DecayChannelVector Channels_;

    /// map of Blatt-Weisskopf barrier factors, key = angular momentum
    BlattWeisskopfMap BlattWeisskopfs_;

    /// Radial size parameter [GeV^-1]
    std::shared_ptr<RealParameter> RadialSize_;

    /// Map of spin projection to DecayTreeVector
    DecayTreeVectorMap DecayTrees_;

};

/// convert to (multiline) string
std::string to_string(const DecayTreeVectorMap& m_dtv_map);

/// \return set of free amplitudes in map of spin projection to decay tree vector
FreeAmplitudeSet freeAmplitudes(const DecayTreeVectorMap& m_dtv_map);

}

#endif

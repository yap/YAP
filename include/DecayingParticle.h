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

#include "AmplitudeComponent.h"
#include "DataAccessor.h"
#include "DecayChannel.h"
#include "Particle.h"

#include <complex>
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

class DecayingParticle : public Particle, public DataAccessor
{
public:

    /// \name Constructors
    /// @{

    /// Constructor
    DecayingParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize);

    /// Copy constructor
    DecayingParticle(const DecayingParticle& other) = delete;

    /// @}

    virtual std::complex<double> amplitude(DataPartition& d, std::shared_ptr<const ParticleCombination> pc) const override;

    /// Check consistency of object
    virtual bool consistent() const override;

    /// Add a DecayChannel and set its parent to this DecayingParticle.
    /// \param c DecayingParticle takes ownership of c, i.e. c will point to nullptr afterwards
    virtual void addChannel(std::unique_ptr<DecayChannel>& c);

    /// Add all possible two-body DecayChannels with #HelicitySpinAmplitudes up to a maximum relative angular momentum
    /// \param A daughter particle
    /// \param B daughter particle
    /// \param L maximum relative angular momentum between A and B * 2
    //virtual void addChannels(std::shared_ptr<Particle> A, std::shared_ptr<Particle> B, unsigned maxTwoL)
    //{ addChannels({A}, {B}, maxTwoL); }

    /// Add all possible two-body DecayChannels with #HelicitySpinAmplitudes up to a maximum relative angular momentum
    /// \param A daughter particle in all possible helicity states
    /// \param B daughter particle in all possible helicity states
    /// \param L maximum relative angular momentum between A and B * 2
    virtual void addChannels(std::shared_ptr<Particle> A, std::shared_ptr<Particle> B, unsigned maxTwoL);

    /// Return final state particles of a channel (vector should be identical for all channels)
    /// \return vector of shared_ptr's to FinalStateParticles of this decaying particle (in channel i)
    /// \param i index of DecayChannel to return FinalStateParticles of.
    std::vector< std::shared_ptr<FinalStateParticle> > finalStateParticles(unsigned i = 0) const;

    /// \name Getters
    /// @{

    /// return channels
    const std::vector< std::unique_ptr<yap::DecayChannel> >& channels() const
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

    /// \name Setters
    /// @{

    /// Set pointer to initial state particle
    void setInitialStateParticle(InitialStateParticle* isp) override;

    /// @}

    /// SpinAmplitudes can be shared among DecayChannels if the QuantumNumbers are equal.
    /// Check if this is the case, and share SpinAmplitudes
    void optimizeSpinAmplitudeSharing();

    /// Print complete decay chain
    void printDecayChain() const
    { printDecayChainLevel(0); }

    /// Print SpinAmplitudes involved in decay chain
    void printSpinAmplitudes(int level = 0);

    // for internal use only
    virtual void setSymmetrizationIndexParents() override;

private:

    void printDecayChainLevel(int level) const;

    /// vector of decay channel objects
    std::vector< std::unique_ptr<DecayChannel> > Channels_;

    /// Radial size parameter [GeV^-1]
    std::shared_ptr<RealParameter> RadialSize_;
};

}

#endif

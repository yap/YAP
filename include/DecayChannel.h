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

#include "AmplitudeComponent.h"
#include "Constants.h"
#include "DataAccessor.h"
#include "DataPoint.h"
#include "Particle.h"
#include "QuantumNumbers.h"
#include "SpinAmplitude.h"

#include <complex>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace yap {

class DecayingParticle;
class FinalStateParticle;
class InitialStateParticle;
class ParticleCombination;

/// \class InitialStateParticle
/// \brief Class implementing a decay channel.
/// \author Johannes Rauch, Daniel Greenwald
class DecayChannel : public AmplitudeComponent, public DataAccessor
{
public:

    /// \class AmplitudePair
    struct AmplitudePair
    {
        AmplitudePair(DecayChannel* dc, std::complex<double> free = Complex_1);
        std::shared_ptr<ComplexCachedDataValue> Fixed;
        std::shared_ptr<ComplexParameter> Free;
    };
    
    /// \name Constructors
    /// @{

    /// N-particle Constructor (at the moment only valid for 2 particles).
    /// DecayChannel inherits ISP from daughters.
    /// \param daughters Vector of shared_ptr's to daughter Particle's
    DecayChannel(ParticleVector daughters);

    /// @}

    /// Calculate complex amplitude
    virtual std::complex<double> amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc, unsigned dataPartitionIndex) const override;

    /// check consistency of object
    virtual bool consistent() const override;

    /// \return vector of shared_ptr's to final-state particles of channel (recursively checked)
    std::vector<std::shared_ptr<FinalStateParticle> > finalStateParticles() const;

    /// \name Getters
    /// @{

    /// Get Daughters
    const ParticleVector& daughters() const
    { return Daughters_; }

    /// Get SpinAmplitude objects
    SpinAmplitudeVector spinAmplitudes();

    /// Get SpinAmplitude objects (const)
    const SpinAmplitudeVector spinAmplitudes() const
    { return const_cast<DecayChannel*>(this)->spinAmplitudes(); }

    /// Get AmplitudePair object corresponding to SpinAmplitude
    AmplitudePair& amplitudes(const std::shared_ptr<SpinAmplitude>& sa);

    /// Get AmplitudePair object corresponding to SpinAmplitude (const)
    const AmplitudePair& amplitudes(const std::shared_ptr<SpinAmplitude>& sa) const
    { return const_cast<DecayChannel*>(this)->amplitudes(sa); }
        
    std::shared_ptr<ComplexParameter> freeAmplitude() const
    { return FreeAmplitude_; }

    /// @}

    virtual ParameterSet ParametersItDependsOn() override
    { return {FreeAmplitude_}; }

    virtual CachedDataValueSet CachedDataValuesItDependsOn() override
    { return {FixedAmplitude_}; }

    /// \return raw pointer to initial state particle through first Daughter
    InitialStateParticle* initialStateParticle() override
    { return Daughters_[0]->initialStateParticle(); }

    /// \return raw pointer to owning DecayingParticle
    DecayingParticle* decayingParticle() const
    { return DecayingParticle_; }

    /// Grant friend status to DecayingParticle to set itself as owner
    friend DecayingParticle;

    /// Grant friend status to InitialStateParticle
    friend class InitialStateParticle;

protected:

    /// set raw pointer to owning DecayingParticle
    void setDecayingParticle(DecayingParticle* dp);

    /// clear SymmetrizationIndices_
    virtual void clearSymmetrizationIndices() override;

    // sets symmetrization index parents
    void setSymmetrizationIndexParents();

    /// \return set of DataAccessors
    virtual DataAccessorSet dataAccessors();

private:

    /// daughters of the decay
    ParticleVector Daughters_;

    /// Map of SpinAmplitude's (key; by shared_ptr)  to vector of AmplitudePair's (value),
    /// with one AmplitudePair per spin projection of parent in SpinAmplitude
    SpinAmplitudeMap<std::vector<AmplitudePair> > Amplitudes_;
    
    /// SpinAmplitude can be shared between several DecayChannels
    std::shared_ptr<SpinAmplitude> SpinAmplitude_;

    /// Independent amplitude for channel's contribution to full model
    std::shared_ptr<ComplexParameter> FreeAmplitude_;

    /// amplitude from spin dynamics, mass shapes, etc.
    std::shared_ptr<ComplexCachedDataValue> FixedAmplitude_;

    /// raw pointer owning DecayingParticle
    DecayingParticle* DecayingParticle_;

};

/// convert to string
std::string to_string(const DecayChannel& dc);

/// << operator
inline std::ostream& operator<< (std::ostream& os, const DecayChannel& dc)
{ os << to_string(dc); return os; }

/// \typedef DecayChannelVector
using DecayChannelVector = std::vector<std::shared_ptr<DecayChannel> >;

}

#endif

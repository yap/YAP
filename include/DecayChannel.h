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
#include "BlattWeisskopf.h"
#include "DataAccessor.h"
#include "DataPoint.h"
#include "Particle.h"
#include "QuantumNumbers.h"

#include <complex>
#include <memory>
#include <string>
#include <vector>

namespace yap {

class DecayingParticle;
class FinalStateParticle;
class InitialStateParticle;
class ParticleCombination;
class SpinAmplitude;

/// \class InitialStateParticle
/// \brief Class implementing a decay channel.
/// \author Johannes Rauch, Daniel Greenwald
class DecayChannel : public AmplitudeComponent, public DataAccessor
{
public:

    /// \name Constructors
    /// @{

    /// N-particle Constructor [at the moment only valid for 2 particles].
    /// DecayChannel inherits ISP from daughters.
    /// \param daughters Vector of shared_ptr's to daughter Particle's
    /// \param spinAmplitude shared_ptr to SpinAmplitude object
    DecayChannel(ParticleVector daughters, std::shared_ptr<SpinAmplitude> spinAmplitude);

    /// 2-particle Constructor
    DecayChannel(std::shared_ptr<Particle> A, std::shared_ptr<Particle> B, std::shared_ptr<SpinAmplitude> spinAmplitude)
        : DecayChannel( {A, B}, spinAmplitude) {}

    /// @}

    /// Calculate complex amplitude
    virtual std::complex<double> amplitude(DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc, unsigned dataPartitionIndex) const override;

    /// check consistency of object
    virtual bool consistent() const override;

    /// \return vector of shared_ptr's to final-state particles of channel (recursively checked)
    std::vector<std::shared_ptr<FinalStateParticle> > finalStateParticles() const;

    /// \name Getters
    /// @{

    /// Get Daughters
    const ParticleVector& daughters() const
    { return Daughters_; }

    /// Get SpinAmplitude object
    std::shared_ptr<SpinAmplitude>& spinAmplitude()
    { return SpinAmplitude_; }

    /// Get SpinAmplitude object (const)
    const SpinAmplitude* spinAmplitude() const
    { return SpinAmplitude_.get(); }

    std::shared_ptr<ComplexParameter> freeAmplitude() const
    { return FreeAmplitude_; }

    /// @}

    virtual ParameterSet ParametersItDependsOn() override
    { return {FreeAmplitude_}; }

    virtual CachedDataValueSet CachedDataValuesItDependsOn() override
    { return {FixedAmplitude_}; }

    /// \return raw pointer to initial state particle through first daughter
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

    /// add symmetrizationIndex to SymmetrizationIndices_,
    /// also add to BlattWeisskopf_ and SpinAmplitude_
    virtual void addSymmetrizationIndex(std::shared_ptr<const ParticleCombination> c) override;

    // sets symmetrization index parents
    void setSymmetrizationIndexParents();

    /// \return set of DataAccessors
    virtual DataAccessorSet dataAccessors();

private:

    /// daughters of the decay
    ParticleVector Daughters_;

    /// Blatt-Weisskopf calculator
    std::shared_ptr<BlattWeisskopf> BlattWeisskopf_;

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

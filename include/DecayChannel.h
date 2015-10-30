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

#ifndef yap_DecayChannel_h
#define yap_DecayChannel_h

#include "AmplitudeComponent.h"
#include "BlattWeisskopf.h"
#include "DataAccessor.h"
#include "DataPoint.h"

#include <complex>
#include <memory>
#include <string>
#include <vector>

namespace yap {

class DecayingParticle;
class FinalStateParticle;
class InitialStateParticle;
class Particle;
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

    /// N-particle Constructor [at the moment only valid for 2 particles]
    DecayChannel(std::vector<std::shared_ptr<Particle> > daughters, std::shared_ptr<SpinAmplitude> spinAmplitude, DecayingParticle* parent);

    /// 2-particle Constructor
    DecayChannel(std::shared_ptr<Particle> daughterA, std::shared_ptr<Particle> daughterB, std::shared_ptr<SpinAmplitude> spinAmplitude, DecayingParticle* parent);

    /// @}

    /// Calculate complex amplitude
    virtual std::complex<double> amplitude(DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc, unsigned dataPartitionIndex) const override;

    /// \return #CalculationStatus of symmetrization index and data-partition index
    /// \param pc shared pointer to #ParticleCombination to check status of
    /// \param dataPartitionIndex index of dataPartitionIndex to check status of
    virtual CalculationStatus calculationStatus(const std::shared_ptr<const ParticleCombination>& pc, unsigned symmetrizationIndex, unsigned dataPartitionIndex) const override;

    /// check consistency of object
    virtual bool consistent() const override;

    /// cast into string
    operator std::string() const;

    /// \return vector of shared_ptr's to final-state particles of channel (recursively checked)
    std::vector<std::shared_ptr<FinalStateParticle> > finalStateParticles() const;

    /// \name Getters
    /// @{

    /// Get Daughters
    const std::vector<std::shared_ptr<Particle> >& daughters() const
    { return Daughters_; }

    /// Get SpinAmplitude object
    std::shared_ptr<SpinAmplitude>& spinAmplitude()
    { return SpinAmplitude_; }

    /// Get SpinAmplitude object (const)
    const SpinAmplitude* spinAmplitude() const
    { return SpinAmplitude_.get(); }

    /// Get parent particle
    DecayingParticle* parent()
    { return Parent_; }

    /// Get parent particle
    const DecayingParticle* parent() const
    { return Parent_; }

    std::shared_ptr<ComplexParameter> freeAmplitude() const
    { return FreeAmplitude_; }

    /// @}

    /// \name Setters
    /// @{

    /// Set pointer to initial state particle
    void setInitialStateParticle(InitialStateParticle* isp) override;

    /// @}

    /// \name SymmetrizationIndex related
    /// @{

    /// add symmetrizationIndex to SymmetrizationIndices_,
    /// also add to BlattWeisskopf_ and SpinAmplitude_
    virtual void addSymmetrizationIndex(std::shared_ptr<const ParticleCombination> c) override;

    /// clear SymmetrizationIndices_
    virtual void clearSymmetrizationIndices() override;

    // for internal use only
    void setSymmetrizationIndexParents();

    /// @}

    // for internal use only
    void addSpinAmplitudeDependencies();

    virtual ParameterSet ParametersItDependsOn() override
    { return {FreeAmplitude_}; }

    virtual CachedDataValueSet CachedDataValuesItDependsOn() override
    { return {FixedAmplitude_}; }

protected:

    /// DecayingParticle this DecayChannel belongs to
    DecayingParticle* Parent_;

    /// 2 daughters of the decay
    std::vector<std::shared_ptr<Particle> > Daughters_;

    /// Blatt-Weisskopf calculator
    std::unique_ptr<BlattWeisskopf> BlattWeisskopf_;

    /// SpinAmplitude can be shared between several DecayChannels
    std::shared_ptr<SpinAmplitude> SpinAmplitude_;

    std::shared_ptr<ComplexParameter> FreeAmplitude_;
    std::shared_ptr<ComplexCachedDataValue> FixedAmplitude_;

};

}

#endif

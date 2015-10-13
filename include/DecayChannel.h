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

#include "AmplitudeComponentDataAccessor.h"
#include "BlattWeisskopf.h"
#include "ParameterSet.h"

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

class DecayChannel : public AmplitudeComponentDataAccessor, public ParameterSet
{
public:
    /// N-particle Constructor [at the moment only valid for 2 particles]
    DecayChannel(std::vector<std::shared_ptr<Particle> > daughters, std::shared_ptr<SpinAmplitude> spinAmplitude, unsigned nDataPartitions = 1);

    /// 2-particle Constructor
    DecayChannel(std::shared_ptr<Particle> daughterA, std::shared_ptr<Particle> daughterB, std::shared_ptr<SpinAmplitude> spinAmplitude, unsigned nDataPartitions = 1);

    /// check consistency of object
    virtual bool consistent() const override;

    virtual CalculationStatus updateCalculationStatus(DataPartition& d, std::shared_ptr<const ParticleCombination> c) const override;

    /// \return breakup momentum
    double breakupMomentum() const;

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

    std::complex<double> freeAmplitude() const
    { return this->at(0)->value(); }

    /// @}

    /// \name Setters
    /// @{

    /// Set pointer to initial state particle
    void setInitialStateParticle(InitialStateParticle* isp) override;

    /// Set parent
    void setParent(DecayingParticle* parent)
    { Parent_ = parent; }

    void setFreeAmplitude(std::complex<double> a)
    { at(0)->setValue(a); }

    /// @}

    /// add symmetrizationIndex to SymmetrizationIndices_,
    /// also add to BlattWeisskopf_ and SpinAmplitude_
    virtual void addSymmetrizationIndex(std::shared_ptr<const ParticleCombination> c) override;

    /// clear SymmetrizationIndices_
    virtual void clearSymmetrizationIndices() override;

    // for internal use only
    void setSymmetrizationIndexParents();

    virtual void precalculate() override;

protected:

    virtual void calcPrecalculate() override
    {}

    /// \return (fixed) Amplitude for decay channel
    virtual std::complex<double> calcAmplitude(DataPartition& d, std::shared_ptr<const ParticleCombination> pc) const override;

    /// 2 daughters of the decay
    std::vector<std::shared_ptr<Particle> > Daughters_;

    /// Blatt-Weisskopf calculator
    BlattWeisskopf BlattWeisskopf_;

    /// SpinAmplitude can be shared between several DecayChannels
    std::shared_ptr<SpinAmplitude> SpinAmplitude_;

    /// DecayingParticle this DecayChannel belongs to
    DecayingParticle* Parent_;

};

}

#endif

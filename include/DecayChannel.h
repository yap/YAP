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

#include <memory>
#include <string>
#include <vector>

namespace yap {

class DecayingParticle;
class FinalStateParticle;
class Particle;
class ParticleCombination;
class SpinAmplitude;

/// \class InitialStateParticle
/// \brief Class implementing a decay channel.
/// \author Johannes Rauch, Daniel Greenwald

class DecayChannel : public AmplitudeComponent, public DataAccessor
{
public:
    /// N-particle Constructor [at the moment only valid for 2 particles]
    DecayChannel(std::vector<std::shared_ptr<Particle> > daughters, std::shared_ptr<SpinAmplitude> spinAmplitude);

    /// 2-particle Constructor
    DecayChannel(std::shared_ptr<Particle> daughterA, std::shared_ptr<Particle> daughterB, std::shared_ptr<SpinAmplitude> spinAmplitude);

    /// \return Amplitude for decay channel
    virtual Amp amplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc) override;

    /// check consistency of object
    virtual bool consistent() const override;

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

    /// Get free amplitude
    Amp freeAmplitude() const
    { return FreeAmplitude_; }

    /// Get parent particle
    DecayingParticle* parent()
    { return Parent_; }

    /// Get parent particle
    const DecayingParticle* parent() const
    { return Parent_; }

    /// @}

    /// \name Setters
    /// @{

    /// Set the free (complex) amplitude
    void setFreeAmplitude(const Amp& amp)
    { FreeAmplitude_ = amp; }

    /// Set parent
    void setParent(DecayingParticle* parent)
    { Parent_ = parent; }

    /// @}

    /// add symmetrizationIndex to SymmetrizationIndices_,
    /// also add to BlattWeisskopf_ and SpinAmplitude_
    virtual void addSymmetrizationIndex(std::shared_ptr<ParticleCombination> c) override;

    /// clear SymmetrizationIndices_
    virtual void clearSymmetrizationIndices() override;

    // for internal use only
    void setSymmetrizationIndexParents();

protected:

    /// 2 daughters of the decay
    std::vector<std::shared_ptr<Particle> > Daughters_;

    /// Blatt-Weisskopf calculator
    BlattWeisskopf BlattWeisskopf_;

    /// SpinAmplitude can be shared between several DecayChannels
    std::shared_ptr<SpinAmplitude> SpinAmplitude_;

    /// free ("fit") amplitude to multiply all others by
    Amp FreeAmplitude_;

    /// DecayingParticle this DecayChannel belongs to
    DecayingParticle* Parent_;

};

}

#endif

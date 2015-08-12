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

#include "DataAccessor.h"

#include "BlattWeisskopf.h"
#include "Particle.h"
#include "SpinAmplitude.h"

#include <memory>

namespace yap {

class DecayingParticle;
typedef std::array<Particle*, 2> Daughters;

/// \class InitialStateParticle
/// \brief Class implementing a decay channel.
/// \author Johannes Rauch, Daniel Greenwald

class DecayChannel : public AmplitudeComponent, DataAccessor
{
public:
    /// Constructor
    DecayChannel(Particle* daughterA, Particle* daughterB, std::shared_ptr<SpinAmplitude> spinAmplitude);

    /// \return Amplitude for decay channel
    virtual Amp amplitude(DataPoint& d) override;

    /// check consistency of object
    virtual bool consistent() const override;

    /// \name Getters
    /// @{

    /// Get Daughters
    const Daughters& daughters() const {return Daughters_;}

    /// Get SpinAmplitude pointer
    const SpinAmplitude* spinAmplitude() const
    { return SpinAmplitude_.get(); }

    /// Get shared SpinAmplitude object
    std::shared_ptr<SpinAmplitude>& sharedSpinAmplitude()
    { return SpinAmplitude_; }

    /// Get free amplitude
    Amp freeAmplitude() const {return FreeAmplitude_;}

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

private:

    /// 2 daughters of the decay
    Daughters Daughters_;

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

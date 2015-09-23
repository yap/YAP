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

#ifndef yap_BlattWeisskopf_h
#define yap_BlattWeisskopf_h

#include "Amp.h"
#include "AmplitudeComponent.h"
#include "CalculationStatus.h"

#include <memory>

namespace yap {

class DecayChannel;
class ParticleCombination;

/// \class BlattWeisskopf
/// \brief Class implementing BlattWeisskopf barrier factors
/// \author Johannes Rauch, Daniel Greenwald

class BlattWeisskopf : public AmplitudeComponent
{
public:

    /// Constructor
    BlattWeisskopf(DecayChannel* decayChannel);

    /// check consistency of object
    virtual bool consistent() const override;

    /// Return DecayChannel this BlattWeisskopf belongs to
    DecayChannel* decayChannel() const {return DecayChannel_;}

    /// Blatt-Weisskopf amplitude at DataPoint
    virtual const Amp& amplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc) override;

private:

    /// DecayChannel this BlattWeisskopf belongs to
    DecayChannel* DecayChannel_;

    Amp Amplitude_;
    CalculationStatus CalculationStatus_;

};

}

#endif

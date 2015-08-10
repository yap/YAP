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

#ifndef yap_BlattWeisskopf_h
#define yap_BlattWeisskopf_h

#include "AmplitudeComponent.h"
#include "DataAccessor.h"

namespace yap {

class DecayChannel;

class BlattWeisskopf : public AmplitudeComponent, public DataAccessor
{
public:
    BlattWeisskopf(DecayChannel* decayChannel) :
        AmplitudeComponent(), DecayChannel_(decayChannel) {;}

    virtual Amp amplitude(DataPoint& d) override;
    virtual bool consistent() const override;

    DecayChannel* decayChannel() const {return DecayChannel_;}

private:
    DecayChannel* DecayChannel_; /// DecayChannel this BlattWeisskopf belongs to
};

}

#endif

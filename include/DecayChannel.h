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
#include <tuple>

namespace yap {

class Resonance;
typedef std::array<Particle*, 2> Daughters;

class DecayChannel : public DataAccessor
{
public:
    DecayChannel(Particle* daughterA, Particle* daughterB, unsigned int L, SpinAmplitude& spinAmplitude);
    ~DecayChannel() {;}

    virtual Amp amplitude(DataPoint& d);
    virtual bool consistent() const;

    const Daughters& daughters() const {return Daughters_;}
    const Particle* daughter(unsigned int i) const {return Daughters_.at(i);}
    const Particle* daughterA() const {return Daughters_[0];}
    const Particle* daughterB() const {return Daughters_[1];}

    unsigned char l() const {return L_;}
    const SpinAmplitude& spinAmplitude() const {return SpinAmplitude_;}
    Amp freeAmplitude() const {return FreeAmplitude_;}
    Resonance* resonance() const {return Resonance_;}

    void setFreeAmplitude(const Amp& amp) {FreeAmplitude_ = amp;}

private:
    Daughters Daughters_;
    unsigned char L_; /// relative angular momentum between daughters
    BlattWeisskopf BlattWeisskopf_;
    SpinAmplitude& SpinAmplitude_; /// SpinAmplitude can be shared between several DecayChannels
    Amp FreeAmplitude_;
    Resonance* Resonance_; /// Resonance this DecayChannel belongs to
};

}

#endif

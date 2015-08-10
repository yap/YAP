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

#ifndef yap_Resonance_h
#define yap_Resonance_h

#include "Particle.h"
#include "DecayChannel.h"
#include "MassShape.h"
#include <vector>

namespace yap {

class FinalStateParticle;

/// \ingroup Particle

class Resonance : public Particle
{
public:
    Resonance(const QuantumNumbers& q, double mass, std::string name, const MassShape& massShape, double radialSize) :
        Particle(q, mass, name), MassShape_(massShape), RadialSize_(radialSize) {;}
    //virtual ~Resonance() {;}

    virtual Amp amplitude(DataPoint& d);
    virtual bool consistent() const;

    const std::vector<const FinalStateParticle*> finalStateParticles(unsigned int channel = 0) const;
    unsigned int nChannels() const {return Channels_.size();}

    void addChannel(const DecayChannel& c) {Channels_.push_back(c);}

private:
    MassShape MassShape_;
    std::vector<yap::DecayChannel> Channels_;
    double RadialSize_;
};

}

#endif

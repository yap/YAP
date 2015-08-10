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

#ifndef yap_FinalStateParticle_h
#define yap_FinalStateParticle_h

#include "Particle.h"

namespace yap {

/// \ingroup Particle

class FinalStateParticle : public Particle
{
public:
    FinalStateParticle(const QuantumNumbers& q, double mass, std::string name, int pdg) :
        Particle(q, mass, name), PDGCode_(pdg) {;}

    virtual Amp amplitude(DataPoint& d) override {return Amp(1.);}
    //virtual bool consistent() const override {return Particle::consistent();}

    int pdgCode() const {return PDGCode_;}

private:
    int PDGCode_; /// PDG code of the particle
};

}

#endif

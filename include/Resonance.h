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

#include "DecayChannel.h"
#include "FinalStateParticle.h"
#include "MassShape.h"
#include <vector>

namespace yap {

/// \ingroup Particle

class Resonance : public Particle {
public:
  Resonance();
  ~Resonance();

  virtual Amp amplitude(DataPoint& d);
  virtual bool checkConsistency() const;

  const std::vector<const FinalStateParticle*> getFinalStateParticles(unsigned int channel = 0) const;

private:
  MassShape MassShape_;
  std::vector<DecayChannel*> Channels_; /// DecayChannel can be shared between several Resonances
  double RadialSize_;
};

}

#endif

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

#include "DecayingParticle.h"
#include "DataAccessor.h"
#include "MassShape.h"

#include <memory>

namespace yap {

class FinalStateParticle;

/// \class Resonance
/// \brief Class for a particle that will decay and has a mass shape
/// \authors Johannes Rauch, Daniel Greenwald
/// \ingroup Particle

class Resonance : public DecayingParticle, public DataAccessor
{
public:
    /// Constructor
    Resonance(const QuantumNumbers& q, double mass, std::string name, double radialSize, MassShape* massShape);

    virtual Amp amplitude(DataPoint& d) override;
    virtual bool consistent() const override;

    /// Return MassShape
    const MassShape* massShape() const {return MassShape_.get();}

private:
    std::unique_ptr<MassShape> MassShape_;

};

}

#endif

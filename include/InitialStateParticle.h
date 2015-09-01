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

#ifndef yap_InitialStateParticle_h
#define yap_InitialStateParticle_h

#include "DecayingParticle.h"

namespace yap {

/// \class InitialStateParticle
/// \brief Class implementing an initial state particle.
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Particle

class InitialStateParticle : public DecayingParticle
{
public:

    /// Constructor
    InitialStateParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize);

    //virtual Amp amplitude(DataPoint& d) override;

    /// Check consistency of object
    virtual bool consistent() const override;

    /// Check consistency of data point with this initial state and its final state particles
    /// \param d DataPoint to check
    /// \todo Flesh out
    virtual bool consisent(const DataPoint& d) const
    { return true; }

    /// Calculate helicity angles for all possible symmetrization indices
    void calculateHelicityAngles(DataPoint& d) const;

private:

    // vector of final state particles
    std::vector<FinalStateParticle*> FinalStateParticles;

};

}

#endif

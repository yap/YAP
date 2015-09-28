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

#ifndef yap_Particle_h
#define yap_Particle_h

#include "Amp.h"
#include "AmplitudeComponent.h"
#include "QuantumNumbers.h"

namespace yap {

class ParticleCombination;

/// \class Particle
/// \brief Abstract Particle base class.
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup Particle Particle-related classes

// keyword virtual is needed to solve diamond problem in DecayingParticle
class Particle : public virtual AmplitudeComponent
{
public:

    /// \name Constructor & clone
    /// @{

    /// Constructor
    Particle(const QuantumNumbers& q, double mass, std::string name);

    /// Clone
    virtual std::shared_ptr<Particle> clone() const = 0;

    /// @}

    /// Check consitency of object
    virtual bool consistent() const override;

    /// Access QuantumNumbers object
    const QuantumNumbers& quantumNumbers() const
    { return QuantumNumbers_; }

    QuantumNumbers& quantumNumbers()
    { return QuantumNumbers_; }

    /// \name Getters
    /// @{

    /// Get mass [GeV]
    double mass() const
    { return Mass_; }

    /// Get name
    std::string name() const
    { return Name_; }

    /// @}

    /// \name Setters
    /// @{

    /// Set mass [GeV]
    void setMass(double m)
    { Mass_ = m; }

    /// @}

    /// for internal use only
    virtual void setSymmetrizationIndexParents() = 0;

private:

    /// Quantum numbers of particle
    QuantumNumbers QuantumNumbers_;

    /// Mass [GeV]
    double Mass_;

    /// Name of particle
    std::string Name_;

};

}

#endif

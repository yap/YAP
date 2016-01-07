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

#ifndef yap_Resonance_h
#define yap_Resonance_h

#include "DecayingParticle.h"
#include "DataAccessor.h"
#include "DataPoint.h"
#include "MassShape.h"

#include <complex>
#include <memory>

namespace yap {

class DecayChannel;
class FinalStateParticle;
class InitialStateParticle;
class ParticleCombination;

/// \class Resonance
/// \brief Class for a particle that will decay and has a mass shape
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Particle

class Resonance : public DecayingParticle
{
public:

    /// Constructor
    /// \param q QuantumNumbers of resonance
    /// \param mass mass of resonance
    /// \param radialSize radial size of resonance
    /// \param massShape unique_ptr to MassShape of resonance
    Resonance(const QuantumNumbers& q, double mass, std::string name, double radialSize, std::unique_ptr<MassShape> massShape);

    /// Calculate complex amplitude
    virtual std::complex<double> amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc, unsigned dataPartitionIndex) const override
    { return DecayingParticle::amplitude(d, pc, dataPartitionIndex) * MassShape_->amplitude(d, pc, dataPartitionIndex); }

    /// Check consistency of object
    virtual bool consistent() const override;

    /// \name Getters
    /// @{

    /// access MassShape
    MassShape& massShape()
    { return *(MassShape_.get()); }

    /// access MassShape (const)
    const MassShape& massShape() const
    { return *(MassShape_.get()); }

    /// @}

    /// \name Setters
    /// @{

    /// set MassShape
    void setMassShape(std::unique_ptr<MassShape> massShape);

    using DecayingParticle::addChannel;

    /// Add a DecayChannel and set its parent to this DecayingParticle.
    /// \param c DecayingParticle takes ownership of c
    void addChannel(std::unique_ptr<DecayChannel> c) override;

    /// add symmetrizationIndex to SymmetrizationIndices_,
    /// also add to MassShape_
    virtual void addSymmetrizationIndex(std::shared_ptr<ParticleCombination> c) override;

    /// @}

    /// clear SymmetrizationIndices_
    virtual void clearSymmetrizationIndices() override;

protected:

    /// \return set of shared_ptr's to DataAccessor's
    virtual DataAccessorSet dataAccessors() override;

private:

    /// MassShape object
    std::shared_ptr<MassShape> MassShape_;

};

}

#endif

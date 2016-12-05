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

#ifndef yap_PoleMass_h
#define yap_PoleMass_h

#include "fwd/DataPoint.h"
#include "fwd/Parameter.h"
#include "fwd/ParticleCombination.h"
#include "fwd/ParticleTable.h"

#include "MassShape.h"

#include <complex>
#include <memory>

namespace yap {

/// \class PoleMass
/// \brief Class for pole-mass resonance shape
/// \author Daniel Greenwald
/// \ingroup MassShapes
///
/// Amplitude is 1 / (mass^2 - s)\n\n
/// mass is complex

class PoleMass : public MassShape
{
public:

    /// Constructor
    /// \param mass Mass of resonance [GeV]
    PoleMass(std::complex<double> mass);

    /// Constructor
    /// \param pde ParticleTableEntry to take mass and width from
    PoleMass(const ParticleTableEntry& pde);
    
    /// Get mass
    std::shared_ptr<ComplexParameter> mass() const
    { return Mass_; }

    /// Check consistency of object
    virtual bool consistent() const override;

    /// update the calculationStatus for a DataPartition
    virtual void updateCalculationStatus(StatusManager& D) const override;
    
    /// \return value for DataPoint and ParticleCombination
    /// \param d DataPoint
    /// \param pc shared_ptr to ParticleCombination
    virtual const std::complex<double> value(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const override;

    using MassShape::calculate;
    
    /// Calculate dynamic amplitude T for particular particle combination and store in each DataPoint in DataPartition
    /// \param D DataPartition to calculate on
    /// \param pc ParticleCombination to calculate for
    /// \param si SymmetrizationIndec to calculate for
    virtual void calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const;

protected:

    /// access cached dynamic amplitude
    const std::shared_ptr<ComplexCachedValue> T() const
    { return T_; }

private:

    /// cached dynamic amplitude
    std::shared_ptr<ComplexCachedValue> T_;

    /// Complex mass [GeV]
    std::shared_ptr<ComplexParameter> Mass_;

};

}

#endif

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

#ifndef yap_BreitWigner_h
#define yap_BreitWigner_h

#include "CachedValue.h"
#include "MassShape.h"

#include <complex>
#include <memory>

namespace yap {

class DataPoint;
class ParticleCombination;

/// \class BreitWigner
/// \brief Class for Breit-Wigner resonance shape
/// \author Daniel Greenwald
/// \ingroup MassShapes
///
/// Amplitude is 1 / (mass^2 - s - i*mass*width)\n\n

class BreitWigner : public MassShape
{
public:

    /// Constructor
    /// \param mass Mass of resonance [GeV]
    /// \param width Width of resonance [GeV]
    BreitWigner(double mass = -1, double width = -1);

    /// Calculate complex amplitude
    /// \param d DataPoint to calculate with
    /// \param pc (shared_ptr to) ParticleCombination to calculate for
    /// \param dataPartitionIndex partition index for parallelization
    virtual std::complex<double> amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc, unsigned dataPartitionIndex) const override;

    /// Set parameters from ParticleTableEntry
    /// \param entry ParticleTableEntry containing information to create mass shape object
    virtual void setParameters(const ParticleTableEntry& entry) override;

    /// \name Getters
    /// @{

    /// Get mass
    std::shared_ptr<RealParameter> mass() const
    { return Mass_; }

    /// Get width
    std::shared_ptr<RealParameter> width() const
    { return Width_; }

    /// @}

    /// \name Bookkeeping related
    /// @{

    virtual bool consistent() const override;

    /// @}

    virtual std::string data_accessor_type() const override
    {return "BreitWigner"; }

protected:

    /// borrow mass from owner
    virtual void borrowParametersFromResonance() override;

    /// set dependency on masses from model
    virtual void setDependenciesFromModel() override;

    std::shared_ptr<RealParameter> Mass_;  ///< [GeV]
    std::shared_ptr<RealParameter> Width_; ///< [GeV]

    std::shared_ptr<ComplexCachedDataValue> T_; ///< Mass-shape amplitude

};

}

#endif

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

#ifndef yap_Flatte_h
#define yap_Flatte_h

#include "fwd/Flatte.h"

#include "fwd/DataPartition.h"
#include "fwd/DataPoint.h"
#include "fwd/FinalStateParticle.h"
#include "fwd/Parameter.h"
#include "fwd/ParticleCombination.h"
#include "fwd/ParticleTable.h"
#include "fwd/StatusManager.h"

#include "MassShapeWithNominalMass.h"

#include <array>
#include <complex>
#include <memory>

namespace yap {

/// \class Flatte
/// \brief Class for Flatte resonance shape
/// \author Daniel Greenwald
/// \ingroup MassShapes
///
/// Amplitude is 1 / (mass^2 - s - i * sum_channels(coupling * phase-space factor)\n\n
/// phase space factor := 2 * breakup-momentum / m; may be complex
class Flatte : public MassShapeWithNominalMass
{
public:

    /// Constructor
    /// \param m mass [GeV]
    Flatte(double m);

    /// Constructor
    /// \param pde ParticleTableEntry to take mass from
    Flatte(const ParticleTableEntry& pde);

    /// Add FlatteChannel
    void add(FlatteChannel fc);

    /// Get FlatteChannel's
    const std::vector<FlatteChannel>& channels() const
    { return FlatteChannels_; }

    /// Check if a DecayChannel is valid for this MassShape; will throw if invalid.
    /// Cheks that decay is to a channel of the Flatte
    virtual void checkDecayChannel(const DecayChannel& c) const override;

    /// update the calculationStatus for a DataPartition
    virtual void updateCalculationStatus(StatusManager& D) const override;
    
    /// \return value for DataPoint and ParticleCombination
    /// \param d DataPoint
    /// \param pc shared_ptr to ParticleCombination
    virtual const std::complex<double> value(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const override;

    using MassShapeWithNominalMass::calculate;
    
    /// Calculate dynamic amplitude T for and store in each DataPoint in DataPartition
    /// \param D DataPartition to calculate on
    /// \param pc ParticleCombination to calculate for
    /// \param si SymmetrizationIndec to calculate for
    virtual void calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const override;


protected:

    /// access cached dynamic amplitude
    const std::shared_ptr<ComplexCachedValue> T() const
    { return T_; }

private:

    /// cached dynamic amplitude
    std::shared_ptr<ComplexCachedValue> T_;

    /// Flatte channels for width calculation
    std::vector<FlatteChannel> FlatteChannels_;

};

/// \struct FlatteChannel
/// \brief Stores information on channel used in calculating Flatte mass shape
/// \ingroup MassShapes
/// \author Daniel Greenwald
struct FlatteChannel {
    /// coupling constant [GeV^2]
    std::shared_ptr<NonnegativeRealParameter> Coupling;
    
    /// Particles of the channel
    std::array<std::shared_ptr<FinalStateParticle>, 2> Particles;
    
    /// constructor
    FlatteChannel(std::shared_ptr<NonnegativeRealParameter> coupling, FinalStateParticle& A, FinalStateParticle& B);
    
    /// constructor
    FlatteChannel(double coupling, FinalStateParticle& A, FinalStateParticle& B);

    /// constructor
    FlatteChannel(double coupling, const ParticleTableEntry& a, const ParticleTableEntry& b);
};

}

#endif

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

#include "fwd/DataPartition.h"
#include "fwd/DataPoint.h"
#include "fwd/Parameter.h"
#include "fwd/ParticleCombination.h"
#include "fwd/StatusManager.h"

#include "MassShapeWithNominalMass.h"

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

    struct FlatteChannel {
        /// [GeV]
        std::shared_ptr<RealParameter> Coupling;

        /// [GeV]
        std::shared_ptr<RealParameter> Mass;

        FlatteChannel(std::shared_ptr<RealParameter> coupling, std::shared_ptr<RealParameter> mass)
            : Coupling(coupling), Mass(mass) {}
    };

    /// Constructor
    Flatte() : MassShapeWithNominalMass() {}

    /// Add FlatteChannel
    void addChannel(std::shared_ptr<RealParameter> coupling, std::shared_ptr<RealParameter> mass);

    /// Add FlatteChannel
    void addChannel(double coupling, double mass);

    /// Get FlatteChannel's
    const std::vector<FlatteChannel>& channels() const
    { return FlatteChannels_; }

    /// Check consistency of object
    virtual bool consistent() const override;

    virtual std::string data_accessor_type() const override
    {return "Flatte"; }

protected:

    /// Calculate dynamic amplitude T for and store in each DataPoint in DataPartition
    /// \param D DataPartition to calculate on
    /// \param pc ParticleCombination to calculate for
    /// \param si SymmetrizationIndec to calculate for
    virtual void calculateT(DataPartition& D, const std::shared_ptr<ParticleCombination>& pc, unsigned si) const override;

    std::vector<FlatteChannel> FlatteChannels_;

};

}

#endif

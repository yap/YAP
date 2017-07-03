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

#ifndef yap_NonrelativisticConstantWidthBreitWigner_h
#define yap_NonrelativisticConstantWidthBreitWigner_h

#include "fwd/BlattWeisskopf.h"
#include "fwd/DataPartition.h"
#include "fwd/DecayChannel.h"
#include "fwd/ParticleCombination.h"
#include "fwd/ParticleTable.h"

#include "ConstantWidthBreitWigner.h"

#include <memory>

namespace yap {

/// \class NonrelativisticConstantWidthBreitWigner
/// \brief Class for Non-Relativistic Breit-Wigner resonance shape
/// \author Daniel Greenwald
/// \ingroup MassShapes
///
/// Amplitude is\n
/// \f$\frac{1}{M_{R} - \sqrt{s} - i*\Gamma}\f$\n
class NonrelativisticConstantWidthBreitWigner : public ConstantWidthBreitWigner
{
public:

    /// Constructor
    /// \param m Mass of resonance [GeV]
    /// \param w Width of resonance [GeV]
    NonrelativisticConstantWidthBreitWigner(double m, double w) : ConstantWidthBreitWigner(m, w) {}

    /// Constructor
    /// \param pde ParticleTableEntry to take mass and width from
    NonrelativisticConstantWidthBreitWigner(const ParticleTableEntry& pde) : ConstantWidthBreitWigner(pde) {}

    using ConstantWidthBreitWigner::calculate;
    
    /// Calculate dynamic amplitude T for and store in each DataPoint in DataPartition
    /// \param D DataPartition to calculate on
    /// \param pc ParticleCombination to calculate for
    /// \param si SymmetrizationIndec to calculate for
    virtual void calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const override;

};

}

#endif

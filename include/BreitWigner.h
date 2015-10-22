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
#include "CalculationStatus.h"
#include "MassShape.h"

#include <complex>

namespace yap {

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
    BreitWigner(double mass = -1, double width = -1);

    /// Calculate complex amplitude
    virtual std::complex<double> amplitude(DataPartition& d, std::shared_ptr<const ParticleCombination> pc) const override;

    /// \name Getters
    /// @{

    /// Get mass
    std::shared_ptr<RealParameter> mass() const
        { return std::dynamic_pointer_cast<RealParameter>(Parameters_[0]); }

    /// Get width
    std::shared_ptr<RealParameter> width() const
        { return std::dynamic_pointer_cast<RealParameter>(Parameters_[1]); }

    /// @}

    /// \name Amplitude related
    /// @{

    /// Calculate MassShape ampltude from squared mass;
    /// A = 1 / [Mass^2 - s - i * Mass * Width]
    /// \return amplitude evaluated at squared mass
    /// \param s squared mass to evaluate at
    virtual std::complex<double> calcAmplitudeS(double s) const override;

    /// @}

    /// \name Bookkeeping related
    /// @{

    virtual bool consistent() const override;

    /// @}

protected:

    std::shared_ptr<ComplexCachedValue> M2iMG_;                  // mass * mass - i * mass * width

    /// set owning resonance, borrow mass from owner
    virtual void borrowParametersFromResonance(Resonance* R) override;

};

}

#endif

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

#ifndef yap_MassShape_h
#define yap_MassShape_h

#include "AmplitudeComponent.h"
#include "DataAccessor.h"
#include "InitialStateParticle.h"
#include "ParameterSet.h"
#include "ParticleFactory.h"

#include <complex>
#include <memory>
#include <vector>

namespace yap {

class Resonance;

class ParticleCombination;

/// \class MassShape
/// \brief Abstract base class for all mass shapes
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup MassShapes Mass Shapes

class MassShape : public AmplitudeComponent, public DataAccessor
{
public:

    /// Constructor
    MassShape(std::vector<std::shared_ptr<ComplexParameter> > pars = {});

    /// Set parameters from ParticleTableEntry
    /// \param entry ParticleTableEntry containing information to create mass shape object
    /// \return Success of action
    virtual bool setParameters(const ParticleTableEntry& entry);

    const ParameterSet& parameters()
    { return Parameters_; }

    /// \name Bookkeeping related
    /// @{

    /// Check consistency of object
    virtual bool consistent() const override
    { return DataAccessor::consistent(); }

    /// @}

    /// \name Friends
    /// @{

    friend class Resonance;

    /// @}

protected:

    /// \name Amplitude related
    /// @{

    /// Calculate MassShape amplitude from squared mass
    /// \return amplitude evaluated at squared mass
    /// \param s squared mass to evaluate at
    virtual std::complex<double> calcAmplitudeS(double s) const = 0;

    /// @}

    /// override in inheriting classes to borrow parameters from Resonance
    virtual void borrowParametersFromResonance(Resonance* R)
    {}

    /// Parameters of the mass shape
    ParameterSet Parameters_;

};

}

#endif

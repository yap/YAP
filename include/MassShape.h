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

#include "Amp.h"
#include "AmplitudeComponent.h"
#include "DataAccessor.h"
#include "ParameterSet.h"

#include <memory>
#include <vector>

namespace yap {

class ParticleCombination;

/// \class MassShape
/// \brief Abstract base class for all mass shapes
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup MassShapes Mass Shapes

class MassShape : public AmplitudeComponent, public DataAccessor, public ParameterSet
{
public:

    /// Constructor
    MassShape();

    /// \name Amplitude related
    /// @{

    using AmplitudeComponent::amplitude;

    /// Calculate MassShape ampltude from squared mass
    /// \return amplitude evaluated at squared mass
    /// \param s squared mass to evaluate at
    virtual Amp amplitude(double s) = 0;

    /// @}

    /// \name Bookkeeping related
    /// @{

    /// Check consistency of object
    virtual bool consistent() const override
    { return DataAccessor::consistent() and ParameterSet::consistent(); }

    /// @}

};

}

#endif

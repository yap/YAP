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

#ifndef yap_AmplitudeComponent_h
#define yap_AmplitudeComponent_h

#include "CalculationStatus.h"
#include "CachedDataValue.h"
#include "DataPoint.h"
#include "Parameter.h"

#include <complex>
#include <memory>

namespace yap {

class ParticleCombination;

/// \name AmplitudeComponent
/// \brief Abstract base class for all objects implementing a (complex) amplitude
/// \author Johannes Rauch, Daniel Greenwald

class AmplitudeComponent
{
public:

    /// \name Constructors and assignment operators
    /// @{

    /// default constructor
    AmplitudeComponent() = default;

    /// default copy constructor
    AmplitudeComponent(const AmplitudeComponent&) = default;

    /// default move constructor
    AmplitudeComponent(AmplitudeComponent&&) = default;

    /// virtual default destructor
    virtual ~AmplitudeComponent() = default;

    /// default copy operator
    AmplitudeComponent& operator=(const AmplitudeComponent&) = default;

    /// default move operator
    AmplitudeComponent& operator=(AmplitudeComponent&&) = default;

    /// @}

    /// Check if AmplitudeComponent is consistent
    virtual bool consistent() const = 0;

    /// \return a set of parameters this AmplitudeComponent depends on
    /// to be overridden in the concrete AmplitudeComponent
    /// \return empty set
    virtual ParameterSet parametersItDependsOn()
    { return {}; }

    /// \return a list of CachedDataValues this AmplitudeComponent depends on
    /// to be overridden in the concrete AmplitudeComponent
    /// \return empty vector
    virtual CachedDataValueSet cachedDataValuesItDependsOn()
    { return {}; }

};

}

#endif

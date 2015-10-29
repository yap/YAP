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

    /// Default constructor
    AmplitudeComponent() {};

    /// Calculate complex amplitude
    virtual std::complex<double> amplitude(DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc, unsigned dataPartitionIndex) const = 0;

    /// Check if AmplitudeComponent is consistent
    virtual bool consistent() const = 0;

    /// \return a set of parameters this AmplitudeComponent depends on
    /// to be overridden in the concrete AmplitudeComponent
    /// \return empty set
    virtual ParameterSet ParametersItDependsOn()
    { return {}; }

    /// \return a list of CachedDataValues this AmplitudeComponent depends on
    /// to be overridden in the concrete AmplitudeComponent
    /// \return empty vector
    virtual CachedDataValueSet CachedDataValuesItDependsOn()
    { return {}; }
};

}

#endif

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
#include "DataPartition.h"
#include "DataPoint.h"

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

    /// \name Constructors, destructor, & operators
    /// @{

    /// Default constructor
    AmplitudeComponent();

    // Defaulted copy constructor
    // Defaulted move constructor
    // Defaulted destructor
    // Defaulted move assignment operator

    /// @}

    /// precalculate intermediary vaules needed for amplitude calculation,
    /// to be called before looping over DataPoints
    virtual void precalculate() {};

    /// Calculate complex amplitude
    virtual const std::complex<double>& amplitude(DataPartition& d, std::shared_ptr<const ParticleCombination> pc) const = 0;

    /// Check if AmplitudeComponent is consistent
    virtual bool consistent() const = 0;

};

}

#endif

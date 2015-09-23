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

#ifndef yap_AmplitudeComponentDataAccessor_h
#define yap_AmplitudeComponentDataAccessor_h

#include "AmplitudeComponent.h"
#include "DataAccessor.h"

#include "Amp.h"
#include "DataPoint.h"

#include <memory>

namespace yap {

class ParticleCombination;

/// \name AmplitudeComponentDataAccessor
/// \brief Abstract base class for all objects implementing an AmplitudeComponent and a DataAccessor
/// This class is meant for handling the caching of Amplitudes
/// \author Johannes Rauch, Daniel Greenwald

// keyword virtual is needed to solve diamond problem in DecayingParticle
class AmplitudeComponentDataAccessor : public virtual AmplitudeComponent, public DataAccessor
{
public:

    /// \name Constructors, destructor, & operators
    /// @{

    /// Default constructor
    AmplitudeComponentDataAccessor(InitialStateParticle* isp, ParticleCombination::Equiv* equiv = &ParticleCombination::equivBySharedPointer);

    // Defaulted copy constructor
    // Defaulted move constructor
    // Defaulted destructor
    // Defaulted move assignment operator

    /// @}

    /// Handle complex amplitude caching
    virtual const Amp& amplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc);

    /// Check if AmplitudeComponentDataAccessor is consistent
    virtual bool consistent() const
    { return DataAccessor::consistent(); }

protected :

    /// Calculate complex amplitude
    virtual Amp calcAmplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc) = 0;

};

}

#endif

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

#include "fwd/AmplitudeComponent.h"
#include "fwd/DataPoint.h"
#include "fwd/Parameter.h"
#include "fwd/ParticleCombination.h"

#include "RecalculableDataAccessor.h"
#include "StaticDataAccessor.h"
#include "VariableStatus.h"

#include <complex>
#include <memory>

namespace yap {

/// Base class for objects that are multiplicatively incorporated into an amplitude
/// \author Daniel Greenwald
class AmplitudeComponent
{
public:
    /// \return whether AmplitudeComponent can calculate for particular ParticleCombination
    virtual const bool validFor(const std::shared_ptr<ParticleCombination>& pc) const = 0;

    /// \return value for DataPoint and ParticleCombination
    /// \param d DataPoint
    /// \param pc shared_ptr to ParticleCombination
    /// must be overloaded in derived class
    virtual const std::complex<double> value(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const = 0;

    /// \return a VariableStatus for this AmplitudeComponent
    virtual const VariableStatus variableStatus() const = 0;
};

/// Base class for AmplitudeComponent's that are also StaticDataAccessor's
/// \author Daniel Greenwald
class StaticAmplitudeComponent : public AmplitudeComponent, public StaticDataAccessor
{
public:
    /// Constructor
    /// \param model owning Model
    /// \param equal ParticleCombination equality struct for determining index assignments
    StaticAmplitudeComponent(Model& m, const ParticleCombinationEqualTo& equal)
        : AmplitudeComponent(), StaticDataAccessor(m, equal) {}

    /// \return whether AmplitudeComponent can calculate for particular ParticleCombination
    virtual const bool validFor(const std::shared_ptr<ParticleCombination>& pc) const
    { return symmetrizationIndices().find(pc) != symmetrizationIndices().end(); }

    /// \return VariableStatus::fixed
    virtual const VariableStatus variableStatus() const override
    { return VariableStatus::fixed; }
};    


/// Base class for AmplitudeComponent's that are also RecalculableDataAccessor's
class RecalculableAmplitudeComponent : public AmplitudeComponent, public RecalculableDataAccessor
{
public:
    /// Constructor
    /// \param equal ParticleCombination equality struct for determining index assignments
    explicit RecalculableAmplitudeComponent(const ParticleCombinationEqualTo& equal)
        : RecalculableDataAccessor(equal) {}

    /// \return whether AmplitudeComponent can calculate for particular ParticleCombination
    virtual const bool validFor(const std::shared_ptr<ParticleCombination>& pc) const
    { return symmetrizationIndices().find(pc) != symmetrizationIndices().end(); }

    /// calls variable_status on parameters of object
    virtual const VariableStatus variableStatus() const override;
};


}

#endif

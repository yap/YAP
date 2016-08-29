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

#ifndef yap_RecalculableDataAccessor_h
#define yap_RecalculableDataAccessor_h

#include "fwd/RecalculableDataAccessor.h"

#include "fwd/DataPartition.h"
#include "fwd/DataPoint.h"
#include "fwd/Parameter.h"
#include "fwd/ParticleCombination.h"
#include "fwd/StatusManager.h"
#include "fwd/VariableStatus.h"

#include "DataAccessor.h"
#include "VariableStatus.h"

#include <complex>
#include <memory>
#include <set>

namespace yap {

/// Base class for all DataAccessor's that will need to be updated
/// \author Daniel Greenwald
class RecalculableDataAccessor : public DataAccessor
{
public:
    /// Constructor
    /// \param equal ParticleCombination equality struct for determining index assignments
    explicit RecalculableDataAccessor(const ParticleCombinationEqualTo& equal)
        : DataAccessor(equal) {}

    /// calculate for every data point in a DataPartition
    /// must be overloaded in derived class
    virtual void calculate(DataPartition& D) const = 0;

    /// update the calculationStatus for a DataPartition
    virtual void updateCalculationStatus(StatusManager& D) const = 0;

    /// set VariableStatus of all Parameters to unchanged (or leave fixed)
    void setParameterFlagsToUnchanged();

    /// \return Parameters of this object
    const ParameterSet& parameters() const
    { return Parameters_; }

protected:

    /// register with Model
    void virtual registerWithModel() override;

    /// add a parameter for evaluating variableStatus
    void addParameter(std::shared_ptr<ParameterBase> p)
    { Parameters_.insert(p); }

private:
    
    /// Parameters of object
    ParameterSet Parameters_;

};

}

#endif

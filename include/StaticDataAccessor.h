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

#ifndef yap_StaticDataAccessor_h
#define yap_StaticDataAccessor_h

#include "DataAccessor.h"
#include "Exceptions.h"
#include "ReportsInitialStateParticle.h"

namespace yap {

class InitialStateParticle;

/// \name StaticDataAccessor
/// \brief Base class for all data accessors that will only write to DataPoint once at initial data loading
/// \author Johannes Rauch, Daniel Greenwald
class StaticDataAccessor : public DataAccessor
{
public:

    /// Constructor
    /// \param isp Raw pointer to owning InitialStateParticle
    /// \param equiv ParticleCombination equivalence struct for determining index assignments
    StaticDataAccessor(InitialStateParticle* isp, ParticleCombination::Equiv* equiv = &ParticleCombination::equivBySharedPointer)
        : DataAccessor(equiv), InitialStateParticle_(isp)
    {
        if (!initialStateParticle())
            throw exceptions::Exception("InitialStateParticle unset", "StaticDataAccessor::StaticDataAccessor");
    }

    /// calculate cachedDataValues and store to DataPoint.
    /// Must be overriden in derived classes.
    virtual void calculate(DataPoint& d) = 0;

    /// does nothing since StaticDataAccessor's never update
    virtual void updateGlobalCalculationStatuses() override {}

    /// include const access to ISP
    using ReportsInitialStateParticle::initialStateParticle;

    /// \return Raw pointer to owning InitialStateParticle
    InitialStateParticle* initialStateParticle() override
    { return InitialStateParticle_; }

protected:

    /// does nothing, since StaticDataAccessor's ignore partitions
    virtual void setNumberOfDataPartitions(unsigned n) override {}

    /// does nothing, since StaticDataAccessor's never update
    virtual void resetCalculationStatuses(unsigned dataPartitionIndex) override {}

    /// does nothing, since StaticDataAccessor's never change
    virtual void setCachedDataValueFlagsToUnchanged(unsigned dataPartitionIndex) override {}

    /// does nothing, since StaticDataAccessor's never change
    virtual void setParameterFlagsToUnchanged() override {}

private:

    InitialStateParticle* InitialStateParticle_;

};

}

#endif

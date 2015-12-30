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

#include "Exceptions.h"

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
    { if (!initialStateParticle()) throw exceptions::InitialStateParticleUnset(); }

    /// Update global calculation statuses of all CachedDataValues
    /// [Does nothing since data is never to be updated!]
    virtual void updateGlobalCalculationStatuses() override
    { }

    /// calculate cachedDataValues and store to DataPoint
    virtual void calculate(DataPoint& d) = 0;

    /// include const access to ISP
    using BelongsToInitialStateParticle::initialStateParticle;

    /// \return Raw pointer to owning InitialStateParticle
    InitialStateParticle* initialStateParticle() override
    { return InitialStateParticle_; }

private:

    InitialStateParticle* InitialStateParticle_;

};

}

#endif

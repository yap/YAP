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

#include "fwd/Model.h"
#include "fwd/ParticleCombination.h"
#include "fwd/StaticDataAccessor.h"
#include "fwd/StatusManager.h"

#include "DataAccessor.h"

namespace yap {

/// \class StaticDataAccessor
/// \brief Base class for all data accessors that will only write to DataPoint once at initial data loading
/// \author Johannes Rauch, Daniel Greenwald
class StaticDataAccessor : public DataAccessor
{
public:

    /// Constructor
    /// \param model owning Model
    /// \param equal ParticleCombination equality struct for determining index assignments
    StaticDataAccessor(Model& m, const ParticleCombinationEqualTo& equal)
        : DataAccessor(equal), Model_(&m) {}

    /// calculate CachedValues, store to DataPoint, and update StatusManager.
    /// Must be overriden in derived classes.
    virtual void calculate(DataPoint& d, StatusManager& sm) const = 0;

    /// \return Raw pointer to owning Model
    const Model* model() const override
    { return Model_; }

protected:

    /// register with Model
    void virtual registerWithModel() override;

    /// add to model's StaticDataAccessors_
    void virtual addToStaticDataAccessors();

    /// access to the model's StaticDataAccessors_
    StaticDataAccessorVector& staticDataAccessors();

private:

    Model* Model_;

};

/// remove expired elements of set
void remove_expired(StaticDataAccessorVector& S);

}

#endif

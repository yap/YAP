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

#ifndef yap_BelongsToInitialStateParticle_h
#define yap_BelongsToInitialStateParticle_h

#include "ParticleCombination.h"

namespace yap {

class InitialStateParticle;

/// \name BelongsToInitialStateParticle
/// \brief Base class for all classes that need to know what
/// InitialStateParticle they belong to and report a vector of
/// ParticleCombination's
/// \author Johannes Rauch, Daniel Greenwald

class BelongsToInitialStateParticle
{
public:

    /// Default constructor
    BelongsToInitialStateParticle() : InitialStateParticle_(nullptr) {}

    /// constructor
    BelongsToInitialStateParticle(InitialStateParticle* isp) : InitialStateParticle_(isp) {}

    /// virtual destructor
    virtual ~BelongsToInitialStateParticle() = default;

    /// default copy constructor
    BelongsToInitialStateParticle(const BelongsToInitialStateParticle& other) = default;

    /// default move constructor
    BelongsToInitialStateParticle(BelongsToInitialStateParticle&& other) = default;

    /// default copy assignment operator
    BelongsToInitialStateParticle& operator=(const BelongsToInitialStateParticle& rhs) = default;

    /// default move assignment operator
    BelongsToInitialStateParticle& operator=(BelongsToInitialStateParticle&& rhs) = default;

    /// set raw pointer to initial state particle
    virtual void setInitialStateParticle(InitialStateParticle* isp)
    { InitialStateParticle_ = isp; }

    /// get raw pointer to initial state particle
    InitialStateParticle* initialStateParticle() const
    { return InitialStateParticle_; }

    /// \return vector of ParticleCombination's
    virtual ParticleCombinationVector particleCombinations() const = 0;

private:

    /// raw pointer to initial state particle
    InitialStateParticle* InitialStateParticle_;

};
}

#endif

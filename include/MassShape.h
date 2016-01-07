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

#ifndef yap_MassShape_h
#define yap_MassShape_h

#include "AmplitudeComponent.h"
#include "DataAccessor.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "ReportsInitialStateParticle.h"

namespace yap {

class Resonance;

/// \class MassShape
/// \brief Abstract base class for all mass shapes
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup MassShapes Mass Shapes
///
/// Inheriting classes (mass shapes) must implement
/// #AmplitudeComponent's amplitude(...) function.

class MassShape : public AmplitudeComponent, public DataAccessor
{
public:

    /// Constructor
    MassShape() : DataAccessor(&ParticleCombination::equivByOrderlessContent), Resonance_(nullptr)
    {}

    /// Set parameters from ParticleTableEntry
    /// Can be overloaded in inheriting classes
    /// \param entry ParticleTableEntry containing information to create mass shape object
    virtual void setParameters(const ParticleTableEntry& entry)
    { }

    /// Check consistency of object
    virtual bool consistent() const override;

    /// get raw pointer to owning resonance
    Resonance* resonance() const
    { return Resonance_; }

    /// include const access to ISP
    using ReportsInitialStateParticle::initialStateParticle;

    /// get raw pointer to initial state particle through resonance
    InitialStateParticle* initialStateParticle() override;

protected:

    /// Grant Resonance friendship, so it can set itself as owner
    friend class Resonance;

    /// Set raw pointer to owning Resonance.
    /// Calls borrowParametersFromResonance()
    void setResonance(Resonance* r)
    { Resonance_ = r; if (Resonance_) borrowParametersFromResonance(); }

    /// override in inheriting classes to borrow parameters from owning resonance
    virtual void borrowParametersFromResonance()
    {}

private:

    /// raw pointer to resonance that owns this mass shape
    Resonance* Resonance_;

};

}

#endif

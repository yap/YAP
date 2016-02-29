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
#include "ParticleFactory.h"

#include <memory>
#include <string>

namespace yap {

class Model;
class ParticleCombination;
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
    MassShape();

    /// Calculate complex amplitude.
    /// Must be overrided in derived classes.
    /// \param d DataPoint to calculate with
    /// \param pc (shared_ptr to) ParticleCombination to calculate for
    /// \param dataPartitionIndex partition index for parallelization
    virtual std::complex<double> amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc, unsigned dataPartitionIndex) const = 0;

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
    using DataAccessor::model;

    /// get raw pointer to Model through resonance
    Model* model() override;

    virtual std::string data_accessor_type() const override
    {return "MassShape"; }

    /// Grant Resonance friendship, so it can set itself as owner
    friend class Resonance;

protected:

    /// calls setDependenciesFromModel
    virtual void addToModel() override;

    /// Set raw pointer to owning Resonance.
    /// Calls borrowParametersFromResonance().
    virtual void setResonance(Resonance* r);

    /// override in inheriting class to set dependencies on values accessible through resonance
    virtual void setDependenciesFromResonance()
    {}

    /// override in inheriting class to set dependencies on values accessible through model
    virtual void setDependenciesFromModel()
    {}

    /// replace resonance's mass
    void replaceResonanceMass(std::shared_ptr<RealParameter> m);

    /// access cached dynamic amplitude
    std::shared_ptr<ComplexCachedDataValue> T()
    { return T_; }

    /// access cached dynamic amplitude (const)
    const std::shared_ptr<ComplexCachedDataValue> T() const
    { return const_cast<MassShape*>(this)->T(); }
    
private:

    /// raw pointer to resonance that owns this mass shape
    Resonance* Resonance_;

    /// cached dynamic amplitude
    std::shared_ptr<ComplexCachedDataValue> T_;
    
};

}

#endif

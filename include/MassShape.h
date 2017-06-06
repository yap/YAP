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

#include "fwd/CachedValue.h"
#include "fwd/DataPartition.h"
#include "fwd/DecayChannel.h"
#include "fwd/DecayingParticle.h"
#include "fwd/DecayTree.h"
#include "fwd/Model.h"
#include "fwd/Parameter.h"
#include "fwd/ParticleCombination.h"
#include "fwd/StatusManager.h"

#include "AmplitudeComponent.h"

#include <complex>
#include <memory>
#include <string>

namespace yap {

/// \class MassShape
/// \brief Abstract base class for all mass shapes
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup MassShapes Mass Shapes
class MassShape : public RecalculableAmplitudeComponent
{
public:

    /// Constructor
    MassShape();

    /// Calculate complex amplitudes for and store in each DataPoint in DataPartition;
    /// calls calculateT, which must be overrided in derived classes
    /// \param D DataPartition to calculate on
    virtual void calculate(DataPartition& D) const override final;

    /// Calculate dynamic amplitude T for particular particle combination and store in each DataPoint in DataPartition
    /// \param D DataPartition to calculate on
    /// \param pc ParticleCombination to calculate for
    /// \param si SymmetrizationIndec to calculate for
    virtual void calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const = 0;
    
    /// Check consistency of object
    virtual bool consistent() const override;

    /// get raw pointer to Model 
    const Model* model() const override
    { return Model_; }

    /// Check if a DecayChannel is valid for this MassShape; will throw if invalid.
    virtual void checkDecayChannel(const DecayChannel& c) const
    {}

    /// Grant DecayingParticle friendship, so it can set itself as
    /// owner and call addDecayChannel
    friend class DecayingParticle;

protected:

    void setModel(const Model& m);

    /// Give MassShape chance to perform operations based on the
    /// addition of a DecayChannel to its owner d
    virtual void addDecayChannel(DecayingParticle& d, std::shared_ptr<DecayChannel> c) {}

private:

    /// raw pointer to model
    const Model* Model_;

};

}

#endif

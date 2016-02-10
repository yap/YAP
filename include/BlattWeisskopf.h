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

#ifndef yap_BlattWeisskopf_h
#define yap_BlattWeisskopf_h

#include "AmplitudeComponent.h"
#include "CachedDataValue.h"
#include "DataAccessor.h"
#include "DataPoint.h"
#include "ReportsInitialStateParticle.h"

#include <complex>
#include <memory>

namespace yap {

class DecayingParticle;
class ParticleCombination;

/// \class BlattWeisskopf
/// \brief Class implementing BlattWeisskopf barrier factors
/// \author Johannes Rauch, Daniel Greenwald

class BlattWeisskopf : public AmplitudeComponent, public DataAccessor
{
public:

    /// \todo this is the inverse square, no?
    /// Calculate square of Blatt-Weisskopf factor (NOT the ratio of two Blatt-Weisskopf factors)
    /// \param l orbital angular momentum
    /// \param z square of (radial size * breakup momentum)
    static double F2(unsigned l, double z);

    /// Constructor
    /// \param L angular momentum of Blatt-Weisskopf barrier factor
    /// \param dp raw pointer to owning DecayingParticle
    BlattWeisskopf(unsigned L, DecayingParticle* dp);

    /// \return angular momentum
    unsigned L() const
    { return L_; }

    /// Calculate complex amplitude
    virtual double amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc, unsigned dataPartitionIndex) const;

    /// check consistency of object
    virtual bool consistent() const override;

    //virtual ParameterSet ParametersItDependsOn() override;

    virtual CachedDataValueSet CachedDataValuesItDependsOn() override
    { return {Fq_r, Fq_ab}; }

    /// include const access to Model
    using ReportsModel::model;

    /// \return raw pointer to Model through owning DecayingParticle
    Model* model() override;

    virtual std::string data_accessor_type() const override
    {return "BlattWeisskopf"; }

    /// grant friend status to DecayingParticle to call addParticleCombination
    friend class DecayingParticle;

private:

    /// raw pointer to owning DecayingParticle
    DecayingParticle* DecayingParticle_;

    /// angular momentum
    unsigned L_;

    /// Blatt-Weisskopf factor at nominal mass
    std::shared_ptr<RealCachedDataValue> Fq_r;

    /// Blatt-Weisskopf factor at data mass
    std::shared_ptr<RealCachedDataValue> Fq_ab;
};

}

#endif

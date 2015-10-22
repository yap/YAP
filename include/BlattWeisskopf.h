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
#include "CalculationStatus.h"
#include "DataAccessor.h"

#include <complex>
#include <memory>

namespace yap {

class DecayChannel;
class ParticleCombination;

/// \class BlattWeisskopf
/// \brief Class implementing BlattWeisskopf barrier factors
/// \author Johannes Rauch, Daniel Greenwald

class BlattWeisskopf : public AmplitudeComponent, public DataAccessor
{
public:

    /// Constructor
    BlattWeisskopf(DecayChannel* decayChannel);

    /// Calculate complex amplitude
    virtual std::complex<double> amplitude(DataPartition& d, const std::shared_ptr<const ParticleCombination>& pc) const override;

    /// check consistency of object
    virtual bool consistent() const override;

    /// Return DecayChannel this BlattWeisskopf belongs to
    DecayChannel* decayChannel() const {return DecayChannel_;}

    /// Calculate square of Blatt-Weisskopf factor (NOT the ratio of two Blatt-Weisskopf factors)
    static double F2(int twoL, double R2, double q2);

    //virtual std::vector<std::shared_ptr<ComplexParameter> > ParametersItDependsOn() override
    //{ return std::vector<std::shared_ptr<ComplexParameter> >(); }

    virtual std::vector<std::shared_ptr<CachedDataValue> > CachedDataValuesItDependsOn() override
    { return std::vector<std::shared_ptr<CachedDataValue> >{Fq_r, Fq_ab}; }

private:

    /// DecayChannel this BlattWeisskopf belongs to
    DecayChannel* DecayChannel_;

    std::shared_ptr<RealCachedDataValue> Fq_r;
    std::shared_ptr<RealCachedDataValue> Fq_ab;
};

}

#endif

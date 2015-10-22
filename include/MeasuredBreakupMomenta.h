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

#ifndef yap_MeasuredBreakupMomenta_h
#define yap_MeasuredBreakupMomenta_h

#include "CachedDataValue.h"
#include "DataAccessor.h"

namespace yap {

/// \class MeasuredBreakupMomenta
/// \brief Calculates, stores and gives access to breakup momenta (using measured masses)
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup SpinAmplitude

class MeasuredBreakupMomenta : public DataAccessor
{
public:

    /// Constructor
    MeasuredBreakupMomenta();

    /// Calculate breakup momenta for all possible symmetrization indices
    void calculate(DataPoint& d);

    /// Access squared breakup momentum
    /// \param d DataPoint to get data from
    /// \param pc ParticleCombination to return breakup momentum of
    double q2(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const
    { return Q2_.value(d, SymmetrizationIndices_.at(pc)); }

    /// Access breakup momentum
    /// \param d DataPoint to get data from
    /// \param pc ParticleCombination to return breakup momentum of
    double q(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const
    { return sqrt(q2(d, pc)); }

    /// Calculate breakup momentum from parent and daughter masses
    static double calcQ2(double m2_R, double m_a, double m_b);

protected:

    /// squared breakup momentum [GeV^2]
    RealCachedDataValue Q2_;
};

}

#endif

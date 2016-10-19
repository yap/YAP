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

#include "fwd/DataPoint.h"
#include "fwd/Model.h"
#include "fwd/ParticleCombination.h"

#include <cmath>
#include <memory>
#include <string>

namespace yap {

/// \namespace measured_breakup_momenta
/// \brief Calculates breakup momenta (using measured masses)
/// \author Johannes Rauch, Daniel Greenwald
namespace measured_breakup_momenta
{
    /// Access squared breakup momentum
    /// \param d DataPoint to get data from
    /// \param pc ParticleCombination to return breakup momentum of
    /// \param m Model to use FourMomenta manager from
    double q2(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc, const Model& m);

    /// Access breakup momentum
    /// \param d DataPoint to get data from
    /// \param pc ParticleCombination to return breakup momentum of
    inline double q(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc, const Model& m)
    { return sqrt(q2(d, pc, m)); }

    /// Calculate breakup momentum from parent and daughter masses
    /// \param m2_R squared mass of parent
    /// \param m_a mass of first daughter
    /// \param m_b mass of second daughter
    double q2(double m2_R, double m_a, double m_b);
}

}

#endif

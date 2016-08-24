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


#ifndef  yap_InvariantMassBinning_h
#define  yap_InvariantMassBinning_h

#include "fwd/CachedValue.h"

#include "StaticDataAccessor.h"

namespace yap {

/// \class InvariantMassBinning
/// \author Paolo Di Giglio.
/// Partitions the invariant-mass range in bins, whose (constant) lower bonuds
/// are stored in the class. Offers a function to calculate which bin the
/// invariant mass of a #ParticleCombination lies in.
class InvariantMassBinning : public StaticDataAccessor
{
public:
    /// Constructor.
    /// The #ParticleCombinationEqualTo is defaulted to #equal_by_orderless_content.
    explicit InvariantMassBinning(Model& m, const std::vector<double>& bins);

    /// Calculate which bin the invatiant mass of the #ParticleCombination's
    /// belong to.
    /// \param d  The #DataPoint to calculate into.
    /// \param sm The #StatusManager to update.
    virtual void calculate(DataPoint& d, StatusManager& sm) const override;
private:

    /// Partitioning of the mass axis.
    /// It contains the \f$n\f$ lower bounds of the bins plus
    /// the upper bound of the last bin:
    /// \f[
    /// [m_0,m_1) \cup [m_1,m_2) \dots [m_{n-1}, m_n).
    /// \f]
    const std::vector<double> Bins_;

    /// The bin the masses of the #ParticleCombination's belong to.
    std::shared_ptr<RealCachedValue> BinNumber_;
};

}

#endif

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

#include "fwd/InvariantMassBinning.h"

#include "fwd/CachedValue.h"
#include "fwd/DataPoint.h"
#include "fwd/Model.h"
#include "fwd/StatusManager.h"

#include "StaticDataAccessor.h"

#include <memory>
#include <vector>

namespace yap {

/// \brief Partitions the invariant-mass range in \f$n\f$ bins, whose (constant) lower
/// edges are stored in `BinLowEdges_`.
/// \details `BinLowEdges_` is a \f$(n+1)\f$-element vector whose entries partition the
/// axis in an closed-open fashion as follows:
/// \f[
/// (-\infty, m_0) \cup [m_0,m_1) \cup [m_1,m_2) \dots [m_{n-1}, m_n) \cup [m_n, +\infty).
/// \f]
/// The `0`-th value of `BinLowEdges_` corresponds to \f$m_0\f$; the `n`-th to \f$m_n\f$.
/// InvariantMassBinning has a function to calculate which bin the invariant mass of a
/// ParticleCombination lies in:
///  * In case of underflow (i.e. \f$m < m_0\f$), the bin will be set to `-1`.
///  * In case of overflow (i.e. \f$m > m_n\f$), the bin will be set to `n`.
///
/// \author Paolo Di Giglio.
class InvariantMassBinning : public StaticDataAccessor
{
public:
    /// Constructor.
    /// \param m         The owning Model.
    /// \param low_edges Low edges of the bins; the last element is the upper edge of the last bin.
    explicit InvariantMassBinning(Model& m, const std::vector<double>& low_edges);

    /// \brief Calculate which bin the invatiant mass of the
    /// ParticleCombination's belong to.
    /// \details The class description explains which bin value is used
    /// in case of overflow or underflow.
    /// \param d  The DataPoint to calculate into.
    /// \param sm The StatusManager to update.
    virtual void calculate(DataPoint& d, StatusManager& sm) const override;

    /// Access the partition.
    const std::vector<double>& lowEdges() const
    { return BinLowEdges_; }

    /// Access the bin number.
    /// \todo Enforce `const`ness on the pointed-to value!!
    const std::shared_ptr<RealCachedValue>& bin() const
    { return Bin_; }

private:

    /// Lower edges of the bins.
    const std::vector<double> BinLowEdges_;

    /// The bin the masses of the ParticleCombination's belong to.
    std::shared_ptr<RealCachedValue> Bin_;
};

/// Check if n corresponds to an underflow value.
/// \param n The index to check agains the bin underflow value.
inline const bool is_underflow(int n)
{ return n < 0; }

/// Check if n corresponds to an overflow value.
/// \param n The index to check agains the bin overflow value.
/// \param B The partition the bin overflow value has to be evaluated from.
inline const bool is_overflow(int n, const InvariantMassBinning& B)
{ return n >= static_cast<int>(B.lowEdges().size()) - 1; }

}

#endif

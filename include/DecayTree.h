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

#ifndef yap_DecayTree_h
#define yap_DecayTree_h

#include "fwd/DataPoint.h"
#include "fwd/DecayTree.h"
#include "fwd/FreeAmplitude.h"
#include "fwd/ParticleCombination.h"
#include "fwd/RecalculableDataAccessor.h"

#include <array>
#include <complex>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace yap {

/// \class DecayTree
/// \brief Class holding vectors of fixed and free amplitudes that define a decay tree
/// \author Johannes Rauch, Daniel Greenwald
class DecayTree
{
public:

    /// \typedef DaughterDecayTreeMap
    /// key is daughter index, value is decay tree of daughter
    using DaughterDecayTreeMap = std::map<unsigned, std::shared_ptr<DecayTree> >;

    /// Constructor
    /// \param two_M (twice) the spin projection of the parent particle
    /// \param two_m array of (twice) the spin projections of the daughters
    /// \param free_amp shared_ptr to ComplexParameter for the free amplitude
    explicit DecayTree(std::shared_ptr<FreeAmplitude> free_amp);

    /// \return amplitude evaluated for DataPoint over all ParticleCombinations of FreeAmplitude_'s DecayChannel
    /// \param d DataPoint
    const std::complex<double> amplitude(const DataPoint& d) const;

    /// \return amplitude evaluated for DataPoint for ParticleCombination
    /// \param d DataPoint
    /// \param pc ParticleCombination
    const std::complex<double> amplitude(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
    { return particleCombinationIndependentAmplitude(d) * particleCombinationDependentAmplitude(d, pc); }

    /// \return FreeAmplitude_
    const std::shared_ptr<FreeAmplitude>& freeAmplitude() const
    { return FreeAmplitude_; }

    const DaughterDecayTreeMap daughterDecayTrees() const
    { return DaughterDecayTrees_; }

    /// grant friend status to DecayChannel to call addDataAccessor
    friend class DecayChannel;

    /// grant friend status to DecayingParticle to call addDataAccessor
    friend class DecayingParticle;

    /// grant friend status to Resonance to call addDataAccessor
    friend class Resonance;

    /// convert to (multiline) string
    std::string asString(std::string offset = "") const;


protected:

    /// return product of all free amplitudes in this decay tree
    /// \param d DataPoint
    const std::complex<double> particleCombinationIndependentAmplitude(const DataPoint& d) const;

    /// return product of all particle-combination-dependent amplitudes in this tree
    /// \param d DataPoint
    /// \param pc shared_ptr<ParticleCombination>
    const std::complex<double> particleCombinationDependentAmplitude(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const;

    /// Set the DecayTree of the i'th daughter
    /// \param i index of daughter to set decay tree for
    /// \param dt shared_ptr to DecayTree to set
    void setDaughterDecayTree(unsigned i, std::shared_ptr<DecayTree> dt);

    /// set daughter spin projection
    /// \param two_m (twice) the spin projection
    void setDaughterSpinProjection(unsigned i, int two_m)
    { DaughtersTwoM_.at(i) = two_m; }

    /// Add a RecalculableDataAccessor
    void addDataAccessor(const RecalculableDataAccessor& rda);

private:

    /// (twice) daughter spin projections
    std::array<int, 2> DaughtersTwoM_;

    /// ComplexParameter of the free amplitude for the decay
    std::shared_ptr<FreeAmplitude> FreeAmplitude_;

    /// vector of RecalculableDataAccessor's
    std::vector<const RecalculableDataAccessor*> RecalculableDataAccessors_;

    /// map of daughter index -> daughter DecayTree
    DaughterDecayTreeMap DaughterDecayTrees_;

};

/// equality operator
inline bool operator==(const DecayTree& lhs, const DecayTree& rhs)
{ return lhs.freeAmplitude() == rhs.freeAmplitude() and lhs.daughterDecayTrees() == rhs.daughterDecayTrees(); }

/// \return Depth of DecayTree
unsigned depth(const DecayTree& DT);

/// \return sum of amplitudes of decay trees in a vector
const std::complex<double> amplitude(const DecayTreeVector& dtv, const DataPoint& d);

}

#endif

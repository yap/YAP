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

#ifndef yap_DataAccessor_h
#define yap_DataAccessor_h

#include "fwd/DataAccessor.h"

#include "fwd/CachedDataValue.h"
#include "fwd/DataPoint.h"
#include "fwd/Model.h"

#include "ParticleCombination.h"

#include <assert.h>
#include <memory>

namespace yap {

/// \name DataAccessor
/// \brief Abstract base class for all objects accessing DataPoint's
/// \author Johannes Rauch, Daniel Greenwald
class DataAccessor
{
public:

    /// Constructor
    /// \param equiv ParticleCombination equatility struct for determining index assignments
    DataAccessor(const ParticleCombination::Equiv& equiv = ParticleCombination::equivBySharedPointer);

    /// \return Equality struct
    const ParticleCombination::Equiv& equiv() const
    { return Equiv_; }

    /// \return index inside DataPoint structure that this DataAccessor accesses
    int index() const
    { return Index_; }

    /// \return index inside row of DataPoint for the requested ParticleCombination
    unsigned symmetrizationIndex(const std::shared_ptr<ParticleCombination>& c) const
    { return SymmetrizationIndices_.at(c); }

    /// \return SymmetrizationIndices_
    const ParticleCombinationMap<unsigned>& symmetrizationIndices() const
    { return SymmetrizationIndices_; }

    /// \return number of independent indices in SymmetrizationIndices_
    const unsigned nSymmetrizationIndices() const
    { return NIndices_; }

    /// \return vector of ParticleCombination's
    const ParticleCombinationVector& particleCombinations() const
    { return ParticleCombinationsCache_; }

    /// print ParticleCombination map
    void printParticleCombinations() const;

    /// \return CachedDataValueSet
    const CachedDataValueSet& cachedDataValues() const
    { return CachedDataValues_; }

    /// \return size of storage in data point (number of real values)
    unsigned size() const
    { return Size_; }

    /// Check consistency of object
    bool consistent() const;

    /// get raw pointer to Model (const)
    virtual const Model* model() const = 0;

    /// \return string denoting DataAccessor type
    /// \todo REMOVE
    virtual std::string data_accessor_type() const = 0;

    /// grant friend status to Model to access CachedDataValues_
    friend Model;

    /// grant friend status to CachedDataValue to call addCachedDataValue
    friend CachedDataValue;

protected:

    /// register with Model
    void virtual addToModel();

    /// add CachedDataValue
    void addCachedDataValue(std::shared_ptr<CachedDataValue> c);

    /// Increase storage
    /// \param n number of elements to increase by
    void increaseSize(unsigned n)
    { Size_ += n; }

    /// add ParticleCombination to SymmetrizationIndices_
    virtual void addParticleCombination(std::shared_ptr<ParticleCombination> pc);

    /// prune SymmetrizationIndices_ to only contain ParticleCombination's tracing back up the ISP
    virtual void pruneSymmetrizationIndices();

    /// set storage index used in DataPoint. Must be unique.
    void setIndex(size_t i)
    { Index_ = i; }

private:

    /// rebuild ParticleCombinations_
    void rebuildParticleCombinations();

    /// Object to check equality of symmetrizations for determining storage indices
    const ParticleCombination::Equiv& Equiv_;

    /// Map of indices for each used symmetrization stored with key = shared_ptr<ParticleCombination>
    ParticleCombinationMap<unsigned> SymmetrizationIndices_;

    /// Number of independent indices stored in SymmetrizationIndices_
    unsigned NIndices_;

    /// Vector of particle combinations. This is a cache kept for performance reasons,
    /// it must always be in sync with SymmetrizationIndices_
    ParticleCombinationVector ParticleCombinationsCache_;

    /// Set of CachedDataValues that have this DataAccessor as an owner
    CachedDataValueSet CachedDataValues_;

    /// number of real values stored per symm. index
    unsigned Size_;

    /// storage index used in DataPoint. Must be unique.
    int Index_;

};

/// remove expired elements of set
void removeExpired(DataAccessorSet& S);

}

#endif

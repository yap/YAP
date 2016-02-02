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

#include "CalculationStatus.h"
#include "CachedDataValue.h"
#include "ReportsInitialStateParticle.h"
#include "ReportsParticleCombinations.h"
#include "ParticleCombination.h"

#include "logging.h"

#include <memory>
#include <set>
#include <vector>

namespace yap {

class DataPoint;

/// \name DataAccessor
/// \brief Base class for all objects accessing DataPoint's
/// \author Johannes Rauch, Daniel Greenwald

class DataAccessor :
    public virtual ReportsInitialStateParticle,
    public virtual ReportsParticleCombinations
{
public:

    /// \name Constructors, destructor, & operators
    /// @{

    /// Constructor
    /// \param equiv ParticleCombination equivalence struct for determining index assignments
    DataAccessor(ParticleCombination::Equiv* equiv = &ParticleCombination::equivBySharedPointer);

    /// Copy constructor is deleted, since we don't need it and implementing it for all deriving classes would be too complicated
    DataAccessor(const DataAccessor& other) = delete;
    /// Copy operator is deleted
    DataAccessor& operator=(const DataAccessor& other) = delete;

    /// Destructor
    virtual ~DataAccessor();

    // Defaulted move constructor
    DataAccessor(DataAccessor&& other) = default;

    // Defaulted move assignment operator
    DataAccessor& operator=(DataAccessor&& other) = default;

    /// @}

    /// \name Access to indices
    /// @{

    /// \return index inside DataPoint structure that this DataAccessor accesses
    int index() const
    { return Index_; }

    /// \return if the given ParticleCombination is in SymmetrizationIndices_
    bool hasParticleCombination(const std::shared_ptr<ParticleCombination>& c) const
    { return SymmetrizationIndices_.find(c) != SymmetrizationIndices_.end(); }

    /// \return if the given ParticleCombination is in SymmetrizationIndices_
    /// \param c ParticleCombination to look for equivalent of
    /// \param equiv ParticleCombination::Equiv object for checking equivalence
    bool hasParticleCombination(const std::shared_ptr<ParticleCombination>& c,
                                const ParticleCombination::Equiv& equiv) const;

    /// \return index inside row of DataPoint for the requested ParticleCombination
    unsigned symmetrizationIndex(const std::shared_ptr<ParticleCombination>& c) const
    { return SymmetrizationIndices_.at(c); }

    /// \return SymmetrizationIndices_
    const ParticleCombinationMap<unsigned>& symmetrizationIndices() const
    { return SymmetrizationIndices_; }

    /// \return maximum index of SymmetrizationIndices_
    /// -1 means empty
    int maxSymmetrizationIndex() const;

    /// \return vector of ParticleCombination's
    ParticleCombinationVector particleCombinations() const override;

    /// print ParticleCombination map
    void printParticleCombinations() const;

    /// @}

    /// \return size of storage in data point (number of real values)
    unsigned size() const
    { return Size_; }

    /// Increase storage
    /// \param n number of elements to increase by
    void increaseSize(unsigned n)
    { Size_ += n; }

    /// Check consistency of object
    bool consistent() const;

    /// grant friend status to InitialStateParticle
    friend class InitialStateParticle;

    /// grant friend status to CachedDataValue to call addCachedDataValue
    friend class CachedDataValue;

protected:

    /// register with InitialStateParticle
    void virtual addToInitialStateParticle();

    /// add CachedDataValue
    void addCachedDataValue(std::shared_ptr<CachedDataValue> c)
    { CachedDataValues_.insert(c); }

    /// add ParticleCombination to SymmetrizationIndices_
    virtual void addParticleCombination(std::shared_ptr<ParticleCombination> pc) override;

    /// prune SymmetrizationIndices_ to only contain ParticleCombination's tracing back up the ISP
    virtual void pruneSymmetrizationIndices();

    /// set storage index used in DataPoint. Must be unique.
    void setIndex(size_t i)
    { Index_ = i; }

    /// set number of data partitions into all members of CachedDataValues_
    virtual void setNumberOfDataPartitions(unsigned n);

    /// Update global calculation statuses of all CachedDataValues
    virtual void updateGlobalCalculationStatuses();

    /// resets CalculationStatus'es for all CachedDataValues_
    virtual void resetCalculationStatuses(unsigned dataPartitionIndex);

    /// set all VariableStatus flags to kUnchanged (or leave at kFixed) for CachedDataValues_
    virtual void setCachedDataValueFlagsToUnchanged(unsigned dataPartitionIndex);

    /// set all VariableStatus flags to kUnchanged (or leave at
    /// kFixed) for all Parameters that CachedDataValues_ depend on
    virtual void setParameterFlagsToUnchanged();

    /// \name Data access
    /// \todo Is this access ever used?
    /// @{

    /// Access a data point's data
    std::vector<double>& data(DataPoint& d, unsigned i) const;

    /// Access a data point's data (const)
    const std::vector<double>& data(const DataPoint& d, unsigned i) const
    { return const_cast<DataAccessor*>(this)->data(d, i); }

    /// @}

private:

    /// Object to check equality of symmetrizations for determining storage indices
    ParticleCombination::Equiv* Equiv_;

    /// Map of indices for each used symmetrization stored with key = shared_ptr<ParticleCombination>
    ParticleCombinationMap<unsigned> SymmetrizationIndices_;

    /// Set of CachedDataValues that have this DataAccessor as an owner
    CachedDataValueSet CachedDataValues_;

    /// number of real values stored per symm. index
    unsigned Size_;

    /// storage index used in DataPoint. Must be unique.
    int Index_;

};

/// \typedef DataAccessorSet
using DataAccessorSet = std::set<DataAccessor*>;
// using DataAccessorSet = std::set<std::shared_ptr<DataAccessor>, std::owner_less<std::shared_ptr<DataAccessor> > >;

/// remove expired elements of set
void removeExpired(DataAccessorSet& S);

std::string data_accessor_type(const DataAccessor* D);

}

#endif

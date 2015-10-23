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
#include "DataPartition.h"
#include "DataPoint.h"
#include "ParticleCombination.h"

#include "logging.h"

#include <complex>
#include <map>
#include <vector>

namespace yap {

class CachedDataValue;
class InitialStateParticle;

/// \name DataAccessor
/// \brief Base class for all objects accessing DataPoint's
/// \author Johannes Rauch, Daniel Greenwald

class DataAccessor
{
public:

    /// \name Constructors, destructor, & operators
    /// @{

    /// Constructor
    /// \param equiv ParticleCombination equivalence struct for determining index assignments
    DataAccessor(ParticleCombination::Equiv* equiv = &ParticleCombination::equivBySharedPointer);

    /// Copy constructor is deleted, since we don't need it and implementing it for all deriving classes would be too complicated
    DataAccessor(const DataAccessor& other) = delete;

    /// Destructor
    virtual ~DataAccessor();

    // Defaulted move constructor
    // Defaulted move assignment operator

    /// @}

    /// \name Data accessor friends
    /// @{

    friend class InitialStateParticle;

    /// @}

    /// \name Access to indices
    /// @{

    /// \return index inside DataPoint structure that this DataAccessor accesses
    unsigned index() const
    { return Index_; }

    /// \return if the given ParticleCombination is in SymmetrizationIndices_
    bool hasSymmetrizationIndex(const std::shared_ptr<const ParticleCombination>& c) const
    { return SymmetrizationIndices_.count(c); }

    /// \return index inside row of DataPoint for the requested symmetrization
    unsigned symmetrizationIndex(const std::shared_ptr<const ParticleCombination>& c) const
    { return SymmetrizationIndices_.at(c); }

    /// \return SymmetrizationIndices_
    const std::map<std::shared_ptr<const ParticleCombination>, unsigned, std::owner_less<std::shared_ptr<const ParticleCombination> > >& symmetrizationIndices() const
    { return SymmetrizationIndices_; }

    /// \return maximum index of SymmetrizationIndices_
    unsigned maxSymmetrizationIndex() const;

    /// \return list of all ParticleCombinations
    std::vector<std::shared_ptr<const ParticleCombination> > particleCombinations() const;

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

    /// add CachedDataValue
    void addCachedDataValue(CachedDataValue* c);

    /// \name Symmetrization functions
    /// @{

    /// add symmetrizationIndex to SymmetrizationIndices_
    virtual void addSymmetrizationIndex(std::shared_ptr<const ParticleCombination> c);

    /// clear SymmetrizationIndices_
    virtual void clearSymmetrizationIndices()
    { SymmetrizationIndices_.clear(); }

    /// @}

    /// \name Data access
    /// @{

    /// Access a data point's data (by friendship)
    std::vector<double>& data(DataPoint& d, unsigned i) const;

#ifdef ELPP_DISABLE_DEBUG_LOGS
    /// Access a data point's data (by friendship) (const)
    const std::vector<double>& data(const DataPoint& d, unsigned i) const
    { return d.Data_[Index_][i]; }
#else
    /// Access a data point's data (by friendship) (const)
    const std::vector<double>& data(const DataPoint& d, unsigned i) const
    { return d.Data_.at(Index_).at(i); }
#endif

    /// Get pointer to the initial state particle
    InitialStateParticle* initialStateParticle() const
    { return InitialStateParticle_; }

    /// @}

    /// \name Setters
    /// @{

    /// Set pointer to initial state particle
    virtual void setInitialStateParticle(InitialStateParticle* isp);

    /// @}



    /// \name calculation statuses
    /// @{

    /// \return #CalculationStatus of symmetrization index and data-partition index
    /// \param pc shared pointer to #ParticleCombination to check status of
    /// \param symmetrizationIndex index of symmetrization to check status of
    /// \param dataPartitionIndex index of dataPartitionIndex to check status of
    virtual CalculationStatus calculationStatus(const std::shared_ptr<const ParticleCombination>& pc, unsigned symmetrizationIndex, unsigned dataPartitionIndex) const;

    /// \return #CalculationStatus of symmetrization index and data-partition index
    /// \param pc_symInd pair of shared pointer to #ParticleCombination and symmetrization index
    /// \param dataPartitionIndex index of dataPartitionIndex to check status of
    CalculationStatus calculationStatus(std::pair<std::shared_ptr<const ParticleCombination>, unsigned> pc_symInd, unsigned dataPartitionIndex = 0) const
    { return calculationStatus(pc_symInd.first, pc_symInd.second, dataPartitionIndex); }

    /// \return #CalculationStatus of symmetrization index and data-partition index
    /// \param pc shared pointer to #ParticleCombination to check status of
    /// \param dataPartitionIndex index of dataPartitionIndex to check status of
    CalculationStatus calculationStatus(const std::shared_ptr<const ParticleCombination>& pc, unsigned dataPartitionIndex = 0) const
    { return calculationStatus(pc, symmetrizationIndex(pc), dataPartitionIndex); }

    /// @}

protected:

    void setIndex(unsigned i)
    { Index_ = i; }

private:

    /// pointer to the initial state particle for access to FourMomenta, HelicityAngles etc.
    InitialStateParticle* InitialStateParticle_;

    /// Object to check equality of symmetrizations for determining storage indices
    ParticleCombination::Equiv* Equiv_;

    /// Map of indices for each used symmetrization stored with key = shared_ptr<ParticleCombination>
    std::map<std::shared_ptr<const ParticleCombination>, unsigned, std::owner_less<std::shared_ptr<const ParticleCombination> > > SymmetrizationIndices_;

    /// Set of CachedDataValues that have this DataAccessor as an owner
    std::set<CachedDataValue*> CachedDataValues_;

    /// number of real values stored per symm. index
    unsigned Size_;

    /// storage index used in DataPoint. Must be unique.
    unsigned Index_;

};

}

#endif

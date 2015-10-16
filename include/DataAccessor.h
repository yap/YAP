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

    /// Copy constructor
    /// \todo Do we need a DataAccessor copy constructor? (Currently copies size, but not CachedDataValue's)
    DataAccessor(const DataAccessor& other);

    /// Destructor
    virtual ~DataAccessor();

    // Defaulted move constructor
    // Defaulted move assignment operator

    /// @}

    /// \name Access to indices
    /// @{

    /// \return index inside DataPoint structure that this DataAccessor accesses
    unsigned index() const
    { return Index_; }

    /// \return if the given ParticleCombination is in SymmetrizationIndices_
    bool hasSymmetrizationIndex(std::shared_ptr<const ParticleCombination> c) const;

    /// \return index inside row of DataPoint for the requested symmetrization
    unsigned symmetrizationIndex(std::shared_ptr<const ParticleCombination> c) const
    { return SymmetrizationIndices_.at(c); }

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
    void addCachedDataValue(CachedDataValue* c)
    { CachedDataValues_.push_back(c); }

    /// \name Symmetrization functions
    /// @{

    /// add symmetrizationIndex to SymmetrizationIndices_
    virtual void addSymmetrizationIndex(std::shared_ptr<const ParticleCombination> c);

    /// clear SymmetrizationIndices_
    virtual void clearSymmetrizationIndices();

    /// @}

    /// \name Data access
    /// @{

    /// Access a data point's data (by friendship)
    std::vector<double>& data(DataPoint& d, unsigned i) const;

#ifdef ELPP_DISABLE_DEBUG_LOGS
    /// Access a data point's data (by friendship) (const)
    const std::vector<double>& data(const DataPoint& d, unsigned i) const
    { return d.Data_[Index_][i]; }

    std::complex<double>& cachedAmplitude(DataPoint& d, unsigned i) const
    { return d.CachedAmplitudes_[Index_][i]; }

    const std::complex<double>& cachedAmplitude(const DataPoint& d, unsigned i) const
    { return d.CachedAmplitudes_[Index_][i]; }
}
#else
    /// Access a data point's data (by friendship) (const)
    const std::vector<double>& data(const DataPoint& d, unsigned i) const
    { return d.Data_.at(Index_).at(i); }

    std::complex<double>& cachedAmplitude(DataPoint& d, unsigned i) const
    { return d.CachedAmplitudes_.at(Index_).at(i); }

    const std::complex<double>& cachedAmplitude(const DataPoint& d, unsigned i) const
    { return d.CachedAmplitudes_.at(Index_).at(i); }
#endif

CalculationStatus calculationStatus(const DataPartition& d, std::shared_ptr<const ParticleCombination> c) const
{ return d.CalculationStatusesDataSet(Index_, symmetrizationIndex(c)); }

CalculationStatus calculationStatus(const DataPartition& d, unsigned symIndex) const
{ return d.CalculationStatusesDataSet(Index_, symIndex); }


/// Get pointer to the initial state particle
InitialStateParticle* initialStateParticle() const;

/// @}

/// \name Setters
/// @{

/// Set CalculationStatus for particleCombination
void setCalculationStatus(DataPartition& d, std::shared_ptr<const ParticleCombination> c, CalculationStatus stat) const
{ d.CalculationStatusesDataSet(index(), symmetrizationIndex(c)) = stat; }

void setCalculationStatus(DataPartition& d, unsigned symIndex, CalculationStatus stat) const
{ d.CalculationStatusesDataSet(Index_, symIndex) = stat; }

/// Set pointer to initial state particle
virtual void setInitialStateParticle(InitialStateParticle* isp);

/// @}

protected:

friend class InitialStateParticle;

void setIndex(unsigned i)
{ Index_ = i; }

/// pointer to the initial state particle for access to FourMomenta, HelicityAngles etc.
InitialStateParticle* InitialStateParticle_;

/// Object to check equality of symmetrizations for determining storage indices
ParticleCombination::Equiv* Equiv_;

/// Map of indices for each used symmetrization stored with key = shared_ptr<ParticleCombination>
std::map<std::shared_ptr<const ParticleCombination>, unsigned, std::owner_less<std::shared_ptr<const ParticleCombination> > > SymmetrizationIndices_;

std::vector<CachedDataValue*> CachedDataValues_;

/// number of real values stored per symm. index
unsigned Size_;

private:

/// storage index used in DataPoint. Must be unique.
unsigned Index_;

};

}

#endif

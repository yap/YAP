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

#ifndef yap_CachedDataValue_h
#define yap_CachedDataValue_h

#include "CalculationStatus.h"
#include "Constants.h"
#include "FourVector.h"
#include "Parameter.h"
#include "VariableStatus.h"

#include <memory>
#include <set>
#include <utility>
#include <vector>

namespace yap {

class CachedDataValue;
class DataAccessor;
class DataPoint;
class ParticleCombination;

/// \typedef CachedDataValueSet
/// \ingroup Data
/// \ingroup Cache
using CachedDataValueSet = std::set<std::shared_ptr<CachedDataValue> >;

/// \typedef CachedDataValueDaughterIndexPair
/// \ingroup Data
/// \ingroup Cache
using CachedDataValueDaughterIndexPair = std::pair<std::shared_ptr<CachedDataValue>, unsigned>;

/// \typedef CachedDataValueDaughterIndexPairSet
/// \ingroup Data
/// \ingroup Cache
using CachedDataValueDaughterIndexPairSet = std::set<CachedDataValueDaughterIndexPair>;

/// \class CachedDataValueBase
/// \brief Class for managing cached values inside a #DataPoint
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
/// \ingroup Cache
class CachedDataValue : public std::enable_shared_from_this<CachedDataValue>
{
public:

    /// check consistency of object
    bool consistent() const;

    /// \name Managing dependencies
    /// @{

    /// add Parameters this CachedDataValue depends on
    void addDependencies(ParameterSet deps)
    { for (auto& dep : deps) addDependency(dep); }

    /// add Parameter this CachedDataValue depends on
    void addDependency(std::shared_ptr<ParameterBase> dep)
    { ParametersItDependsOn_.insert(dep); }

    /// add CachedDataValue's this CachedDataValue depends on
    void addDependencies(CachedDataValueSet deps)
    { for (auto& dep : deps) addDependency(dep); }

    /// add CachedDataValue's of a daughter this CachedDataValue depends on
    void addDependencies(CachedDataValueDaughterIndexPairSet deps)
    { for (auto& dep : deps) addDependency(dep); }

    /// add CachedDataValue this CachedDataValue depends on
    void addDependency(std::shared_ptr<CachedDataValue> dep)
    { CachedDataValuesItDependsOn_.insert(dep); }

    /// add CachedDataValue of a daughter this CachedDataValue depends on
    /// \param pcDaughterIndex Index of the daughter in the list of ParticleCombination's daughters
    void addDependency(std::shared_ptr<CachedDataValue> dep, unsigned pcDaughterIndex)
    { DaughterCachedDataValuesItDependsOn_.insert(std::make_pair(dep, pcDaughterIndex)); }

    /// add CachedDataValue of a daughter this CachedDataValue depends on
    void addDependency(CachedDataValueDaughterIndexPair dep)
    { DaughterCachedDataValuesItDependsOn_.insert(dep); }

    /// remove dependency
    void removeDependency(std::shared_ptr<ParameterBase> dep);

    /// remove dependencies
    void removeDependencies(ParameterSet deps)
    { for (auto& dep : deps) removeDependency(dep); }

    /// check for dependency
    bool dependsOn(std::shared_ptr<ParameterBase> dep) const
    { return ParametersItDependsOn_.find(dep) != ParametersItDependsOn_.end(); }

    /// check for dependency
    bool dependsOn(std::shared_ptr<CachedDataValue> dep) const
    { return CachedDataValuesItDependsOn_.find(dep) != CachedDataValuesItDependsOn_.end(); }

    /// @}

    /// \name Getters
    /// @{

    /// \return raw pointer owning DataAccessor
    DataAccessor* owner() const
    { return Owner_; }

    /// overload and hide #CachedValue::calculationStatus
    /// \return #CalculationStatus of symmetrization index and data-partition index
    /// \param pc shared pointer to #ParticleCombination to check status of
    /// \param dataPartitionIndex index of dataPartitionIndex to check status of
    virtual CalculationStatus calculationStatus(const std::shared_ptr<ParticleCombination>& pc, unsigned symmetrizationIndex,  unsigned dataPartitionIndex) const
#ifdef ELPP_DISABLE_DEBUG_LOGS
    { return CalculationStatus_[dataPartitionIndex][symmetrizationIndex]; }
#else
    { return CalculationStatus_.at(dataPartitionIndex).at(symmetrizationIndex); }
#endif

    /// get global CalculationStatus
    /// \return #CalculationStatus of symmetrization index and data-partition index
    /// \param symmetrization index
    CalculationStatus globalCalculationStatus(unsigned symmetrizationIndex)
    { return GlobalCalculationStatus_[symmetrizationIndex]; }

    /// get global CalculationStatus
    /// \return #CalculationStatus of symmetrization index and data-partition index
    /// \param pc shared pointer to #ParticleCombination to check status of
    CalculationStatus globalCalculationStatus(const std::shared_ptr<ParticleCombination>& pc);

    /// \return VariableStatus for symmetrization index and data-partition index
    /// \param symmetrizationIndex index of symmetrization to check status of
    /// \param dataPartitionIndex index of dataPartitionIndex to check status of
    VariableStatus variableStatus(unsigned symmetrizationIndex, unsigned dataPartitionIndex) const
#ifdef ELPP_DISABLE_DEBUG_LOGS
    { return VariableStatus_[dataPartitionIndex][symmetrizationIndex]; }
#else
    { return VariableStatus_.at(dataPartitionIndex).at(symmetrizationIndex); }
#endif

    /// Get value from #DataPoint for particular symmetrization
    /// \param index index of value to get from within cached value (must be less than #Size_)
    /// \param d #DataPoint to get value from
    /// \param symmetrizationIndex index of symmetrization to grab from
    /// \return Value of CachedDataValue inside the data point
    double value(unsigned index, const DataPoint& d, unsigned symmetrizationIndex) const;

    /// \return Size of cached value (number of real elements)
    virtual unsigned size() const
    { return Size_; }

    /// \return number of Symmetrizations
    unsigned numberOfSymmetrizations() const
    { return GlobalCalculationStatus_.size(); }

    /// \return number DataPartitions
    unsigned numberOfDataPartitions() const
    { return CalculationStatus_.size(); }

    /// @}

    /// \name Setters
    /// @{

    /// set VariableStatus for symmetrization index and data-partition index
    /// \param stat VariableStatus to set to
    /// \param symmetrizationIndex index of symmetrization to set status of
    /// \param dataPartitionIndex index of dataPartitionIndex to set status of
    void setVariableStatus(VariableStatus stat, unsigned symmetrizationIndex, unsigned dataPartitionIndex)
    { VariableStatus_[dataPartitionIndex][symmetrizationIndex] = stat; }

    /// set all variable statuses
    /// \param stat VariableStatus to set to
    /// \param dataPartitionIndex index of dataPartitionIndex to set status of
    void setVariableStatus(VariableStatus stat, unsigned dataPartitionIndex)
    { for (VariableStatus& s : VariableStatus_[dataPartitionIndex]) s = stat; }


    /// set CalculationStatus for symmetrization index and data-partition index
    /// \param stat VariableStatus to set to
    /// \param symmetrizationIndex index of symmetrization to set status of
    /// \param dataPartitionIndex index of dataPartitionIndex to set status of
    void setCalculationStatus(CalculationStatus stat, unsigned symmetrizationIndex, unsigned dataPartitionIndex)
    { CalculationStatus_[dataPartitionIndex][symmetrizationIndex] = stat; }

    /// set CalculationStatus for ParticleCombination and data-partition index
    /// \param stat VariableStatus to set to
    /// \param ParticleCombination to set status of
    /// \param dataPartitionIndex index of dataPartitionIndex to set status of
    void setCalculationStatus(CalculationStatus stat, const std::shared_ptr<ParticleCombination>& pc, unsigned dataPartitionIndex);

    /// set all calculation statuses
    /// \param stat CalculationStatus to set to
    /// \param dataPartitionIndex index of dataPartitionIndex to set status of
    void setCalculationStatus(CalculationStatus stat, unsigned dataPartitionIndex)
    { for (CalculationStatus& s : CalculationStatus_[dataPartitionIndex]) s = stat; }


    /// set global CalculationStatus for symmetrization index
    /// \param stat VariableStatus to set to
    /// \param symmetrizationIndex index of symmetrization to set status of
    void setGlobalCalculationStatus(CalculationStatus stat,  unsigned symmetrizationIndex)
    { GlobalCalculationStatus_[symmetrizationIndex] = stat; }

    /// set global CalculationStatus for symmetrization index
    /// \param stat VariableStatus to set to
    /// \param ParticleCombination to set status of
    void setGlobalCalculationStatus(CalculationStatus stat, const std::shared_ptr<ParticleCombination>& pc);

    /// Set value into #DataPoint for particular symmetrization
    /// (No update to VariableStatus or CalculationStatus is made!)
    /// \param index index of value to get from within cached value (must be less than #Size_)
    /// \param val Value to set to
    /// \param d #DataPoint to update
    /// \param symmetrizationIndex index of symmetrization to apply to
    void setValue(unsigned index, double val, DataPoint& d, unsigned symmetrizationIndex) const;

    /// @}

    /// update the global calculation status, depending on everything it depends on
    /// \param pc Shared pointer (const reference) to a Particle combination
    void updateGlobalCalculationStatus(const std::shared_ptr<ParticleCombination>& pc);

    /// update the global calculation status, depending on everything it depends on
    /// \param pc Shared pointer (const reference) to a Particle combination
    /// \param symmetrizationIndex index of symmetrization to apply to
    void updateGlobalCalculationStatus(const std::shared_ptr<ParticleCombination>& pc, unsigned symmetrizationIndex);

    /// reset the CalculationStatus'es for the dataPartitionIndex to the GlobalCalculationStatus_
    /// to be called before calculating the amplitude for a new dataPoint
    void resetCalculationStatus(unsigned dataPartitionIndex)
    { CalculationStatus_[dataPartitionIndex] = GlobalCalculationStatus_; }

    /// grant friend status to DataAccessor to set itself owner
    friend class DataAccessor;

protected:

    /// Constructor (protected)
    CachedDataValue(unsigned size, ParameterSet pars = {}, CachedDataValueSet vals = {});

    /// set the owning DataAccessor
    void setDataAccessor(DataAccessor* owner);

    /// resize status vectors for number of symmetrizations
    void setNumberOfSymmetrizations(unsigned n);

    /// resize status vectors for number of data partitions
    void setNumberOfDataPartitions(unsigned n);

private:

    DataAccessor* Owner_;       ///< Owning #DataAccessor
    int Position_;              ///< Position of first element of cached value within data vector
    unsigned Size_;             ///< Size of cached value (number of real elements)

    ParameterSet ParametersItDependsOn_;
    CachedDataValueSet CachedDataValuesItDependsOn_;
    CachedDataValueDaughterIndexPairSet DaughterCachedDataValuesItDependsOn_;

    /// CalculationStatus'es for the current DataPoint
    /// first index is for data partion
    /// second index is for symmetrization
    std::vector<std::vector<CalculationStatus> > CalculationStatus_;

    /// index is for symmetrization
    std::vector<CalculationStatus> GlobalCalculationStatus_;

    /// first index is for data partion
    /// second index is for symmetrization
    std::vector<std::vector<VariableStatus> > VariableStatus_;

};

/// \class RealCachedDataValue
/// \brief Class for managing a single real cached value inside a #DataPoint
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
/// \ingroup Cache
class RealCachedDataValue : public CachedDataValue
{
public:

    /// create shared_ptr to RealCachedDataValue
    /// \param owner #DataAccessor to which this cached value belongs
    /// \param pars set of shared pointers to Parameters cached value depends on
    /// \param vals set of shared pointers to CachedValues cached value depends on
    static std::shared_ptr<RealCachedDataValue> create(DataAccessor* da, ParameterSet pars = {}, CachedDataValueSet vals = {});

    /// Set value into #DataPoint for particular symmetrization, and
    /// update VariableStatus for symm. and partition index
    /// \param val Value to set to
    /// \param d #DataPoint to update
    /// \param symmetrizationIndex index of symmetrization to apply to
    /// \param dataPartitionIndex index of data partition being worked on
    void setValue(double val, DataPoint& d, unsigned symmetrizationIndex, unsigned dataPartitionIndex);

    /// Get value from #DataPoint for particular symmetrization
    /// \param d #DataPoint to get value from
    /// \param symmetrizationIndex index of symmetrization to grab from
    /// \return Value of CachedDataValue inside the data point
    double value(const DataPoint& d, unsigned symmetrizationIndex) const
    { return CachedDataValue::value(0, d, symmetrizationIndex); }

private:

    /// Constructor
    RealCachedDataValue(ParameterSet pars = {}, CachedDataValueSet vals = {})
        : CachedDataValue(1, pars, vals) {}

};

/// \class ComplexCachedDataValue
/// \brief Class for managing a complex cached value inside a #DataPoint
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
/// \ingroup Cache
class ComplexCachedDataValue : public CachedDataValue
{
public:

    /// create shared pointer to ComplexCachedDataValue
    /// \param owner #DataAccessor to which this cached value belongs
    /// \param pars set of shared pointers to Parameters cached value depends on
    /// \param vals set of shared pointers to CachedValues cached value depends on
    static std::shared_ptr<ComplexCachedDataValue> create(DataAccessor* da, ParameterSet pars = {}, CachedDataValueSet vals = {});

    /// Set value into #DataPoint for particular symmetrization, and
    /// update VariableStatus for symm. and partition index
    /// \param val Value to set to
    /// \param d #DataPoint to update
    /// \param symmetrizationIndex index of symmetrization to apply to
    /// \param dataPartitionIndex index of data partition being worked on
    void setValue(std::complex<double> val, DataPoint& d, unsigned symmetrizationIndex, unsigned dataPartitionIndex)
    { setValue(real(val), imag(val), d, symmetrizationIndex, dataPartitionIndex); }

    /// Set value into #DataPoint for particular symmetrization, and
    /// update VariableStatus for symm. and partition index
    /// \param val_re real part of value to set to
    /// \param val_im imaginary part of value to set to
    /// \param d #DataPoint to update
    /// \param symmetrizationIndex index of symmetrization to apply to
    /// \param dataPartitionIndex index of data partition being worked on
    void setValue(double val_re, double val_im, DataPoint& d, unsigned symmetrizationIndex, unsigned dataPartitionIndex);

    /// Get value from #DataPoint for particular symmetrization
    /// \param d #DataPoint to get value from
    /// \param symmetrizationIndex index of symmetrization to grab from
    /// \return Value of CachedDataValue inside the data point
    std::complex<double> value(const DataPoint& d, unsigned  symmetrizationIndex) const
    { return std::complex<double>(CachedDataValue::value(0, d, symmetrizationIndex), CachedDataValue::value(1, d, symmetrizationIndex)); }

private:

    /// Constructor (protected)
    /// see #create for details
    ComplexCachedDataValue(ParameterSet pars = {}, CachedDataValueSet vals = {})
        : CachedDataValue(2, pars, vals) {}


};

/// \class FourVectorCachedDataValue
/// \brief Class for managing a four-vector cached value inside a #DataPoint
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
/// \ingroup Cache
class FourVectorCachedDataValue : public CachedDataValue
{
public:

    /// create shared pointer to ComplexCachedDataValue
    /// \param owner #DataAccessor to which this cached value belongs
    /// \param pars set of shared pointers to Parameters cached value depends on
    /// \param vals set of shared pointers to CachedValues cached value depends on
    static std::shared_ptr<FourVectorCachedDataValue> create(DataAccessor* da, ParameterSet pars = {}, CachedDataValueSet vals = {});

    /// Set value into #DataPoint for particular symmetrization, and
    /// update VariableStatus for symm. and partition index
    /// \param val Value to set to
    /// \param d #DataPoint to update
    /// \param symmetrizationIndex index of symmetrization to apply to
    /// \param dataPartitionIndex index of data partition being worked on
    void setValue(FourVector<double> val, DataPoint& d, unsigned symmetrizationIndex, unsigned dataPartitionIndex);

    /// Get value from #DataPoint for particular symmetrization
    /// \param d #DataPoint to get value from
    /// \param symmetrizationIndex index of symmetrization to grab from
    /// \return Value of CachedDataValue inside the data point
    FourVector<double> value(const DataPoint& d, unsigned  symmetrizationIndex) const
    {
        return FourVector<double>( { CachedDataValue::value(0, d, symmetrizationIndex),
                                     CachedDataValue::value(1, d, symmetrizationIndex),
                                     CachedDataValue::value(2, d, symmetrizationIndex),
                                     CachedDataValue::value(3, d, symmetrizationIndex)
                                   });
    }

private:

    /// Constructor (protected)
    /// see #create for details
    FourVectorCachedDataValue(ParameterSet pars = {}, CachedDataValueSet vals = {})
        : CachedDataValue(4, pars, vals) {}


};

}

#endif

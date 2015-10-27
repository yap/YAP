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
#include "DataAccessor.h"
#include "DataPoint.h"
#include "Parameter.h"
#include "ParticleCombination.h"
#include "VariableStatus.h"

#include <memory>
#include <set>
#include <utility>
#include <vector>

namespace yap {

class CachedDataValue;

/// \typedef CachedDataValueSet
/// \ingroup Data
/// \ingroup Cache
using CachedDataValueSet = std::set<std::shared_ptr<CachedDataValue> >;

/// \class CachedDataValue
/// \brief Class for managing cached values inside a #DataPoint
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
/// \ingroup Cache

class CachedDataValue
{
public:
    /// Constructor
    /// \param owner #DataAccessor to which this cached value belongs
    /// \param size Length of cached value (number of real elements)
    /// \param pars set of shared pointers to Parameters cached value depends on
    /// \param vals set of shared pointers to CachedValues cached value depends on
    CachedDataValue(DataAccessor* owner, unsigned size, ParameterSet pars = {}, CachedDataValueSet vals = {});

    /// \name CachedDataValue friends
    /// @{

    friend class DataAccessor;
    friend class InitialStateParticle;

    /// @}

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

    /// add CachedDataValue this CachedDataValue depends on
    void addDependency(std::shared_ptr<CachedDataValue> dep)
    { CachedDataValuesItDependsOn_.insert(dep); }

    /// remove dependency
    void removeDependency(std::shared_ptr<ParameterBase> dep);

    /// @}

    /// \name Getters
    /// @{

    DataAccessor* owner() const
    { return Owner_; }

    /// overload and hide #CachedValue::calculationStatus
    /// \return #CalculationStatus of symmetrization index and data-partition index
    /// \param pc shared pointer to #ParticleCombination to check status of
    /// \param symmetrizationIndex index of symmetrization to check status of
    /// \param dataPartitionIndex index of dataPartitionIndex to check status of
    CalculationStatus calculationStatus(const std::shared_ptr<const ParticleCombination>& pc, unsigned symmetrizationIndex, unsigned dataPartitionIndex)
    { return CalculationStatus_[dataPartitionIndex][symmetrizationIndex]; }

    /// overload and hide #CachedValue::calculationStatus
    /// \return #CalculationStatus of symmetrization index and data-partition index
    /// \param pc_symInd pair of shared pointer to #ParticleCombination and symmetrization index
    /// \param dataPartitionIndex index of dataPartitionIndex to check status of
    CalculationStatus calculationStatus(std::pair<std::shared_ptr<const ParticleCombination>, unsigned> pc_symInd, unsigned dataPartitionIndex)
    { return calculationStatus(pc_symInd.first, pc_symInd.second, dataPartitionIndex); }

    /// overload and hide #CachedValue::calculationStatus
    /// \return #CalculationStatus of symmetrization index and data-partition index
    /// \param pc shared pointer to #ParticleCombination to check status of
    /// \param dataPartitionIndex index of dataPartitionIndex to check status of
    CalculationStatus calculationStatus(const std::shared_ptr<const ParticleCombination>& pc, unsigned dataPartitionIndex)
    { return calculationStatus(pc, Owner_->symmetrizationIndex(pc), dataPartitionIndex); }


    /// get global CalculationStatus
    /// \return #CalculationStatus of symmetrization index and data-partition index
    /// \param symmetrization index
    CalculationStatus globalCalculationStatus(unsigned symmetrizationIndex)
    { return GlobalCalculationStatus_[symmetrizationIndex]; }

    /// get global CalculationStatus
    /// \return #CalculationStatus of symmetrization index and data-partition index
    /// \param pc shared pointer to #ParticleCombination to check status of
    CalculationStatus globalCalculationStatus(const std::shared_ptr<const ParticleCombination>& pc)
    { return GlobalCalculationStatus_[Owner_->symmetrizationIndex(pc)]; }


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
    double value(unsigned index, const DataPoint& d, unsigned symmetrizationIndex) const
    { return d.Data_[Owner_->index()][symmetrizationIndex][Position_ + index];}

    /// \return Size of cached value (number of real elements)
    virtual unsigned size() const
    { return Size_; }

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
    void setCalculationStatus(CalculationStatus stat, const std::shared_ptr<const ParticleCombination>& pc, unsigned dataPartitionIndex)
    { CalculationStatus_[dataPartitionIndex][Owner_->symmetrizationIndex(pc)] = stat; }

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
    void setGlobalCalculationStatus(CalculationStatus stat, const std::shared_ptr<const ParticleCombination>& pc)
    { GlobalCalculationStatus_[Owner_->symmetrizationIndex(pc)] = stat; }


    /// Set value into #DataPoint for particular symmetrization
    /// (No update to VariableStatus or CalculationStatus is made!)
    /// \param index index of value to get from within cached value (must be less than #Size_)
    /// \param val Value to set to
    /// \param d #DataPoint to update
    /// \param symmetrizationIndex index of symmetrization to apply to
    void setValue(unsigned index, double val, DataPoint& d, unsigned symmetrizationIndex) const
    { d.Data_[Owner_->index()][symmetrizationIndex][Position_ + index] = val; }

    /// @}

    void updateGlobalCalculationStatus(const std::shared_ptr<const ParticleCombination>& pc);
    void resetCalculationStatus(unsigned dataPartitionIndex);

protected:
    DataAccessor* Owner_;       ///< Owning #DataAccessor
    int Position_;              ///< Position of first element of cached value within data vector
    unsigned Size_;             ///< Size of cached value (number of real elements)

    ParameterSet ParametersItDependsOn_;
    CachedDataValueSet CachedDataValuesItDependsOn_;

    /// CalculationStatus'es for the current DataPoint
    /// first index is for data partion
    /// second index is for symmetrization
    std::vector<std::vector<CalculationStatus> > CalculationStatus_;

    /// index is for symmetrization
    std::vector<CalculationStatus> GlobalCalculationStatus_;

    /// first index is for data partion
    /// second index is for symmetrization
    std::vector<std::vector<VariableStatus> > VariableStatus_;

    /// resize status vectors for number of symmetrizations
    void setNumberOfSymmetrizations(unsigned n);

    /// resize status vectors for number of data partitions
    void setNumberOfDataPartitions(unsigned n);

};

/// \class RealCachedValue
/// \brief Class for managing a single real cached value inside a #DataPoint
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
/// \ingroup Cache

class RealCachedDataValue : public CachedDataValue
{
public:

    /// Constructor
    /// \param owner #DataAccessor to which this cached value belongs
    /// \param pars set of shared pointers to Parameters cached value depends on
    /// \param vals set of shared pointers to CachedValues cached value depends on
    RealCachedDataValue(DataAccessor* owner, ParameterSet pars = {}, CachedDataValueSet vals = {})
        : CachedDataValue(owner, 1, pars, vals)
    {}

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
    double value(const DataPoint& d, unsigned  symmetrizationIndex) const
    { return CachedDataValue::value(0, d, symmetrizationIndex); }
};

/// \class ComplexCachedValue
/// \brief Class for managing a complex cached value inside a #DataPoint
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
/// \ingroup Cache

class ComplexCachedDataValue : public CachedDataValue
{
public:

    /// Constructor
    /// \param owner #DataAccessor to which this cached value belongs
    /// \param pars set of shared pointers to Parameters cached value depends on
    /// \param vals set of shared pointers to CachedValues cached value depends on
    ComplexCachedDataValue(DataAccessor* owner, ParameterSet pars = {}, CachedDataValueSet vals = {})
        : CachedDataValue(owner, 2, pars, vals)
    {}

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
};

}

#endif

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
#include "Parameter.h"
#include "VariableStatus.h"

#include <memory>
#include <set>
#include <vector>

namespace yap {

/// \class CachedDataValueBase
/// \brief Class for managing cached values inside a #DataPoint
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data

class CachedDataValueBase
{
public:
    /// Constructor
    /// \param ParametersItDependsOn vector of shared pointers to Parameters cached value depends on
    /// \param CachedValuesItDependsOn vector of shared pointers to CachedValues cached value depends on
    CachedDataValueBase(std::vector<std::shared_ptr<Parameter> > ParametersItDependsOn = {},
                        std::vector<std::shared_ptr<CachedDataValueBase> > CachedDataValuesItDependsOn = {});

    /// add Parameters this CachedDataValueBase depends on
    void addDependencies(std::vector<std::shared_ptr<Parameter> > deps)
    { for (auto& dep : deps) addDependency(dep); }

    /// add Parameter this CachedDataValueBase depends on
    void addDependency(std::shared_ptr<Parameter> dep)
    { ParametersItDependsOn_.insert(dep); }

    /// add CachedDataValueBase's this CachedDataValueBase depends on
    void addDependencies(std::vector<std::shared_ptr<CachedDataValueBase> > deps)
    { for (auto& dep : deps) addDependency(dep); }

    /// add CachedDataValueBase this CachedDataValueBase depends on
    void addDependency(std::shared_ptr<CachedDataValueBase> dep)
    { CachedDataValuesItDependsOn_.insert(dep); }

    /// overload and hide #CachedValueBase::calculationStatus
    /// \return #CalculationStatus of symmetrization index and data-partition index
    /// \param symmetrizationIndex index of symmetrization to check status of
    /// \param dataPartitionIndex index of dataPartitionIndex to check status of
    /// \todo Use ParticleCombination instead of symmetrization index?
    /// This would allow us to use dependencies from outside the
    /// owning DataAccessor.
    CalculationStatus calculationStatus(unsigned symmetrizationIndex, unsigned dataPartitionIndex);

    /// \return VariableStatus for symmetrization index and data-partition index
    /// \param symmetrizationIndex index of symmetrization to check status of
    /// \param dataPartitionIndex index of dataPartitionIndex to check status of
    VariableStatus variableStatus(unsigned symmetrizationIndex, unsigned dataPartitionIndex) const
    { return VariableStatus_[symmetrizationIndex][dataPartitionIndex]; }

    /// set VariableStatus for symmetrization index and data-partition index
    /// \param symmetrizationIndex index of symmetrization to check status of
    /// \param dataPartitionIndex index of dataPartitionIndex to check status of
    /// \param stat VariableStatus to set to
    void setVariableStatus(unsigned symmetrizationIndex, unsigned dataPartitionIndex, VariableStatus stat)
    { VariableStatus_[symmetrizationIndex][dataPartitionIndex] = stat; }

protected:
    std::set<std::shared_ptr<Parameter> > ParametersItDependsOn_;
    std::set<std::shared_ptr<CachedDataValueBase> > CachedDataValuesItDependsOn_;

    /// first index is for data partion
    /// second index is for symmetrization
    std::vector<std::vector<CalculationStatus> > CalculationStatus_;

    /// first index is for data partion
    /// second index is for symmetrization
    std::vector<std::vector<VariableStatus> > VariableStatus_;

};

// /// \class RealCachedValue
// /// \brief extension of #CachedValue for a real numbers
// /// \author Johannes Rauch, Daniel Greenwald
// /// \ingroup Parameters

// class RealCachedValue : public CachedValue
// {
// public:

//     /// Constructor
//     /// \param ParametersItDependsOn vector of shared pointers to Parameters cached value depends on
//     /// \param CachedValuesItDependsOn vector of shared pointers to CachedValues cached value depends on
//     RealCachedValue(std::vector<std::shared_ptr<Parameter> > ParametersItDependsOn = {},
//                     std::vector<std::shared_ptr<CachedValueBase> > CachedValuesItDependsOn = {})
//         : CachedValue(ParametersItDependsOn, CachedValuesItDependsOn)
//     {}

//     /// Replaces & hides #CachedValue::value with return type double
//     const double value() const
//     { return real(CachedValue_); }

//     /// get number of real components in cached value (real number = 1 real components)
//     const unsigned size() const
//     { return 1; }

//     /// Overloading & hides #CachedValue::setValue with argument double
//     void setValue(double val)
//     { CachedValue::setValue(val * Complex_1); }
// };

}

#endif

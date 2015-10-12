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

#ifndef yap_CachedValue_h
#define yap_CachedValue_h

#include "Parameter.h"

#include <memory>

namespace yap {

class CachedValueBase {
public:
    /// Constructor
    CachedValueBase(std::vector<std::shared_ptr<ParameterBase> > ParametersItDependsOn,
        std::vector<std::shared_ptr<CachedValueBase> > CachedValuesItDependsOn) :
      ParametersItDependsOn_(ParametersItDependsOn),
      CachedValuesItDependsOn_(CachedValuesItDependsOn),
      VariableStatus_(kChanged)
      CalculationStatus_(kUncalculated),
    {}

    /// \name getters
    /// @{

    /// get CalculationStatus of ith DataPartition
    CalculationStatus calculationStatus(unsigned iDataPartition) const
    { return CalculationStatuses_.at(iDataPartition); }

    VariableStatus variableStatus() const
    { return VariableStatus_; }

    /// @}

    /// \name setters
    /// @{

    /// get CalculationStatus of ith DataPartition
    void setCalculationStatus(CalculationStatus stat)
    { CalculationStatus_ = stat; }

    void setVariableStatus(VariableStatus stat)
    { VariableStatus_ = stat; }

    /// @}

    CalculationStatus cache()
    {
        result = CalculationStatus_;
        for (auto& p : ParametersItDependsOn_) {
            // \todo check p
        }

        for (auto& v : CachedValuesItDependsOn_) {
            if (v.cache() == kUncalculated)
                result = kUncalculated;
        }

        return result;
    }


private:
    std::vector<std::shared_ptr<ParameterBase> > ParametersItDependsOn_;
    std::vector<std::shared_ptr<CachedValueBase> > CachedValuesItDependsOn_;
    VariableStatus VariableStatus_;
    CalculationStatus CalculationStatus_;
};


template <typename T>
class CachedValue : public CachedValueBase {
public:
    /// Constructor
    CachedValue(std::vector<std::shared_ptr<ParameterBase> > ParametersItDependsOn,
          std::vector<std::shared_ptr<CachedValueBase> > CachedValuesItDependsOn) :
      CachedValueBase(ParametersItDependsOn, CachedValuesItDependsOn)
    {}

    /// \name getters
    /// @{

    /// get CalculationStatus of ith DataPartition
    const T& value() const
    { return CachedValueValue_; }

    /// @}

    /// \name setters
    /// @{

    void setValue(T val)
    { CachedValueValue_ = val; }

    /// @}


private:
    T CachedValueValue_;
};

}

#endif

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

class CachedValue
{
public:
    /// Constructor
    CachedValue(std::vector<std::shared_ptr<ParameterBase> > ParametersItDependsOn = {},
                    std::vector<std::shared_ptr<CachedValueBase> > CachedValuesItDependsOn = {}) :
        ParametersItDependsOn_(ParametersItDependsOn),
        CachedValuesItDependsOn_(CachedValuesItDependsOn),
        CalculationStatus_(kUncalculated)
    {}

    CachedValue(ParameterSet ParametersItDependsOn = {},
                    std::vector<std::shared_ptr<CachedValueBase> > CachedValuesItDependsOn = {}) :
        ParametersItDependsOn_(ParametersItDependsOn),
        CachedValuesItDependsOn_(CachedValuesItDependsOn),
        CalculationStatus_(kUncalculated)
    {}

    /// \name getters
    /// @{

    /// get CalculationStatus of ith DataPartition
    const T& value() const
    { return CachedValueValue_; }

    /// get CalculationStatus of ith DataPartition
    CalculationStatus calculationStatus(unsigned iDataPartition) const
    { return CalculationStatuses_.at(iDataPartition); }

    VariableStatus variableStatus() const
    { return VariableStatus_; }

    /// @}

    /// \name setters
    /// @{

    void setValue(std::complex<double> val)
    {
        if (val != CachedValue_) {
            CachedValue_ = val;
            CalculationStatus_ = kCalculated;
        }
    }

    /// get CalculationStatus of ith DataPartition
    void setCalculationStatus(CalculationStatus stat)
    { CalculationStatus_ = stat; }

    /// @}

    CalculationStatus cache()
    {
        result = CalculationStatus_;
        for (auto& p : ParametersItDependsOn_) {
            if (p->variableStatus() == kChanged)
                result = kUncalculated;
        }

        for (auto& v : CachedValuesItDependsOn_) {
            if (v.cache() == kUncalculated)
                result = kUncalculated;
        }

        return result;
    }

private:
    std::complex<double> CachedValue_;
    std::vector<std::shared_ptr<ParameterBase> > ParametersItDependsOn_;
    std::vector<std::shared_ptr<CachedValueBase> > CachedValuesItDependsOn_;
    CalculationStatus CalculationStatus_;
};

}

#endif

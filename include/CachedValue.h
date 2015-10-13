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
#include "ParameterSet.h"

#include <memory>

namespace yap {

class CachedValue
{
public:
    /// Constructor
    CachedValue(std::vector<std::shared_ptr<Parameter> > ParametersItDependsOn = {},
                    std::vector<std::shared_ptr<CachedValue> > CachedValuesItDependsOn = {}) :
        ParametersItDependsOn_(ParametersItDependsOn),
        CachedValuesItDependsOn_(CachedValuesItDependsOn),
        CalculationStatus_(kUncalculated)
    {}

    /// \name getters
    /// @{

    /// get
    const std::complex<double>& value() const
    { return CachedValue_; }

    /// get CalculationStatus of ith DataPartition
    CalculationStatus calculationStatus() const
    { return CalculationStatus_; }

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

    void addDependencies(std::vector<std::shared_ptr<Parameter> > dep)
    { ParametersItDependsOn_.insert(ParametersItDependsOn_.end(), dep.begin(), dep.end()); }

    void addDependencies(std::shared_ptr<Parameter> dep)
    { ParametersItDependsOn_.push_back(dep); }

    void addDependencies(std::vector<std::shared_ptr<CachedValue> > dep)
    { CachedValuesItDependsOn_.insert(CachedValuesItDependsOn_.end(), dep.begin(), dep.end()); }

    void addDependencies(std::shared_ptr<CachedValue> dep)
    { CachedValuesItDependsOn_.push_back(dep); }

    CalculationStatus cache()
    {
        for (auto& p : ParametersItDependsOn_) {
            if (p->variableStatus() == kChanged)
                CalculationStatus_ = kUncalculated;
        }

        for (auto& v : CachedValuesItDependsOn_) {
            if (v->cache() == kUncalculated)
                CalculationStatus_ = kUncalculated;
        }

        return CalculationStatus_;
    }

private:
    std::complex<double> CachedValue_;
    std::vector<std::shared_ptr<Parameter> > ParametersItDependsOn_;
    std::vector<std::shared_ptr<CachedValue> > CachedValuesItDependsOn_;
    CalculationStatus CalculationStatus_;
};

}

#endif

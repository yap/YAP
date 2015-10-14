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

#include "Constants.h"
#include "Parameter.h"
#include "ParameterSet.h"

#include <memory>
#include <set>

namespace yap {

/// \class CachedValue
/// \brief Class for managing cached values (as complex numbers)
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Parameters

class CachedValue
{
public:
    /// Constructor
    /// \param ParametersItDependsOn vector of shared pointers to Parameters cached value depends on
    /// \param CachedValuesItDependsOn vector of shared pointers to CachedValues cached value depends on
    CachedValue(std::vector<std::shared_ptr<Parameter> > ParametersItDependsOn = {},
                std::vector<std::shared_ptr<CachedValue> > CachedValuesItDependsOn = {});

    /// \name getters
    /// @{

    /// get value
    const std::complex<double>& value() const
    { return CachedValue_; }

    /// get number of real components in cached value (complex number = 2 real components)
    const unsigned size() const
    { return 2; }

    /// @}

    /// \name setters
    /// @{

    /// set value and set CalculationStatus_ to kCalculated
    void setValue(std::complex<double> val);

    /// @}

    /// add Parameters this CachedValue depends on
    void addDependencies(std::vector<std::shared_ptr<Parameter> > dep);

    /// add Parameter this CachedValue depends on
    void addDependency(std::shared_ptr<Parameter> dep)
    { ParametersItDependsOn_.insert(dep); }

    /// add CachedValues this CachedValue depends on
    void addDependencies(std::vector<std::shared_ptr<CachedValue> > dep);

    /// add CachedValue this CachedValue depends on
    void addDependency(std::shared_ptr<CachedValue> dep)
    { CachedValuesItDependsOn_.insert(dep); }

    /// update (depending on Parameters and CachedValues it depends) and return CalculationStatus_
    CalculationStatus calculationStatus();

    /// set VariableStatus of ParametersItDependsOn_ to kUnchanged (or leave at kFixed)
    void finishedPrecalculation();

protected:
    std::complex<double> CachedValue_;
    std::set<std::shared_ptr<Parameter> > ParametersItDependsOn_;
    std::set<std::shared_ptr<CachedValue> > CachedValuesItDependsOn_;
    CalculationStatus CalculationStatus_;
};

/// \class RealCachedValue
/// \brief extension of #CachedValue for a real numbers
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Parameters

class RealCachedValue : public CachedValue
{
public:

    /// Constructor
    /// \param ParametersItDependsOn vector of shared pointers to Parameters cached value depends on
    /// \param CachedValuesItDependsOn vector of shared pointers to CachedValues cached value depends on
    RealCachedValue(std::vector<std::shared_ptr<Parameter> > ParametersItDependsOn = {},
                    std::vector<std::shared_ptr<CachedValue> > CachedValuesItDependsOn = {})
        : CachedValue(ParametersItDependsOn, CachedValuesItDependsOn)
    {}

    /// Replaces & hides #CachedValue::value with return type double
    const double value() const
    { return real(CachedValue_); }

    /// get number of real components in cached value (real number = 1 real components)
    const unsigned size() const
    { return 1; }

    /// Overloading & hides #CachedValue::setValue with argument double
    void setValue(double val)
    { CachedValue::setValue(val * Complex_1); }
};

}

#endif

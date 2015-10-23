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

#include "CalculationStatus.h"
#include "Constants.h"
#include "Parameter.h"
#include "VariableStatus.h"

#include <memory>
#include <set>
#include <vector>

namespace yap {

/// \class CachedValueBase
/// \brief Base class for cached value managers
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Parameters

class CachedValueBase
{
public:
    /// Constructor
    /// \param ParametersItDependsOn set of shared pointers to Parameters cached value depends on
    CachedValueBase(ParameterSet ParametersItDependsOn = {});

    /// add Parameters this CachedValueBase depends on
    void addDependencies(ParameterSet deps)
    { for (auto& dep : deps) addDependency(dep); }

    /// add Parameter this CachedValueBase depends on
    void addDependency(std::shared_ptr<ParameterBase> dep)
    { ParametersItDependsOn_.insert(dep); }

    /// remove Parameter from list of dependencies
    void removeDependency(std::shared_ptr<ParameterBase> dep);

    /// update (depending on Parameters and CachedValueBase's it
    /// depends) and return CalculationStatus_
    CalculationStatus calculationStatus();

    /// set VariableStatus of members of #ParametersItDependsOn_ to
    /// kUnchanged (or leave at kFixed).
    void finishedPrecalculation();

protected:
    ParameterSet ParametersItDependsOn_;
    CalculationStatus CalculationStatus_;

};

/// \class ComplexCachedValue
/// \brief Class for managing cached values (as complex numbers)
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Parameters

class ComplexCachedValue : public CachedValueBase, public ComplexParameter
{
public:
    /// Constructor
    /// \param ParametersItDependsOn set of shared pointers to Parameters cached value depends on
    ComplexCachedValue(ParameterSet ParametersItDependsOn = {})
        : CachedValueBase(ParametersItDependsOn), ComplexParameter()
    {}

    /// set value and set CalculationStatus_ to kCalculated
    void setValue(std::complex<double> val)
    { ComplexParameter::setValue(val); CalculationStatus_ = kCalculated; }

};

/// \class RealCachedValue
/// \brief extension of #ComplexCachedValue for a real numbers
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Parameters

class RealCachedValue : public CachedValueBase, public RealParameter
{
public:

    /// Constructor
    /// \param ParametersItDependsOn set of shared pointers to Parameters cached value depends on
    RealCachedValue(ParameterSet ParametersItDependsOn = {})
        : CachedValueBase(ParametersItDependsOn), RealParameter()
    {}

    /// Overloading & hides #ComplexCachedValue::setValue with argument double
    void setValue(double val)
    { RealParameter::setValue(val); CalculationStatus_ = kCalculated; }

};

}

#endif

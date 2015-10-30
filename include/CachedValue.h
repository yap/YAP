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
/// \defgroup Cache

class CachedValueBase
{
public:

    /// Constructor
    /// \param pars set of shared pointers to Parameters cached value depends on
    CachedValueBase(ParameterSet pars = {});

    /// add Parameters this CachedValueBase depends on
    void addDependencies(ParameterSet deps)
    { for (auto& dep : deps) addDependency(dep); }

    /// add Parameter this CachedValueBase depends on
    void addDependency(std::shared_ptr<ParameterBase> dep)
    { ParametersItDependsOn_.insert(dep); }

    /// remove Parameter from list of dependencies
    void removeDependency(std::shared_ptr<ParameterBase> dep);

    /// remove dependencies
    void removeDependencies(ParameterSet deps)
    { for (auto& dep : deps) removeDependency(dep); }

    /// update (depending on Parameters and CachedValueBase's it
    /// depends) and return CalculationStatus_
    CalculationStatus calculationStatus();

protected:

    ParameterSet ParametersItDependsOn_;
    CalculationStatus CalculationStatus_;

};

/// \typedef CachedValueSet
/// \ingroup Parameters
/// \ingroup Cache
using CachedValueSet = std::set<std::shared_ptr<CachedValueBase> >;

/// \class CachedValue
/// \brief Class for linking CachedValueBase and Parameter
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Parameters
/// \ingroup Cache

template <typename T>
class CachedValue : public CachedValueBase, public Parameter<T>
{
public:

    /// Constructor
    /// \param pars set of shared pointers to Parameters cached value depends on
    CachedValue(ParameterSet pars = {})
        : CachedValueBase(pars), Parameter<T>()
    {}

    /// set value and set CalculationStatus_ to kCalculated
    void setValue(T val)
    { Parameter<T>::setValue(val); CalculationStatus_ = kCalculated; }

};

/// \typedef ComplexCachedValue
/// \ingroup Parameters
/// \ingroup Cache
using ComplexCachedValue = CachedValue<std::complex<double> >;

/// \typedef RealCachedValue
/// \ingroup Parameters
/// \ingroup Cache
using RealCachedValue = CachedValue<double>;

}

#endif

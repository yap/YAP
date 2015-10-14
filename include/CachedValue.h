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

#include "CachedValueBase.h"
#include "Constants.h"
#include "Parameter.h"

#include <memory>

namespace yap {

/// \class CachedValue
/// \brief Class for managing cached values (as complex numbers)
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Parameters

class CachedValue : public CachedValueBase, public Parameter
{
public:
    /// Constructor
    /// \param ParametersItDependsOn vector of shared pointers to Parameters cached value depends on
    /// \param CachedValuesItDependsOn vector of shared pointers to CachedValues cached value depends on
    CachedValue(std::vector<std::shared_ptr<Parameter> > ParametersItDependsOn = {})
        : CachedValueBase(ParametersItDependsOn), Parameter()
    {}

    /// set value and set CalculationStatus_ to kCalculated
    void setValue(std::complex<double> val)
    { Parameter::setValue(val); CalculationStatus_ = kCalculated; }

};

/// \class RealCachedValue
/// \brief extension of #CachedValue for a real numbers
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Parameters

class RealCachedValue : public CachedValueBase, public RealParameter
{
public:

    /// Constructor
    /// \param ParametersItDependsOn vector of shared pointers to Parameters cached value depends on
    RealCachedValue(std::vector<std::shared_ptr<Parameter> > ParametersItDependsOn = {})
        : CachedValueBase(ParametersItDependsOn), RealParameter()
    {}

    /// Overloading & hides #CachedValue::setValue with argument double
    void setValue(double val)
    { RealParameter::setValue(val); CalculationStatus_ = kCalculated; }

};

}

#endif

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

#ifndef yap_CachedValueBase_h
#define yap_CachedValueBase_h

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
    /// \param ParametersItDependsOn vector of shared pointers to Parameters cached value depends on
    CachedValueBase(std::vector<std::shared_ptr<Parameter> > ParametersItDependsOn = {});

    /// add Parameters this CachedValueBase depends on
    void addDependencies(std::vector<std::shared_ptr<Parameter> > deps)
    { for (auto& dep : deps) addDependency(dep); }

    /// add Parameter this CachedValueBase depends on
    void addDependency(std::shared_ptr<Parameter> dep)
    { ParametersItDependsOn_.insert(dep); }

    /// update (depending on Parameters and CachedValueBase's it
    /// depends) and return CalculationStatus_
    CalculationStatus calculationStatus();

    /// set VariableStatus of members of #ParametersItDependsOn_ to
    /// kUnchanged (or leave at kFixed).
    void finishedPrecalculation();

protected:
    std::set<std::shared_ptr<Parameter> > ParametersItDependsOn_;
    CalculationStatus CalculationStatus_;

};

}

#endif

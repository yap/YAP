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

#ifndef yap_Parameter_h
#define yap_Parameter_h

#include "CalculationStatus.h"
#include "VariableStatus.h"

#include <complex>
#include <vector>

namespace yap {

class Parameter
{
public:
    /// Constructor
    Parameter(unsigned nDataPartitions = 1) :
        VariableStatus_(kChanged),
        CalculationStatuses_(nDataPartitions, kUncalculated)
    {}

    /// \name getters
    /// @{

    const std::complex<double>& value() const
    { return ParameterValue_; }

    /// get CalculationStatus of ith DataPartition
    CalculationStatus calculationStatus(unsigned iDataPartition) const
    { return CalculationStatuses_.at(iDataPartition); }

    VariableStatus variableStatus() const
    { return VariableStatus_; }

    /// @}

    /// \name setters
    /// @{

    void setValue(std::complex<double> val)
    { ParameterValue_ = val; }

    /// get CalculationStatus of ith DataPartition
    void setCalculationStatus(unsigned iDataPartition, CalculationStatus stat)
    { CalculationStatuses_.at(iDataPartition) = stat; }

    void setVariableStatus(VariableStatus stat)
    { VariableStatus_ = stat; }

    /// @}

private:
    std::complex<double> ParameterValue_;
    VariableStatus VariableStatus_;
    /// One CalculationStatus per DataPartition
    std::vector<CalculationStatus> CalculationStatuses_;
};

}

#endif

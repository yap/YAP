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

#include "logging.h"
#include "VariableStatus.h"

#include <complex>
#include <memory>
#include <set>
#include <vector>

namespace yap {

/// \class ParameterBase
/// \brief Class holding basic properties of a parameter, but not a value!
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup Parameters

class ParameterBase
{
public:

    // \param size Number of real components in variable
    ParameterBase() : VariableStatus_(kChanged)
    {}

    /// \return VariableStatus
    VariableStatus variableStatus() const
    { return VariableStatus_; }

    // /// \return number of real parameters
    // virtual unsigned size() const
    // { return Size_; }

    /// set VariableStatus
    void setVariableStatus(VariableStatus stat)
    { VariableStatus_ = stat; }

protected:

    // unsigned Size_;
    VariableStatus VariableStatus_;

};

/// \typedef ParameterVector
/// \ingroup Parameters
/// Use when order of parameters must be kept constant
using ParameterVector = std::vector<std::shared_ptr<ParameterBase> >;

/// \typedef ParameterSet
/// \ingroup Parameters
/// Use when enforcement of no duplicates is required, and order does not matter
using ParameterSet = std::set<std::shared_ptr<ParameterBase> >;

/// \class Parameter
/// \brief Template class holding also a value for a parameter
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Parameters

template <typename T>
class Parameter : public ParameterBase
{
public:

    /// Default constructor
    Parameter() : ParameterBase()
    { }

    /// Value-assigning constructor
    Parameter(T t) : ParameterBase(), ParameterValue_(t)
    {}

    /// \return value of parameter
    const T& value() const
    { return ParameterValue_; }

    /// set complex value
    void setValue(T val)
    {
        if (ParameterValue_ != val) {
            if (VariableStatus_ == kFixed) {
                FLOG(ERROR) << "Error, trying to change the value of a fixed parameter! Abort.";
                return;
            }
            ParameterValue_ = val;
            VariableStatus_ = kChanged;
        }
    }

protected:

    T ParameterValue_;

};

/// \typedef ComplexParameter
/// \ingroup Parameters
using ComplexParameter = Parameter<std::complex<double> >;

/// \typedef RealParameter
/// \ingroup Parameters
using RealParameter = Parameter<double>;

/// \typedef ComplexParameterVector
/// \ingroup Parameters
using ComplexParameterVector = std::vector<std::shared_ptr<ComplexParameter> >;

/// \typedef RealParameterVector
/// \ingroup Parameters
using RealParameterVector = std::vector<std::shared_ptr<RealParameter> >;

}

#endif

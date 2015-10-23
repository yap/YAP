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
#include "Constants.h"
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
    ParameterBase(/*unsigned size*/) : /*Size_(size),*/ VariableStatus_(kChanged)
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

/// \class Parameter
/// \brief Template class holding also a value for a parameter
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Parameters

template <typename T/*, unsigned size*/>
class Parameter : public ParameterBase
{
public:

    /// Default constructor
    Parameter() : ParameterBase()
    { }

    /// Value-assigning constructor
    Parameter(T t) : ParameterBase(/*size*/), ParameterValue_(t)
    {}

    /// \return value of parameter
    const T& value() const
    { return ParameterValue_; }

    /// set complex value
    void setValue(T val)
    {
        if (ParameterValue_ != val) {
            ParameterValue_ = val;
            VariableStatus_ = kChanged;
        }
    }

protected:
    T ParameterValue_;

};

/// \typedef ComplexParameter
using ComplexParameter = Parameter<std::complex<double>/*, 2u*/>;

/// \typedef RealParameter
using RealParameter = Parameter<double/*, 1u*/>;

/// \typedef ParameterVector
/// Use when order of parameters must be kept constant
using ParameterVector = std::vector<std::shared_ptr<ParameterBase> >;

/// \typedef ParameterSet
/// Use when enforcement of no duplicates is required, and order does not matter
using ParameterSet = std::set<std::shared_ptr<ParameterBase> >;

// /// \class ComplexParameter
// /// \brief Complex Parameter
// /// \author Johannes Rauch, Daniel Greenwald
// /// \ingroup Parameters

// class ComplexParameter
// {
// public:
//     /// Constructor
//     ComplexParameter(double real = 0, double imag = 0) :
//         ParameterValue_(real, imag),
//         VariableStatus_(kChanged)
//     {}

//     /// Constructor
//     ComplexParameter(std::complex<double> val) :
//         ParameterValue_(val),
//         VariableStatus_(kChanged)
//     {}

//     /// \name getters
//     /// @{

//     /// \return complex value
//     const std::complex<double>& value() const
//     { return ParameterValue_; }

//     /// \return VariableStatus
//     VariableStatus variableStatus() const
//     { return VariableStatus_; }

//     /// \return number of real parameters
//     virtual unsigned size() const
//     { return 2; }

//     /// @}

//     /// \name setters
//     /// @{

//     /// set complex value
//     void setValue(std::complex<double> val)
//     {
//         if (ParameterValue_ != val) {
//             ParameterValue_ = val;
//             VariableStatus_ = kChanged;
//         }
//     }

//     /// set VariableStatus
//     void setVariableStatus(VariableStatus stat)
//     { VariableStatus_ = stat; }

//     /// @}

// protected:
//     std::complex<double> ParameterValue_;
//     VariableStatus VariableStatus_;

// };

// /// \class RealParameter
// /// \brief Real Parameter
// /// \author Johannes Rauch, Daniel Greenwald
// /// \ingroup Parameters

// class RealParameter : public ComplexParameter
// {
// public:

//     /// Constructor
//     RealParameter(double real = 0) :
//         ComplexParameter(real * Complex_1)
//     {}

//     /// Replace & hide #Parameters::value returning double
//     double value() const
//     { return real(ParameterValue_); }

//     /// \return number of real components in parameter
//     unsigned size() const override
//     { return 1; }

//     /// Overloads & hides #Parameters::setValue taking double argument
//     void setValue(double val)
//     { ComplexParameter::setValue(val * Complex_1); }

// };

}

#endif

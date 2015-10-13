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

/// \class Parameter
/// \brief Complex Parameter
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup Parameters

class Parameter
{
public:
    /// Constructor
    Parameter(double real = 0, double imag = 0) :
        ParameterValue_(real, imag),
        VariableStatus_(kChanged)
    {}

    Parameter(std::complex<double> val) :
        ParameterValue_(val),
        VariableStatus_(kChanged)
    {}

    /// \name getters
    /// @{

    const std::complex<double>& value() const
    { return ParameterValue_; }

    VariableStatus variableStatus() const
    { return VariableStatus_; }

    /// @}

    /// \name setters
    /// @{

    virtual void setValue(std::complex<double> val)
    {
        if (ParameterValue_ != val) {
            ParameterValue_ = val;
            VariableStatus_ = kChanged;
        }
    }

    void setVariableStatus(VariableStatus stat)
    { VariableStatus_ = stat; }

    /// @}

protected:
    std::complex<double> ParameterValue_;
    VariableStatus VariableStatus_;

};

/// \class RealParameter
/// \brief Real Parameter
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Parameters

class RealParameter : public Parameter
{
public:

    RealParameter(double real = 0) :
        Parameter(real, 0)
    {}

    double realValue() const
    { return ParameterValue_.real(); }

    virtual void setValue(std::complex<double> val) override
    {
        val = (val.real()); // sets imag to 0
        Parameter::setValue(val);
    }

    void setValue(double val)
    { Parameter::setValue(std::complex<double>(val)); }

};

}

#endif

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

#include "Exceptions.h"
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
    ParameterBase() : VariableStatus_(VariableStatus::changed)
    {}

    /// \return VariableStatus
    VariableStatus variableStatus() const
    { return VariableStatus_; }

    /// set VariableStatus
    void setVariableStatus(VariableStatus stat)
    { VariableStatus_ = stat; }

private:

    /// Status of variable
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

    /// virtual destructor defaulted
    ~Parameter() = default;

    /// copy constructor defaulted
    Parameter(const Parameter&) = default;

    /// move constructor defaulted
    Parameter(Parameter&&) = default;

    /// copy assignment defaulted
    Parameter& operator=(const Parameter&) = default;

    /// move assignment defaulted
    Parameter& operator=(Parameter&&) = default;

    /// \return value of parameter
    virtual typename std::conditional<std::is_fundamental<T>::value, T, const T&>::type
    value() const
    { return ParameterValue_; }

    /// set value
    virtual void setValue(T val)
    {
        if (variableStatus() == VariableStatus::fixed)
            throw exceptions::ParameterIsFixed("", "Parameter::setValue");
        if (ParameterValue_ == val)
            return;
        ParameterValue_ = val;
        setVariableStatus(VariableStatus::changed);
    }

protected:

    /// Value stored
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

/// \class ComplexComponentParameter
/// \brief Abstract base allowing access to the components of a ComplexParameter as a RealParameter
/// \author Daniel Greenwald
class ComplexComponentParameter : public RealParameter
{
public:
    /// Constructor
    /// \param par shared_ptr to ComplexParameter to access
    ComplexComponentParameter(std::shared_ptr<ComplexParameter> par)
        : RealParameter(), Parent_(par)
    {
        if (!Parent_)
            throw exceptions::Exception("Parent unset", "RealSubparameter::RealSubparameter");
    }

    /// \return value of parameter by accessing parent
    double value() const override
    { return component(Parent_->value()); }

    /// set value by accessing parent
    void setValue(double val) override
    {
        if (variableStatus() == VariableStatus::fixed)
            throw exceptions::ParameterIsFixed("", "ComplexComponentParameter::setValue");
        if (value() == val)
            return;
        Parent_->setValue(setComponent(Parent_->value(), val));
        setVariableStatus(VariableStatus::changed);
    }

    /// \return shared_ptr to parent
    std::shared_ptr<ComplexParameter> parent()
    { return Parent_; }

protected:

    /// \return component value from whole value
    virtual double component(const std::complex<double>& c) const = 0;

    /// \return complex value with component changed
    virtual std::complex<double> setComponent(const std::complex<double>& c, double v) const = 0;

private:

    std::shared_ptr<ComplexParameter> Parent_;

};

/// \class RealComponentParameter
/// \brief RealParameter accessing real component of ComplexParameter
/// \author Daniel Greenwald
class RealComponentParameter : public ComplexComponentParameter
{
public:

    /// Constructor
    /// \param par shared_ptr to ComplexParameter to access
    RealComponentParameter(std::shared_ptr<ComplexParameter> par)
        : ComplexComponentParameter(par) {}

protected:

    /// \return component value from whole value
    double component(const std::complex<double>& c) const override
    { return real(c); }

    /// \return complex value with component changed
    virtual std::complex<double> setComponent(const std::complex<double>& c, double v) const
    { return std::complex<double>(v, imag(c)); }

};

/// \class ImaginaryComponentParameter
/// \brief ImaginaryParameter accessing imaginary component of ComplexParameter
/// \author Daniel Greenwald
class ImaginaryComponentParameter : public ComplexComponentParameter
{
public:

    /// Constructor
    /// \param par shared_ptr to ComplexParameter to access
    ImaginaryComponentParameter(std::shared_ptr<ComplexParameter> par)
        : ComplexComponentParameter(par) {}

protected:

    /// \return component value from whole value
    double component(const std::complex<double>& c) const override
    { return imag(c); }

    /// \return complex value with component changed
    virtual std::complex<double> setComponent(const std::complex<double>& c, double v) const
    { return std::complex<double>(real(c), v); }

};

}

#endif

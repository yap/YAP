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

#include "fwd/Parameter.h"

#include "Exceptions.h"
#include "VariableStatus.h"

#include <algorithm>
#include <complex>
#include <memory>
#include <numeric>

namespace yap {

/// \class ParameterBase
/// \brief Class holding basic properties of a parameter, but not a value!
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup Parameters
class ParameterBase
{
protected:

    /// Constructor
    constexpr explicit ParameterBase() noexcept : VariableStatus_(VariableStatus::changed) {}

public:

    /// \return VariableStatus
    VariableStatus& variableStatus() noexcept
    { return VariableStatus_; }

    /// \return VariableStatus (const)
    constexpr const VariableStatus variableStatus() const noexcept
    { return VariableStatus_; }

    /// \return number of real elements in parameter
    virtual const size_t size() const = 0;

    /// Set value from vector
    virtual const VariableStatus setValue(const std::vector<double>& V) = 0;

private:

    /// Status of variable
    VariableStatus VariableStatus_;

};

/// \return VariableStatus::changed if any element in range is changed; VariableStatus::unchanged otherwise
template <typename IterType>
constexpr VariableStatus variable_status(IterType first, IterType last)
{
    return std::any_of(first, last,
                       [](const std::shared_ptr<ParameterBase>& p)
                       {return p->variableStatus() == VariableStatus::changed;})
        ? VariableStatus::changed : VariableStatus::unchanged;
}

/// \class Parameter
/// \brief Template class holding also a value for a parameter
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Parameters
/// \tparam T type to be stored in Parameter
template <typename T>
class Parameter : public ParameterBase
{
public:

    /// The underlying type of the parameter.
    using parameter_type = T;

    /// Default constructor
    Parameter() = default;

    /// Value-assigning constructor
    constexpr explicit Parameter(T t) : ParameterBase(), ParameterValue_(t) {}

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
    virtual typename std::conditional<std::is_fundamental<T>::value, const T, const T&>::type
    value() const
    { return ParameterValue_; }

    using ParameterBase::setValue;

    /// set value
    virtual const VariableStatus setValue(typename std::conditional<std::is_fundamental<T>::value, const T, const T&>::type val)
    {
        if (variableStatus() == VariableStatus::fixed)
            throw exceptions::ParameterIsFixed("", "Parameter::setValue");
        if (ParameterValue_ == val)
            return variableStatus();
        ParameterValue_ = val;
        return variableStatus() = VariableStatus::changed;
    }

    /// set value by operator
    Parameter& operator=(typename std::conditional<std::is_fundamental<T>::value, const T, const T&>::type t)
    { setValue(t); return *this; }


private:

    /// Value stored
    T ParameterValue_;

};

/// \return string of Parameter
template <typename T>
inline std::string to_string(const Parameter<T>& P)
{ using std::to_string; return to_string(P.value()) + " (" + to_string(P.variableStatus()) + ")"; }

/// @brief Multiplication-equal operator for parameters.
/// @param par The parameter to update.
/// @param a   The multiplication factor.
/// @tparam P      The parameter type.
/// @tparam enable Control template parameter to abort the compilation if _P_ doesn't inherit from Parameter<>. 
template <typename P,
          typename enable = typename std::enable_if<std::is_base_of<Parameter<typename P::parameter_type>, P>::value>::type>
inline P& operator*=(P& par, typename P::parameter_type a)
{
    par.setValue(a * par.value());
    return par;
}

/// @brief Division-equal operator for parameters.
/// @param par The parameter to update.
/// @param a   The division factor.
/// @tparam P      The parameter type.
/// @tparam enable Control template parameter to abort the compilation if _P_ doesn't inherit from Parameter<>. 
template <typename P,
          typename enable = typename std::enable_if<std::is_base_of<Parameter<typename P::parameter_type>, P>::value>::type>
inline P& operator/=(P& par, typename P::parameter_type a)
{ return (par *= (1. / a)); }

/// @brief Addition-equal operator for parameters.
/// @param par The parameter to update.
/// @param a   The addition factor.
/// @tparam P      The parameter type.
/// @tparam enable Control template parameter to abort the compilation if _P_ doesn't inherit from Parameter<>. 
template <typename P,
          typename enable = typename std::enable_if<std::is_base_of<Parameter<typename P::parameter_type>, P>::value>::type>
inline P& operator+=(P& par, typename P::parameter_type a)
{
    par.setValue(a + par.value());
    return par;
}

/// @brief Subtraction-equal operator for parameters.
/// @param par The parameter to update.
/// @param a   The subtraction factor.
/// @tparam P      The parameter type.
/// @tparam enable Control template parameter to abort the compilation if _P_ doesn't inherit from Parameter<>. 
template <typename P,
          typename enable = typename std::enable_if<std::is_base_of<Parameter<typename P::parameter_type>, P>::value>::type>
inline P& operator-=(P& par, typename P::parameter_type a)
{ return (par += -a); }

/// \class RealParameter
/// \ingroup Parameters
class RealParameter : public Parameter<double>
{
public:

    /// constructor
    explicit RealParameter(double t = 0) : Parameter(t) {}

    /// \return size = 1
    const size_t size() const override
    { return 1; }

    using Parameter::setValue;

    /// Set value from vector
    const VariableStatus setValue(const std::vector<double>& V) override
    { return setValue(V[0]); }

    using Parameter::operator=;
};

/// \class NonnegativeRealParameter
/// \ingroup Parameters
class NonnegativeRealParameter : public RealParameter
{
public:
    /// Constructor
    explicit NonnegativeRealParameter(double t = 0) : RealParameter(t)
    {
        if (value() < 0)
            throw exceptions::Exception("value is negative", "NonnegativeRealParameter::NonnegativeRealParameter");
    }

    using RealParameter::setValue;

    /// Checks if value is non-negative
    const VariableStatus setValue(double t) override
    {
        if (t < 0)
            throw exceptions::Exception("value is negative", "NonnegativeRealParameter::setValue");
        return RealParameter::setValue(t);
    }

    using RealParameter::operator=;
};

/// \class PositiveRealParameter
/// \ingroup Parameters
class PositiveRealParameter : public NonnegativeRealParameter
{
public:
    /// Constructor
    explicit PositiveRealParameter(double t = 0) : NonnegativeRealParameter(t)
    {
        if (t <= 0)
            throw exceptions::Exception("value is not positive", "PositiveRealParameter::PositiveRealParameter");
    }
    
    using NonnegativeRealParameter::setValue;

    /// Checks if value is positive
    const VariableStatus setValue(double t) override
    {
        if (t <= 0)
            throw exceptions::Exception("value is not positive", "NonnegativeRealParameter::setValue");
        return RealParameter::setValue(t);
    }

    using NonnegativeRealParameter::operator=;
};

/// \class ComplexParameter
/// \ingroup Parameters
class ComplexParameter : public Parameter<std::complex<double> >
{
public:

    /// constructor
    explicit ComplexParameter(const std::complex<double>& t = 0) : Parameter(t) {}

    /// \return size = 2
    const size_t size() const override
    { return 2; }

    using Parameter::setValue;

    /// Set value from vector
    const VariableStatus setValue(const std::vector<double>& V) override
    { return setValue(std::complex<double>(V[0], V[1])); }

    using Parameter::operator=;
};


/// \class ComplexComponentParameter
/// \brief Abstract base allowing access to the components of a ComplexParameter as a RealParameter
/// \author Daniel Greenwald
class ComplexComponentParameter : public RealParameter
{
public:
    /// Constructor
    /// \param par shared_ptr to ComplexParameter to access
    explicit ComplexComponentParameter(std::shared_ptr<ComplexParameter> par) : RealParameter(), Parent_(par)
    { if (!Parent_) throw exceptions::Exception("Parent unset", "ComplexComponentParameter::ComplexComponentParameter"); }
        
    /// \return value of parameter by accessing parent
    const double value() const override
    { return component(Parent_->value()); }

    using RealParameter::setValue;

    /// set value by accessing parent
    const VariableStatus setValue(double val) override
    { return Parent_->setValue(setComponent(Parent_->value(), val)); }
        
    /// \return shared_ptr to parent
    std::shared_ptr<ComplexParameter> parent()
    { return Parent_; }

    using RealParameter::operator=;

protected:

    /// \return component value from whole value
    virtual double component(const std::complex<double>& c) const = 0;

    /// \return complex value with component changed
    virtual std::complex<double> setComponent(const std::complex<double>& c, double v) const = 0;

private:

    /// ComplexParameter this object points to a component of
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
    explicit RealComponentParameter(std::shared_ptr<ComplexParameter> par)
        : ComplexComponentParameter(par) {}

    using ComplexComponentParameter::operator=;

protected:

    /// \return component value from whole value
    double component(const std::complex<double>& c) const override
    { return real(c); }

    /// \return complex value with component changed
    virtual std::complex<double> setComponent(const std::complex<double>& c, double v) const override
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
    explicit ImaginaryComponentParameter(std::shared_ptr<ComplexParameter> par)
        : ComplexComponentParameter(par) {}

    using ComplexComponentParameter::operator=;

protected:

    /// \return component value from whole value
    double component(const std::complex<double>& c) const override
    { return imag(c); }

    /// \return complex value with component changed
    virtual std::complex<double> setComponent(const std::complex<double>& c, double v) const override
    { return std::complex<double>(real(c), v); }

};

/// \return number of real elements in ParameterVector
inline const size_t size(const ParameterVector& V)
{ return std::accumulate(V.begin(), V.end(), size_t(0), [](size_t& s, const ParameterVector::value_type & p) {return s += p->size();}); }

/// set value into parameter from iterator, calling parameter's setValue(vector<double>) function
/// \param P Parameter to set into
/// \param first Iterator to take values from
template <class InputIt>
std::vector<double>::const_iterator set_value(ParameterBase& P, InputIt first)
{ std::vector<double> V(first, first + P.size()); P.setValue(V); return first += (P.size() - 1); }

/// set values in ParameterVector from iterators
template <class InputIt>
void set_values(ParameterVector::iterator first_par, ParameterVector::iterator last_par, InputIt first_val, InputIt last_val)
{
    while (first_par != last_par and first_val != last_val) {
        first_val = set_value(**first_par, first_val);
        ++first_par;
        ++first_val;
    }
    if (first_par != last_par)
        throw exceptions::Exception("insufficient number of values provided", "set_values");
}

/// set values in ParameterVector from values in vector<double>
inline void set_values(ParameterVector& pars, const std::vector<double>& vals)
{ set_values(pars.begin(), pars.end(), vals.begin(), vals.end()); }

}

#endif

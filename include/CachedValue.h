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

#include "fwd/CachedValue.h"
#include "fwd/CalculationStatus.h"
#include "fwd/FourVector.h"
#include "fwd/Parameter.h"
#include "fwd/ParticleCombination.h"
#include "fwd/StatusManager.h"
#include "fwd/VariableStatus.h"

#include "DataAccessor.h"
#include "DataPoint.h"

#include <complex>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace yap {

/// \class CachedValue
/// \brief Class for managing cached values inside a #DataPoint
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
/// \ingroup Cache
class CachedValue : public std::enable_shared_from_this<CachedValue>
{

protected:

    /// Constructor (protected)
    /// \param size number of real elements in cached value
    /// \param da DataAccessor it to belong to
    CachedValue(unsigned size, DataAccessor& da);

public:

    /// \struct stores calculation and variable statuses for a CachedValue
    struct Status {

        /// constructor
        Status();

        /// Calculation status
        CalculationStatus Calculation;

        /// Variable status
        VariableStatus Variable;

        /// assignment of Calculation
        Status& operator=(const CalculationStatus& s)
        { Calculation = s; return *this; }

        /// assignment of Variable;
        /// does not change Variable if it is fixed!
        Status& operator=(const VariableStatus& s);
    };

    /// \name Getters
    /// @{

    /// \return raw pointer owning DataAccessor
    DataAccessor* owner() const
    { return Owner_; }

    /// \return index within owner
    const int index() const
    { return Index_; }

    /// Get value from #DataPoint for particular symmetrization
    /// \param index index of value to get from within cached value (must be less than #Size_)
    /// \param d #DataPoint to get value from
    /// \param sym_index index of symmetrization to grab from
    /// \return Value of CachedValue inside the data point
    inline const double value(unsigned index, const DataPoint& d, unsigned sym_index) const
    { return d.Data_[Owner_->index()][sym_index][Position_ + index]; }

    /// \return Size of cached value (number of real elements)
    virtual const unsigned size() const
    { return Size_; }

    /// @}

    /// \name Setters
    /// @{

    /// Set value into #DataPoint for particular symmetrization
    /// (No update to VariableStatus or CalculationStatus is made!)
    /// \param index index of value to get from within cached value (must be less than #Size_)
    /// \param val Value to set to
    /// \param d #DataPoint to update
    /// \param sym_index index of symmetrization to apply to
    void setValue(unsigned index, double val, DataPoint& d, unsigned sym_index) const
    { d.Data_[Owner_->index()][sym_index][Position_ + index] = val; }

    /// @}

    /// grant friend status to DataAccessor to set itself owner
    friend class DataAccessor;

protected:

    /// add to the Owner_
    void addToDataAccessor();

    /// set index
    void setIndex(int i)
    { Index_ = i; }

    /// set position
    void setPosition(int p)
    { Position_ = p; }

private:

    /// Owning DataAccessor
    DataAccessor* Owner_;

    /// index within owner
    int Index_;

    /// Position of first element of cached value within data vector
    int Position_;

    /// Size of cached value (number of real elements)
    unsigned Size_;

};

/// equality operator for checking the CalculationStatus
inline bool operator==(const CachedValue::Status& S, const CalculationStatus& s)
{ return S.Calculation == s; }

/// inequality operator for checking the CalculationStatus
inline bool operator!=(const CachedValue::Status& S, const CalculationStatus& s)
{ return S.Calculation != s; }

/// equality operator for checking the VariableStatus
inline bool operator==(const CachedValue::Status& S, const VariableStatus& s)
{ return S.Variable == s; }

/// inequality operator for checking the VariableStatus
inline bool operator!=(const CachedValue::Status& S, const VariableStatus& s)
{ return S.Variable != s; }

/// streaming operator for CachedValue::Status
std::string to_string(const CachedValue::Status& S);

/// \class RealCachedValue
/// \brief Class for managing a single real cached value inside a #DataPoint
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
/// \ingroup Cache
class RealCachedValue : public CachedValue
{
private:

    /// Constructor
    RealCachedValue(DataAccessor& da) : CachedValue(1, da) {}

public:

    /// create shared_ptr to RealCachedValue
    /// \param owner #DataAccessor to which this cached value belongs
    static std::shared_ptr<RealCachedValue> create(DataAccessor& da);

    /// Set value into #DataPoint for particular symmetrization, and
    /// update VariableStatus for symm. and partition index
    /// \param val Value to set to
    /// \param d #DataPoint to update
    /// \param sym_index index of symmetrization to apply to
    /// \param sm StatusManager
    void setValue(double val, DataPoint& d, unsigned sym_index, StatusManager& sm) const;

    /// Set value into #DataPoint for particular symmetrization;
    /// does not update any statuses
    /// \param val Value to set to
    /// \param d #DataPoint to update
    /// \param sym_index index of symmetrization to apply to
    void setValue(double val, DataPoint& d, unsigned sym_index) const
    { CachedValue::setValue(0, val, d, sym_index); }

    /// Get value from #DataPoint for particular symmetrization
    /// \param d #DataPoint to get value from
    /// \param sym_index index of symmetrization to grab from
    /// \return Value of CachedValue inside the data point
    const double value(const DataPoint& d, unsigned sym_index) const
    { return CachedValue::value(0, d, sym_index); }

};

/// \class ComplexCachedValue
/// \brief Class for managing a complex cached value inside a #DataPoint
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
/// \ingroup Cache
class ComplexCachedValue : public CachedValue
{
private:

    /// Constructor (protected)
    /// see #create for details
    ComplexCachedValue(DataAccessor& da) : CachedValue(2, da) {}

public:

    /// create shared pointer to ComplexCachedValue
    /// \param owner #DataAccessor to which this cached value belongs
    static std::shared_ptr<ComplexCachedValue> create(DataAccessor& da);

    /// Set value into #DataPoint for particular symmetrization, and
    /// update VariableStatus for symm. and partition index
    /// \param val Value to set to
    /// \param d #DataPoint to update
    /// \param sym_index index of symmetrization to apply to
    /// \param sm StatusManager
    void setValue(const std::complex<double>& val, DataPoint& d, unsigned sym_index, StatusManager& sm) const
    { setValue(real(val), imag(val), d, sym_index, sm); }

    /// Set value into #DataPoint for particular symmetrization, and
    /// update VariableStatus for symm. and partition index
    /// \param val_re real part of value to set to
    /// \param val_im imaginary part of value to set to
    /// \param d #DataPoint to update
    /// \param sym_index index of symmetrization to apply to
    /// \param sm StatusManager
    void setValue(double val_re, double val_im, DataPoint& d, unsigned sym_index, StatusManager& sm) const;

    /// Set value into #DataPoint for particular symmetrization;
    /// does not update any statuses
    /// \param val Value to set to
    /// \param d #DataPoint to update
    /// \param sym_index index of symmetrization to apply to
    void setValue(const std::complex<double>& val, DataPoint& d, unsigned sym_index) const
    { setValue(real(val), imag(val), d, sym_index); }

    /// Set value into #DataPoint for particular symmetrization
    /// does not update any statuses
    /// \param val_re real part of value to set to
    /// \param val_im imaginary part of value to set to
    /// \param d #DataPoint to update
    /// \param sym_index index of symmetrization to apply to
    void setValue(double val_re, double val_im, DataPoint& d, unsigned sym_index) const
    { CachedValue::setValue(0, val_re, d, sym_index); CachedValue::setValue(1, val_im, d, sym_index); }

    /// Get value from #DataPoint for particular symmetrization
    /// \param d #DataPoint to get value from
    /// \param sym_index index of symmetrization to grab from
    /// \return Value of CachedValue inside the data point
    const std::complex<double> value(const DataPoint& d, unsigned  sym_index) const
    { return std::complex<double>(CachedValue::value(0, d, sym_index), CachedValue::value(1, d, sym_index)); }

};

/// \class FourVectorCachedValue
/// \brief Class for managing a four-vector cached value inside a #DataPoint
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
/// \ingroup Cache
class FourVectorCachedValue : public CachedValue
{
private:

    /// Constructor (protected)
    /// see #create for details
    FourVectorCachedValue(DataAccessor& da) : CachedValue(4, da) {}

public:

    /// create shared pointer to ComplexCachedValue
    /// \param owner #DataAccessor to which this cached value belongs
    static std::shared_ptr<FourVectorCachedValue> create(DataAccessor& da);

    /// Set value into #DataPoint for particular symmetrization, and
    /// update VariableStatus for symm. and partition index
    /// \param val Value to set to
    /// \param d #DataPoint to update
    /// \param sym_index index of symmetrization to apply to
    /// \param sm StatusManager
    void setValue(const FourVector<double>& val, DataPoint& d, unsigned sym_index, StatusManager& sm) const;

    /// Set value into #DataPoint for particular symmetrization;
    /// does not update any statuses.
    /// \param val Value to set to
    /// \param d #DataPoint to update
    /// \param sym_index index of symmetrization to apply to
    void setValue(const FourVector<double>& val, DataPoint& d, unsigned sym_index) const;

    /// Get value from #DataPoint for particular symmetrization
    /// \param d #DataPoint to get value from
    /// \param sym_index index of symmetrization to grab from
    /// \return Value of CachedValue inside the data point
    const FourVector<double> value(const DataPoint& d, unsigned  sym_index) const;

};

}

#endif

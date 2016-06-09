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

#ifndef yap_CachedDataValue_h
#define yap_CachedDataValue_h

#include "fwd/CachedDataValue.h"
#include "fwd/CalculationStatus.h"
#include "fwd/FourVector.h"
#include "fwd/Parameter.h"
#include "fwd/ParticleCombination.h"
#include "fwd/StatusManager.h"
#include "fwd/VariableStatus.h"

#include "DataAccessor.h"
#include "DataPoint.h"

#include <memory>
#include <set>
#include <utility>
#include <vector>

namespace yap {

/// \struct DaughterCachedDataValue
/// \brief Stores a shared_ptr to a CachedDataValue and the index of
/// the daughter within a Particle Combination to pass to it
/// \ingroup Data
/// \ingroup Cache
struct DaughterCachedDataValue {
    /// constructor
    DaughterCachedDataValue(std::shared_ptr<CachedDataValue> cdv, size_t index)
        : CDV(cdv), Daughter(index) {}

    /// CachedDataValue
    std::shared_ptr<CachedDataValue> CDV;

    /// index of daughter to pass to it
    size_t Daughter;
};

/// less-than comparison operator
inline bool operator<(const DaughterCachedDataValue& A, const DaughterCachedDataValue& B)
{ return (A.CDV < B.CDV) or (A.Daughter < B.Daughter); }

/// \class CachedDataValueBase
/// \brief Class for managing cached values inside a #DataPoint
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
/// \ingroup Cache
class CachedDataValue : public std::enable_shared_from_this<CachedDataValue>
{
public:

    /// \struct stores calculation and variable statuses for a CachedDataValue
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
    /// \return Value of CachedDataValue inside the data point
    inline const double value(unsigned index, const DataPoint& d, unsigned sym_index) const
    {
#ifdef ELPP_DISABLE_DEBUG_LOGS
        return d.Data_[Owner_->index()][sym_index][Position_ + index];
#else
        return d.Data_.at(Owner_->index()).at(sym_index).at(Position_ + index);
#endif
    }

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
    void setValue(unsigned index, double val, DataPoint& d, unsigned sym_index) const;

    /// @}

    /// grant friend status to DataAccessor to set itself owner
    friend class DataAccessor;

protected:

    /// Constructor (protected)
    CachedDataValue(unsigned size, ParameterSet pars = {}, CachedDataValueSet vals = {});

    /// set the owning DataAccessor
    void setDataAccessor(DataAccessor* owner);

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
inline bool operator==(const CachedDataValue::Status& S, const CalculationStatus& s)
{ return S.Calculation == s; }

/// inequality operator for checking the CalculationStatus
inline bool operator!=(const CachedDataValue::Status& S, const CalculationStatus& s)
{ return S.Calculation != s; }

/// equality operator for checking the VariableStatus
inline bool operator==(const CachedDataValue::Status& S, const VariableStatus& s)
{ return S.Variable == s; }

/// inequality operator for checking the VariableStatus
inline bool operator!=(const CachedDataValue::Status& S, const VariableStatus& s)
{ return S.Variable != s; }

/// streaming operator for CachedDataValue::Status
std::ostream& operator<<(std::ostream& str, const CachedDataValue::Status& S);

/// \class RealCachedDataValue
/// \brief Class for managing a single real cached value inside a #DataPoint
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
/// \ingroup Cache
class RealCachedDataValue : public CachedDataValue
{
public:

    /// create shared_ptr to RealCachedDataValue
    /// \param owner #DataAccessor to which this cached value belongs
    /// \param pars set of shared pointers to Parameters cached value depends on
    /// \param vals set of shared pointers to CachedValues cached value depends on
    static std::shared_ptr<RealCachedDataValue> create(DataAccessor* da, ParameterSet pars = {}, CachedDataValueSet vals = {});

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
    { CachedDataValue::setValue(0, val, d, sym_index); }

    /// Get value from #DataPoint for particular symmetrization
    /// \param d #DataPoint to get value from
    /// \param sym_index index of symmetrization to grab from
    /// \return Value of CachedDataValue inside the data point
    const double value(const DataPoint& d, unsigned sym_index) const
    { return CachedDataValue::value(0, d, sym_index); }

private:

    /// Constructor
    RealCachedDataValue(ParameterSet pars = {}, CachedDataValueSet vals = {})
        : CachedDataValue(1, pars, vals) {}

};

/// \class ComplexCachedDataValue
/// \brief Class for managing a complex cached value inside a #DataPoint
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
/// \ingroup Cache
class ComplexCachedDataValue : public CachedDataValue
{
public:

    /// create shared pointer to ComplexCachedDataValue
    /// \param owner #DataAccessor to which this cached value belongs
    /// \param pars set of shared pointers to Parameters cached value depends on
    /// \param vals set of shared pointers to CachedValues cached value depends on
    static std::shared_ptr<ComplexCachedDataValue> create(DataAccessor* da, ParameterSet pars = {}, CachedDataValueSet vals = {});

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
    { CachedDataValue::setValue(0, val_re, d, sym_index); CachedDataValue::setValue(1, val_im, d, sym_index); }

    /// Get value from #DataPoint for particular symmetrization
    /// \param d #DataPoint to get value from
    /// \param sym_index index of symmetrization to grab from
    /// \return Value of CachedDataValue inside the data point
    const std::complex<double> value(const DataPoint& d, unsigned  sym_index) const
    { return std::complex<double>(CachedDataValue::value(0, d, sym_index), CachedDataValue::value(1, d, sym_index)); }

private:

    /// Constructor (protected)
    /// see #create for details
    ComplexCachedDataValue(ParameterSet pars = {}, CachedDataValueSet vals = {})
        : CachedDataValue(2, pars, vals) {}
};

/// \class FourVectorCachedDataValue
/// \brief Class for managing a four-vector cached value inside a #DataPoint
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
/// \ingroup Cache
class FourVectorCachedDataValue : public CachedDataValue
{
public:

    /// create shared pointer to ComplexCachedDataValue
    /// \param owner #DataAccessor to which this cached value belongs
    /// \param pars set of shared pointers to Parameters cached value depends on
    /// \param vals set of shared pointers to CachedValues cached value depends on
    static std::shared_ptr<FourVectorCachedDataValue> create(DataAccessor* da, ParameterSet pars = {}, CachedDataValueSet vals = {});

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
    /// \return Value of CachedDataValue inside the data point
    const FourVector<double> value(const DataPoint& d, unsigned  sym_index) const;

private:

    /// Constructor (protected)
    /// see #create for details
    FourVectorCachedDataValue(ParameterSet pars = {}, CachedDataValueSet vals = {})
        : CachedDataValue(4, pars, vals) {}


};

}

#endif

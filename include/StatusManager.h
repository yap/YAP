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

#ifndef yap_StatusManager_h
#define yap_StatusManager_h

#include "fwd/CalculationStatus.h"
#include "fwd/VariableStatus.h"

#include "CachedValue.h"
#include "DataAccessor.h"

namespace yap {

class StatusManager
{
public:

    /// constructor
    /// \param sDA DataAccessorSet to construct StatusManager for
    StatusManager(const DataAccessorSet& sDA);

    /// \name direct access to individual statuses
    /// @{

    /// Access status
    /// \return CachedValue::Status
    /// \param da_index Index of DataAccessor
    /// \param cdv_index Index of CachedValue
    /// \param sym_index Index of symmetrization
    CachedValue::Status& status(size_t da_index, size_t cdv_index, size_t sym_index)
    { return Statuses_[da_index][cdv_index][sym_index]; }

    /// retrieve status (const)
    /// \return CachedValue::Status (const)
    /// \param da_index Index of DataAccessor
    /// \param cdv_index Index of CachedValue
    /// \param sym_index Index of symmetrization
    const CachedValue::Status& status(size_t da_index, size_t cdv_index, size_t sym_index) const
    { return const_cast<StatusManager*>(this)->status(da_index, cdv_index, sym_index); }

    /// Access status
    /// \return CachedValue::Status
    /// \param cdv CachedValue
    /// \param sym_index Index of symmetrization
    CachedValue::Status& status(const CachedValue& cdv, size_t sym_index);

    /// retrieve status (const)
    /// \return CachedValue::Status (const)
    /// \param cdv CachedValue
    /// \param sym_index Index of symmetrization
    const CachedValue::Status& status(const CachedValue& cdv, size_t sym_index) const
    { return const_cast<StatusManager*>(this)->status(cdv, sym_index); }

    /// @}

    /// set all statuses for a particular CachedValue
    /// \param cdv CachedValue
    /// \param stat status to set to
    /// \tparam T type of status
    template <class T>
    void set(const CachedValue& cdv, const T& stat)
    {
        for (auto& s : Statuses_[cdv.owner()->index()][cdv.index()])
            s = stat;
    }

    /// set all statuses for all CachedValue's of a DataAccessor
    /// \param da DataAccessor
    /// \param stat status to set to
    /// \tparam T type of status
    template <class T>
    void set(const DataAccessor& da, const T& stat)
    {
        for (auto& v : Statuses_[da.index()])
            for (auto& s : v)
                s = stat;
    }

    /// set all statuses to a value
    /// \param stat value to set all statuses to
    /// \tparam T type of status
    template <class T>
    void setAll(const T& stat)
    {
        for (auto& v1 : Statuses_)
            for (auto& v2 : v1)
                for (auto& s : v2)
                    s = stat;
    }

    /// copy all calculation statuses from another manager
    /// \param sm StatusManager to copy from
    void copyCalculationStatuses(const StatusManager& sm);

private:

    /// vector of Status;
    /// first index is for DataAccessor;
    /// second index is for CachedValue
    /// third index is for SymmetrizationIndex
    std::vector<std::vector<std::vector<CachedValue::Status> > > Statuses_;

};

}
#endif

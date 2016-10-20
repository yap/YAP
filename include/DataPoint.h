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

#ifndef yap_DataPoint_h
#define yap_DataPoint_h

#include "fwd/DataAccessor.h"

#include <string>
#include <vector>

namespace yap {

/// \class DataPoint
/// \brief Class for holding data and cached values per data point for fast calculation
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup Data Data-related classes
class DataPoint
{
public:

    /// Constructor
    /// \param dataSet DataSet this DataPoint belongs to
    DataPoint(const DataAccessorSet& dataAccessorSet);

    /// \return number of data accessor rows
    size_t nDataAccessors() const
    { return Data_.size(); }

    /// \return size of data point
    unsigned bytes() const;

    /// check that two DataPoint's have same internal structure
    friend bool equalStructure(const DataPoint& A, const DataPoint& B);

    /// check that two DataPoint's are equal
    friend bool operator==(const DataPoint& lhs, const DataPoint& rhs)
    { return lhs.Data_ == rhs.Data_; }

    /// grant friend status to CachedValue to access Data_
    friend class CachedValue;

private:

    /// Data storage for all DataAccessors
    /// first index is for the DataAccessor
    /// second index is for the symmeterization state (as known by the DataAccessor)
    /// third index is internal to the DataAccessor
    std::vector<std::vector<double> > Data_;

};


}

#endif

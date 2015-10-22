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

#ifndef yap_DataPartition_h
#define yap_DataPartition_h

#include "CalculationStatus.h"
#include "DataPoint.h"

#include <memory>
#include <vector>

namespace yap {

typedef std::vector<DataPoint>::iterator DataIterator;

/// \class DataPartition
/// \brief Class defining a partition of the DataSet
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data

class DataPartition
{
public:

    /// Constructor
    DataPartition(const DataPoint& dataPoint, DataIterator begin, DataIterator end, unsigned spacing = 1);

    /// increment and
    /// \return if still in range
    bool increment();

    /// \return current DataPoint
    DataPoint& dataPoint()
    { return *CurrentPosition_; }

    /// \return current DataPoint (const)
    const DataPoint& dataPoint() const
    { return *CurrentPosition_; }

    /// \return DataPartition's index
    unsigned index() const
    { return DataPartitionIndex_; }

private:
    unsigned DataPartitionIndex_;

    DataIterator CurrentPosition_;

    DataIterator Begin_;
    DataIterator End_;
    unsigned Spacing_;

};

}

#endif

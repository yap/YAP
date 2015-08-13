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

#ifndef yap_DataSet_h
#define yap_DataSet_h

#include <vector>

namespace yap {

class DataPoint;

/// \class DataSet
/// \brief Class holding a set of DataPoint objects.
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data

class DataSet
{
public:

    /// Add data point via move
    /// \param d DataPoint to move into DataSet
    /// \return Success of action
    bool addDataPoint(DataPoint&& d);

    /// Add data point via copy
    /// \param d DataPoint to copy into DataSet
    /// \return Success of action
    bool addDataPoint(const DataPoint& d);

    /// Check if data point is consisent with data set
    bool consisent(const DataPoint& d) const;

private:

    /// Vector of data points
    std::vector<DataPoint> DataPoints_;
};

}

#endif

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

#include "fwd/DataPoint.h"
#include "fwd/FourVector.h"
#include "fwd/Model.h"
#include "fwd/StatusManager.h"

#include "DataPartition.h"

#include <vector>

namespace yap {

/// \class DataSet
/// \brief Class holding a set of DataPoint objects.
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
class DataSet : public DataPartitionBlock
{
public:

    /// Constructor
    DataSet(const Model& m);

    /// Check if data point is consisent with data set
    bool consistent(const DataPoint& d) const;

    /// add empty data points
    /// \param n number of points to add
    void addEmptyDataPoints(size_t n);

    /// \typedef value_type
    /// needed by std::back_inserter and std::inserter
    using value_type = std::vector<FourVector<double> >;

    /// creates a new #DataPoint and calls #Model::setFinalStateMomenta on it
    /// \param P the momenta to be set
    /// \param sm StatusManager to use when calling #Model::setFinalStateMomenta
    /// \return the created #DataPoint
    const DataPoint createDataPoint(const std::vector<FourVector<double> >& P, StatusManager& sm);

    /// creates a new #DataPoint and #setFinalStateMomenta to P
    /// \param P the momenta to be set
    /// \return the created #DataPoint
    const DataPoint createDataPoint(const std::vector<FourVector<double> >& P)
    { return createDataPoint(P, *this); }

    /// creates a DataPoint from a vector of FourVector's using #createDataPoint
    /// and pushes it back in the vector data points
    /// \param P vector of FourVector's to create DataPoint from
    void push_back(const std::vector<FourVector<double> >& P)
    { DataPoints_.push_back(createDataPoint(P)); }

    /// checks consistency of DataPoint, and pushes it back in the vector of data points
    /// \param d DataPoint to copy into DataSet
    void push_back(const DataPoint& d);

    /// checks consistency of DataPoint, and pushes it back in the vector of data points
    /// \param d DataPoint to move into DataSet
    void push_back(DataPoint&& d);

    /// creates a DataPoint from a vector of FourVector's using #createDataPoint
    /// and inserts it into the vector of DataPoint's at a specified position
    /// \param pos DataIterator of position in DataSet to insert into
    /// \param P vector of FourVector to create DataPoint from
    DataIterator insert(DataIterator pos, const std::vector<FourVector<double> >& P)
    { return dataIterator(DataPoints_.insert(rawIterator(pos), createDataPoint(P)), pos.partition()); }

    /// checks consistency of DataPoint and inserts it into vector of
    /// data points at specifief position
    /// \param pos DataIterator of position in DataSet to insert into
    /// \param d DataPoint to copy into DataSet
    DataIterator insert(const DataIterator& pos, const DataPoint& d);

    /// checks consistency of DataPoint and inserts it into vector of
    /// data points at specifief position
    /// \param pos DataIterator of position in DataSet to insert into
    /// \param d DataPoint to move into DataSet
    DataIterator insert(const DataIterator& pos, DataPoint&& d);

    /// clear the data set
    void clear()
    { DataPoints_.clear(); }

    /// removes last element added to data set
    void pop_back()
    { DataPoints_.pop_back(); }

    /// remove specified element from data set
    /// \param pos iterator to element to remove
    DataIterator erase(const DataIterator& pos)
    { return dataIterator(DataPoints_.erase(rawIterator(pos)), pos.partition()); }

    /// remove specified elements from data set
    /// \param first iterator to first element to remove
    /// \param last iterator beyond last element to remove
    DataIterator erase(const DataIterator& first, const DataIterator& last);

    /// \return iterator to front of set
    const DataIterator& begin() const override
    { return const_cast<DataSet*>(this)->setBegin(const_cast<DataPointVector*>(&DataPoints_)->begin()); }

    /// \return iterator to end of set
    const DataIterator& end() const override
    { return const_cast<DataSet*>(this)->setEnd(const_cast<DataPointVector*>(&DataPoints_)->end()); }

    /// access by index
    DataPoint& operator[](size_t i)
    { return DataPoints_[i]; }

    /// access by index (with check)
    DataPoint& at(size_t i)
    { return DataPoints_.at(i); }

    /// access front
    const DataPoint& front() const
    { return DataPoints_.front(); }

    /// access back
    const DataPoint& back() const
    { return DataPoints_.back(); }

    /// \return number of data points
    const size_t size() const override
    { return DataPoints_.size(); }

    /// \return maximum possible number data points that can be stored
    const size_t max_size() const
    { return DataPoints_.max_size(); }

    /// \return number of data points that can be held in currently allocated storage
    const size_t capacity() const
    { return DataPoints_.capacity(); }

    /// \return whether DataPoints_ is empty
    const bool empty() const
    { return DataPoints_.empty(); }

    /// reserve storage space
    void reserve(size_t n)
    { DataPoints_.reserve(n); }

    /// call shrink to fit on DataPoints_
    void shrink_to_fit()
    { DataPoints_.shrink_to_fit(); }

    /// \return raw pointer to associated model
    const Model* model() const
    { return Model_; }

    /// grant friend status to DataPartition to access non-const dataPoints()
    friend DataPartition;

protected:

    /// non-const access to DataPoints_
    DataPointVector& dataPoints()
    { return DataPoints_; }

private:
    /// vector of data points contained in set
    DataPointVector DataPoints_;

    /// Associated model
    const Model* Model_;

};

}

#endif

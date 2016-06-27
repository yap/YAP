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

#include "fwd/FourVector.h"
#include "fwd/Model.h"

#include "DataPoint.h"
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

    /// Copy constructor
    DataSet(const DataSet& other);

    /// Move constructor
    DataSet(DataSet&& other);

    /// Copy assignment operator
    DataSet& operator=(const DataSet& other);

    /// Move assignment operator
    DataSet& operator=(DataSet&& other);

    /// Swap
    void swap(DataSet& other);

    /// Check if data point is consisent with data set
    bool consistent(const DataPoint& d) const;

    /// add empty data points
    /// \param n number of points to add
    void addEmptyDataPoints(size_t n);

	/// creates a new #DataPoint and #setFinalStateMomenta to P
	/// \param P the momenta to be set
	/// \return the created #DataPoint
	const DataPoint createDataPoint(const std::vector<FourVector<double> >& P);

	void push_back(const std::vector<FourVector<double> >& P);
	void push_back(std::vector<FourVector<double> >&& P);
	void insert(DataIterator pos, const std::vector<FourVector<double> >& P);
	void insert(DataIterator pos, std::vector<FourVector<double> >&& P);

    /// set four momenta of data point
    /// \param P vector of FourVectors of final-state momenta
    /// \param sm StatusManager to update
    void setFinalStateMomenta(DataPoint& d, const std::vector<FourVector<double> >& P, StatusManager& sm);

    /// set four momenta of data point
    /// \param P vector of FourVectors of final-state momenta
    void setFinalStateMomenta(DataPoint& d, const std::vector<FourVector<double> >& P);

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

    /// access back
    DataPoint& back()
    { return DataPoints_.back(); }

    size_t size() const
    { return DataPoints_.size(); }

    /// const access to DataPoints_
    const DataPointVector& points() const
    { return DataPoints_; }

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

/// swap
inline void swap(DataSet& A, DataSet& B)
{ A.swap(B); }

}

namespace std {

/// swap
template <>
inline void swap<yap::DataSet>(yap::DataSet& A, yap::DataSet& B)
{ A.swap(B); }

}

#endif

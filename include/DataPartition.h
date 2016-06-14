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

#include "fwd/DataPartition.h"

#include "fwd/DataPoint.h"
#include "fwd/DataSet.h"

#include "StatusManager.h"

#include <algorithm>
#include <iterator>
#include <memory>

namespace yap {

/// \class DataIterator
/// \brief Class for iterating over a #DataPartition
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
class DataIterator : public std::iterator<std::random_access_iterator_tag, DataPoint>
{
public:

    /// pre-increment operator
    DataIterator& operator++();

	/// post-increment operator
	/// \todo return object or reference?
	DataIterator operator++(int);

	/// addition assignment operator
	DataIterator& operator+=(DataIterator::difference_type n)
	{ this->Iterator_ += n; return *this;}

	/// addition operator
	friend DataIterator operator+(const DataIterator& lhs, DataIterator::difference_type n)
	{ return DataIterator(*lhs.Partition_, lhs.Iterator_ + n); }

	/// addition operator (make it commutative)
	friend DataIterator operator+(DataIterator::difference_type n, const DataIterator& rhs)
	{ return rhs + n; }

	/// pre-decrement operator
	DataIterator& operator--();

	/// post-decrement operator
	/// \todo return object or reference?
	DataIterator operator--(int);

	/// subtraction assignment operator
	DataIterator& operator-=(DataIterator::difference_type n)
	{ this->Iterator_ -= n; return *this;}

	/// subraction operator
	const DataIterator operator-(DataIterator::difference_type n) const
	{ return DataIterator(*Partition_, Iterator_ - n); }

	/// subraction operator
	friend const DataIterator::difference_type operator-(const DataIterator& lhs, const DataIterator& rhs)
	{ return lhs.Iterator_ - rhs.Iterator_; }

    /// dereference operator
    DataPoint& operator*()
    { return *Iterator_; }

    /// dereference operator (const)
    const DataPoint& operator*() const
    { return *Iterator_; }

	/// pointer operator
	DataPoint operator->()
	{ return *(this->Iterator_).operator->(); }

    /// check ownership
    bool ownedBy(const DataPartition& dp) const
    { return Partition_ == &dp; }

	/// less-than operator
	friend const bool operator<(const DataIterator& lhs, const DataIterator& rhs)
	{ return lhs.Iterator_ < rhs.Iterator_; }

	/// greater-than operator
	friend const bool operator>(const DataIterator& lhs, const DataIterator& rhs)
	{ return lhs.Iterator_ > rhs.Iterator_; }

	/// less-than-or-equal operator
	friend const bool operator<=(const DataIterator& lhs, const DataIterator& rhs)
	{ return !(lhs > rhs); }

	/// greater-than-or-equal operator
	friend const bool operator>=(const DataIterator& lhs, const DataIterator& rhs)
	{ return !(lhs < rhs); }

	/// equality operator
	friend const bool operator==(const DataIterator& lhs, const DataIterator& rhs)
	{ return lhs.Iterator_ == rhs.Iterator_; }

	/// inequality operator
	friend const bool operator!=(const DataIterator& lhs, const DataIterator& rhs)
	{ return !(lhs == rhs); }

	/// access operator
	DataPoint operator[](DataIterator::difference_type n) const
	{ return this->Iterator_[n]; }

    /// grant friend status to DataPartition to access Iterator_
    friend DataPartition;

protected:

    /// constructor with defaulted iterator
    /// \param p owning DataPartition
    DataIterator(const DataPartition& p)
        : Partition_(&p) {}

    /// constructor
    /// \param p owning DataPartition
    /// \param it vector<DataPoint> iterator to contain
    DataIterator(const DataPartition& p, DataPointVector::iterator it)
        : Partition_(&p), Iterator_(it) {}

    /// owning DataPartition
    const DataPartition* Partition_;

    /// iterator within vector<DataPoint>
    DataPointVector::iterator Iterator_;

};


/// \class DataPartition
/// \brief Class defining a partition of the DataSet
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
class DataPartition : public StatusManager
{
protected:

    /// Constructor
    /// \param sm StatusManager to copy StatusManager structure from
    /// \param begin vector<DataPoint>::iterator of start
    /// \param end vector<DataPoint>::iterator of end
    DataPartition(const StatusManager& sm, DataPointVector::iterator begin, DataPointVector::iterator end)
        : StatusManager(sm), Begin_(*this, begin), End_(*this, end) {}

    /// constructor taking a DataAccessorSet
    /// \param sDA DataAccessorSet to initialize StatusManager from
    DataPartition(const DataAccessorSet& sDA)
        : StatusManager(sDA), Begin_(*this), End_(*this) {}

public:

    /// virtual destructor (defaulted)
    virtual ~DataPartition() = default;

    /// copy constructor (defaulted)
    DataPartition(const DataPartition&) = default;

    /// move constructor (defaulted)
    DataPartition(DataPartition&&) = default;

    /// copy assignment operator (defaulted)
    DataPartition& operator=(const DataPartition&) = default;

    /// move assignment operator (defaulted)
    DataPartition& operator=(DataPartition&&) = default;

    /// \return begin iterator
    virtual const DataIterator& begin() const
    { return Begin_; }

    /// \return end iterator
    virtual const DataIterator& end() const
    { return End_; }

    /// grant friend status to DataIterator to call increment
    friend DataIterator;

    /// grant friend status to Model to call ...
    friend class Model;

protected:

    /// increment iterator;
    /// \attention Must be overloaded in derived classes
    virtual void increment(DataIterator& it, DataIterator::difference_type n = 1) const
    { it = End_; }

    /// decrement iterator;
    /// \attention Must be overloaded in derived classes
    virtual void decrement(DataIterator& it, DataIterator::difference_type n = 1) const
    { it = End_; }

    /// \return vector<DataPoint> iterator inside DataIterator
    DataPointVector::iterator& rawIterator(DataIterator& it) const
    { return it.Iterator_; }

    /// set begin
    const DataIterator& setBegin(DataPointVector::iterator it)
    { Begin_ = DataIterator(*this, it); return Begin_; }

    /// set end
    const DataIterator& setEnd(DataPointVector::iterator it)
    { End_ = DataIterator(*this, it); return End_; }

    /// get non-const begin from DataSet
    static DataPointVector::iterator begin(DataSet& ds);

    /// get non-const end from DataSet
    static DataPointVector::iterator end(DataSet& ds);

private:

    /// begin DataIterator
    DataIterator Begin_;

    /// end DataIterator
    DataIterator End_;

};

/// \class DataPartitionBlock
/// \brief A contiguous block of data
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
class DataPartitionBlock : public DataPartition
{
public:

    /// Constructor
    /// \param sm StatusManager to copy StatusManager structure from
    /// \param begin vector<DataPoint>::iterator of start
    /// \param end vector<DataPoint>::iterator of end
    DataPartitionBlock(const StatusManager& sm, DataPointVector::iterator begin, DataPointVector::iterator end)
        : DataPartition(sm, begin, end) {}

    /// \return DataParitionVector covering DataSet as contiguous blocks
    /// \param dataSet The dataSet
    /// \param n number of partitions to divide the dataSet into
    static DataPartitionVector create(DataSet& dataSet, unsigned n);

    /// \return DataParitionVector covering DataSet as contiguous blocks of specified size
    /// \param dataSet The dataSet
    /// \param s maximum size of partitions to divide the dataSet into
    static DataPartitionVector createBySize(DataSet& dataSet, size_t s);

protected:

    /// constructor taking a DataAccessorSet
    /// \param sDA DataAccessorSet to initialize StatusManager from
    DataPartitionBlock(const DataAccessorSet& sDA)
        : DataPartition(sDA) {}

    /// increment DataIterator
    /// \param it DataIterator to iterate
    virtual void increment(DataIterator& it, DataIterator::difference_type n = 1) const override;

    /// decrement DataIterator
    /// \param it DataIterator to iterate
    virtual void decrement(DataIterator& it, DataIterator::difference_type n = 1) const override;
};

/// \class DataPartitionWeave
/// \brief A set of data spaced over the range [B,E) with spacing S = [B+0S, B+1S, B+2S, B+3S, ..., E)
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
class DataPartitionWeave : public DataPartition
{
public:

    /// Constructor
    /// \param sm StatusManager to copy StatusManager structure from
    /// \param begin vector<DataPoint>::iterator of start
    /// \param end vector<DataPoint>::iterator of end
    /// \param spacing Spacing between consecutively evaluated points
    DataPartitionWeave(const StatusManager& sm, DataPointVector::iterator begin, DataPointVector::iterator end, unsigned spacing)
        : DataPartition(sm, begin, end), Spacing_(spacing) {}

    /// \return DataParitionVector covering DataSet as a weave
    /// \param dataSet The dataSet
    /// \param n number of partitions to divide the dataSet into
    static DataPartitionVector create(DataSet& dataSet, unsigned n);

protected:

    /// constructor taking a DataAccessorSet
    /// \param sDA DataAccessorSet to initialize StatusManager from
    /// \param spacing for weave
    DataPartitionWeave(const DataAccessorSet& sDA, unsigned spacing)
        : DataPartition(sDA), Spacing_(spacing) {}

    /// increment DataIterator
    /// \param it DataIterator to iterate
    virtual void increment(DataIterator& it, DataIterator::difference_type n = 1) const override;

	/// decrement DataIterator
	/// \param it DataIterator to iterate
	virtual void decrement(DataIterator& it, DataIterator::difference_type n = 1) const override;

    /// spacing between data points for the weaving
    unsigned Spacing_;
};

}

#endif

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

#include "DataPoint.h"
#include "DataSet.h"
#include "logging.h"

#include <algorithm>
#include <iterator>
#include <vector>

namespace yap {

class DataPartitionBase;

/// \class DataIterator
/// \brief Class for iterating over a #DataPartitionBase
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data

class DataIterator : virtual public std::forward_iterator_tag
{
public:
    DataIterator& operator++();

    DataPoint& operator*()
    { return *Iterator_; }

    const DataPoint& operator*() const
    { return *Iterator_; }

    bool operator!=(const DataIterator& it) const
    { return Iterator_ != it.Iterator_; }

    bool ownedBy(DataPartitionBase* dpb) const
    { return Owner_ == dpb; }

    friend DataPartitionBase;

protected:

    // Protected constructor
    DataIterator(DataPartitionBase* p, std::vector<DataPoint>::iterator& it)
        : Owner_(p), Iterator_(it)
    {}

    DataPartitionBase* Owner_;
    std::vector<DataPoint>::iterator Iterator_;
};


/// \class DataPartitionBase
/// \brief Class defining a partition of the DataSet
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data

class DataPartitionBase
{
public:
    /// Constructor
    DataPartitionBase(std::vector<DataPoint>::iterator begin, std::vector<DataPoint>::iterator end);

    const DataIterator& begin() const
    { return Begin_; }

    const DataIterator& end() const
    { return End_; }

    /// \return DataPartition's index
    unsigned index() const
    { return DataPartitionIndex_; }

    friend DataIterator;

protected:

    /// increment and
    /// \return if still in range
    virtual void increment(DataIterator& it) = 0;

    std::vector<DataPoint>::iterator& rawIterator(DataIterator& it)
    { return it.Iterator_; }

    DataIterator Begin_;
    DataIterator End_;

    unsigned DataPartitionIndex_;

};

/// \class DataPartitionWeave
/// \brief A set of data spaced over the range [B,E) with spacing S = [B+0S, B+1S, B+2S, B+3S, ..., E)
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data

class DataPartitionWeave : public DataPartitionBase
{
public:
    /// Constructor
    DataPartitionWeave(std::vector<DataPoint>::iterator begin, std::vector<DataPoint>::iterator end,
                       unsigned spacing)
        : DataPartitionBase(begin, end),
          Spacing_(spacing)
    {}

protected:
    virtual void increment(DataIterator& it) override;

    unsigned Spacing_;
};

/// \class DataPartitionBlock
/// \brief A contiguous block of data
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data

class DataPartitionBlock : public DataPartitionWeave
{
public:
    /// Constructor
    DataPartitionBlock(std::vector<DataPoint>::iterator begin, std::vector<DataPoint>::iterator end)
        : DataPartitionWeave(begin, end, 1)
    {}
};

// DataPartitionCreators

/// create DataPartitionWeave'es
/// \param dataSet The dataSet
/// \param nPartitions number of partitions to divide the dataSet into
inline std::vector<DataPartitionBase*> createDataPartitionsWeave(DataSet& dataSet, unsigned nPartitions)
{
    DEBUG("Partition dataSet of size " << dataSet.size() << " into " << nPartitions << " interweaved partitions");

    std::vector<DataPartitionBase*> partitions;
    partitions.reserve(nPartitions);

    for (unsigned i = 0; i < nPartitions; ++i) {
        DEBUG("Create DataPartitionWeave with size " << unsigned(0.5 + double(std::distance(dataSet.begin() + i, dataSet.end())) / nPartitions));
        partitions.push_back(new DataPartitionWeave(dataSet.begin() + i, dataSet.end(), nPartitions));
    }

    return partitions;
}

/// create DataPartitionBlock's
/// \param dataSet The dataSet
/// \param maxBlockSize maximum size of DataPoints a block can have
inline std::vector<DataPartitionBase*> createDataPartitionsBlockBySize(DataSet& dataSet, unsigned maxBlockSize)
{
    DEBUG("Partition dataSet of size " << dataSet.size() << " into blocks with a maximum size of " << maxBlockSize);

    std::vector<DataPartitionBase*> partitions;
    unsigned dataSetSize = dataSet.size();

    maxBlockSize = std::min(dataSetSize, maxBlockSize);
    maxBlockSize = std::max(dataSetSize / unsigned(double(dataSetSize) / maxBlockSize + 0.5), maxBlockSize);

    partitions.reserve(dataSetSize / maxBlockSize + 1);

    auto begin = dataSet.begin();
    auto end   = dataSet.begin() + maxBlockSize;

    while (true) {
        DEBUG("Create DataPartitionBlock with size " << std::distance(begin, end));
        partitions.push_back(new DataPartitionBlock(begin, end));

        if (end >= dataSet.end())
            break;

        begin = end;
        end += maxBlockSize;
        if (end >= dataSet.end())
            end = dataSet.end();
    }

    return partitions;
}

/// create DataPartitionBlock's
/// \param dataSet The dataSet
/// \param nPartitions number of blocks to divide the dataSet into
inline std::vector<DataPartitionBase*> createDataPartitionsBlock(DataSet& dataSet, unsigned nPartitions)
{
    DEBUG("Partition dataSet of size " << dataSet.size() << " into " << nPartitions << " partition blocks");

    unsigned partitionSize(dataSet.size());

    if (nPartitions < dataSet.size())
        partitionSize = 0.5 + double(dataSet.size()) / nPartitions;

    partitionSize = std::max(1u, partitionSize);

    return createDataPartitionsBlockBySize(dataSet, partitionSize);
}


}

#endif

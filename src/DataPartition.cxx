#include "DataPartition.h"

#include "logging.h"
#include "make_unique.h"

namespace yap {

//-------------------------
DataIterator& DataIterator::operator++()
{
    Owner_->increment(*this);
    return *this;
}

//-------------------------
DataPartitionBase::DataPartitionBase(std::vector<DataPoint>::iterator begin, std::vector<DataPoint>::iterator end) :
    Begin_(this, begin),
    End_(this, end),
    DataPartitionIndex_(0)
{

}

//-------------------------
void DataPartitionWeave::increment(DataIterator& it)
{
    // check that iterator belongs to this partition
    if (!it.ownedBy(this)) {
        LOG(FATAL) << "DataPartition::increment - called with DataIterator belonging to different partition!";
        it = End_;
        return;
    }

    for (unsigned i = 0; i < Spacing_ && it != End_; ++i)
        ++rawIterator(it);
}

//-------------------------
std::vector<std::unique_ptr<DataPartitionBase> > createDataPartitionsWeave(DataSet& dataSet, unsigned nPartitions)
{
    LOG(INFO) << "Partition dataSet of size " << dataSet.size() << " into " << nPartitions << " interweaved partitions";

    std::vector<std::unique_ptr<DataPartitionBase> > partitions;
    partitions.reserve(nPartitions);

    for (unsigned i = 0; i < nPartitions; ++i) {
        LOG(INFO) << "Create DataPartitionWeave with size " << unsigned(0.5 + double(std::distance(dataSet.begin() + i, dataSet.end())) / nPartitions);
        partitions.push_back(std::make_unique<DataPartitionWeave>(dataSet.begin() + i, dataSet.end(), nPartitions));
    }

    return partitions;
}

//-------------------------
std::vector<std::unique_ptr<DataPartitionBase> > createDataPartitionsBlockBySize(DataSet& dataSet, unsigned maxBlockSize)
{
    LOG(INFO) << "Partition dataSet of size " << dataSet.size() << " into blocks with a maximum size of " << maxBlockSize;

    std::vector<std::unique_ptr<DataPartitionBase> > partitions;
    unsigned dataSetSize = dataSet.size();

    maxBlockSize = std::min(dataSetSize, maxBlockSize);
    maxBlockSize = std::max(dataSetSize / unsigned(double(dataSetSize) / maxBlockSize + 0.5), maxBlockSize);

    partitions.reserve(dataSetSize / maxBlockSize + 1);

    auto begin = dataSet.begin();
    auto end   = dataSet.begin() + maxBlockSize;

    while (true) {
        LOG(INFO) << "Create DataPartitionBlock with size " << std::distance(begin, end);
        partitions.push_back(std::make_unique<DataPartitionBlock>(begin, end));

        if (end >= dataSet.end())
            break;

        begin = end;
        end += maxBlockSize;
        if (end >= dataSet.end())
            end = dataSet.end();
    }

    return partitions;
}

//-------------------------
std::vector<std::unique_ptr<DataPartitionBase> > createDataPartitionsBlock(DataSet& dataSet, unsigned nPartitions)
{
    LOG(INFO) << "Partition dataSet of size " << dataSet.size() << " into " << nPartitions << " partition blocks";

    unsigned partitionSize(dataSet.size());

    if (nPartitions < dataSet.size())
        partitionSize = 0.5 + double(dataSet.size()) / nPartitions;

    partitionSize = std::max(1u, partitionSize);

    return createDataPartitionsBlockBySize(dataSet, partitionSize);
}

}

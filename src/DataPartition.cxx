#include "DataPartition.h"

#include "DataPoint.h"
#include "DataSet.h"
#include "Exceptions.h"
#include "logging.h"
#include "make_unique.h"

namespace yap {

//-------------------------
DataIterator& DataIterator::operator++()
{
    Partition_->increment(*this);
    return *this;
}

//-------------------------
DataIterator DataIterator::operator++(int)
{
	DataIterator it(*Partition_, Iterator_);
    Partition_->increment(*this);
    return it;
}

//-------------------------
DataIterator& DataIterator::operator--()
{
    Partition_->decrement(*this);
    return *this;
}

//-------------------------
DataIterator DataIterator::operator--(int)
{
	DataIterator it(*Partition_, Iterator_);
    Partition_->decrement(*this);
    return it;
}

//-------------------------
DataPointVector::iterator DataPartition::begin(DataSet& ds)
{
    return ds.dataPoints().begin();
}

//-------------------------
DataPointVector::iterator DataPartition::end(DataSet& ds)
{
    return ds.dataPoints().end();
}

//-------------------------
void DataPartitionWeave::increment(DataIterator& it, DataIterator::difference_type n) const
{
	auto distanceToEnd = end() - it;
	rawIterator(it) += ((distanceToEnd < n * Spacing_) ? (n * Spacing_) : distanceToEnd);
}

//-------------------------
void DataPartitionWeave::decrement(DataIterator& it, DataIterator::difference_type n) const
{
	auto distanceFromBegin = it - begin();
	rawIterator(it) -= ((distanceFromBegin > n * Spacing_) ? (n * Spacing_) : distanceFromBegin);
}

//-------------------------
DataPartitionVector DataPartitionWeave::create(DataSet& dataSet, unsigned n)
{
    if (n == 0)
        throw exceptions::Exception("number of partitions is zero", "DataParitionWeave::create");

    auto N = dataSet.points().size();

    LOG(INFO) << "Partitioning data set of size " << N << " into " << n << " interwoven partitions";

    DataPartitionVector P;
    P.reserve(n);

    for (unsigned i = 0; i < n; ++i) {
        LOG(INFO) << "Creating DataPartitionWeave with size " << std::ceil(1.*(N - i) / n);
        P.push_back(std::make_unique<DataPartitionWeave>(dataSet, begin(dataSet) + i, end(dataSet), n));
    }

    return P;
}

//-------------------------
void DataPartitionBlock::increment(DataIterator& it, DataIterator::difference_type n) const
{
    rawIterator(it) += n; 
}

//-------------------------
void DataPartitionBlock::decrement(DataIterator& it, DataIterator::difference_type n) const
{
    rawIterator(it) -= n;
}

//-------------------------
DataPartitionVector DataPartitionBlock::create(DataSet& dataSet, unsigned n)
{
    if (n == 0)
        throw exceptions::Exception("number of partitions is zero", "DataParitionBlock::create");

    auto N = dataSet.points().size();

    if (n > N)
        n = N;

    LOG(INFO) << "Partitioning data set of size " << N << " into " << n << " contiguous blocks";

    unsigned p_size = std::round(N / n);

    DataPartitionVector P;
    P.reserve(n);

    auto it_b = begin(dataSet);

    for (unsigned i = 0; i < n - 1; ++i) {
        auto it_e = it_b + p_size;
        LOG(INFO) << "Creating DataPartitionBlock with size " << std::distance(it_b, it_e);
        P.push_back(std::make_unique<DataPartitionBlock>(dataSet, it_b, it_e));
        it_b = it_e;
    }
    P.push_back(std::make_unique<DataPartitionBlock>(dataSet, it_b, end(dataSet)));

    return P;
}

//-------------------------
DataPartitionVector DataPartitionBlock::createBySize(DataSet& dataSet, size_t s)
{
    if (s == 0)
        throw exceptions::Exception("block size is zero", "DataPartitionBlock::createBySize");

    auto N = dataSet.points().size();

    if (s > N)
        s = N;

    LOG(INFO) << "Partitioning data set of size " << N << " into blocks with a maximum size of " << s;

    DataPartitionVector P;
    P.reserve(std::ceil(N / s));

    auto it_b = begin(dataSet);
    while (it_b != end(dataSet)) {
        auto it_e = it_b + s;
        LOG(INFO) << "Creating DataPartitionBlock with size " << std::distance(it_b, it_e);
        P.push_back(std::make_unique<DataPartitionBlock>(dataSet, it_b, it_e));
        it_b = it_e;
    }

    return P;
}

}

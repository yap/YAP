#include "DataPartition.h"

#include "DataPoint.h"
#include "DataSet.h"
#include "Exceptions.h"
#include "logging.h"
#include "make_unique.h"

namespace yap {

//-------------------------
DataIterator& DataIterator::operator+=(DataIterator::difference_type n)
{
    return Partition_->increment(*this, n);
}

//-------------------------
const DataIterator::difference_type operator-(const DataIterator& lhs, const DataIterator& rhs)
{
    if (lhs.Partition_ != rhs.Partition_)
        throw exceptions::Exception("DataIterator's belong to different DataPartition's", "operator-");
    return lhs.Partition_->difference(lhs.Iterator_, rhs.Iterator_);
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
DataIterator& DataPartitionWeave::increment(DataIterator& it, DataIterator::difference_type n) const
{
    rawIterator(it) += n * Spacing_;
    return it;
}

//-------------------------
const DataIterator::difference_type DataPartition::difference(const DataPointVector::iterator& lhs, const DataPointVector::iterator& rhs) const
{
    return static_cast<DataIterator::difference_type>(lhs - rhs);
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
        P.push_back(std::unique_ptr<DataPartitionBlock>(new DataPartitionBlock(dataSet, it_b, it_e)));
        it_b = it_e;
    }
    P.push_back(std::unique_ptr<DataPartitionBlock>(new DataPartitionBlock(dataSet, it_b, end(dataSet))));

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
        P.push_back(std::unique_ptr<DataPartitionBlock>(new DataPartitionBlock(dataSet, it_b, it_e)));
        it_b = it_e;
    }

    return P;
}

}

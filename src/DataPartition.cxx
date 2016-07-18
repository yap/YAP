#include "DataPartition.h"

#include "DataPoint.h"
#include "DataSet.h"
#include "Exceptions.h"

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

    DataPartitionVector P;
    P.reserve(n);

    for (unsigned i = 0; i < n; ++i)
        P.push_back(new DataPartitionWeave(dataSet, begin(dataSet) + i, end(dataSet), n));

    return P;
}

//-------------------------
DataPartitionVector DataPartitionBlock::create(DataSet& dataSet, unsigned n)
{
    if (n == 0)
        throw exceptions::Exception("number of partitions is zero", "DataParitionBlock::create");

    auto N = dataSet.size();
    n = std::min<unsigned>(n, N);

    unsigned p_size = std::round(N / n);

    DataPartitionVector P;
    P.reserve(n);

    auto it_b = begin(dataSet);

    for (unsigned i = 0; i < n - 1; ++i) {
        auto it_e = it_b + p_size;
        P.push_back(new DataPartitionBlock(dataSet, it_b, it_e));
        it_b = it_e;
    }
    P.push_back(new DataPartitionBlock(dataSet, it_b, end(dataSet)));

    return P;
}

//-------------------------
DataPartitionVector DataPartitionBlock::createBySize(DataSet& dataSet, size_t s)
{
    if (s == 0)
        throw exceptions::Exception("block size is zero", "DataPartitionBlock::createBySize");

    auto N = dataSet.size();
    s = std::min(s, N);

    DataPartitionVector P;
    P.reserve(std::ceil(N / s));

    auto it_b = begin(dataSet);
    while (it_b != end(dataSet)) {
        auto it_e = it_b + s;
        P.push_back(new DataPartitionBlock(dataSet, it_b, it_e));
        it_b = it_e;
    }

    return P;
}

}

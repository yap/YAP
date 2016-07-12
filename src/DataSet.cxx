#include "DataSet.h"

#include "DataPoint.h"
#include "Exceptions.h"
#include "Model.h"

namespace yap {

//-------------------------
DataSet::DataSet(const Model& m) :
    DataPartitionBlock(m.dataAccessors()),
    Model_(&m)
{
}

//-------------------------
bool DataSet::consistent(const DataPoint& d) const
{
    return DataPoints_.empty() or equalStructure(DataPoints_.front(), d);
}

//-------------------------
void DataSet::addEmptyDataPoints(size_t n)
{
    if (n == 0)
        return;

    if (!model())
        throw exceptions::Exception("Model unset or deleted", "DataSet::addEmptyDataPoints");

    // create first data point, either via DataAccessorSet constructor or copy constructor
    DataPoints_.emplace_back(DataPoints_.empty() ? model()->dataAccessors() : DataPoints_.back());

    // create remaining data points via copy constructor
    for (size_t i = 1; i < n; ++i)
        DataPoints_.emplace_back(DataPoints_.back());
}

//-------------------------
const DataPoint DataSet::createDataPoint(const std::vector<FourVector<double> >& P, StatusManager& sm)
{
    DataPoint d = DataPoints_.empty() ? DataPoint(model()->dataAccessors()) : DataPoints_.back();
    model()->setFinalStateMomenta(d, P, sm);
    return d;
}

//-------------------------
void DataSet::push_back(const DataPoint& d)
{
    if (!consistent(d))
        throw exceptions::InconsistentDataPoint("DataSet::push_back");
    DataPoints_.push_back(d);
}

//-------------------------
void DataSet::push_back(DataPoint&& d)
{
    if (!consistent(d))
        throw exceptions::InconsistentDataPoint("DataSet::push_back");
    DataPoints_.push_back(std::move(d));
}

//-------------------------
DataIterator DataSet::insert(const DataIterator& pos, const DataPoint& d)
{
    if (!consistent(d))
        throw exceptions::InconsistentDataPoint("DataSet::push_back");
    return dataIterator(DataPoints_.insert(rawIterator(pos), d), pos.partition());
}

//-------------------------
DataIterator DataSet::insert(const DataIterator& pos, DataPoint&& d)
{
    if (!consistent(d))
        throw exceptions::InconsistentDataPoint("DataSet::push_back");
    return dataIterator(DataPoints_.insert(rawIterator(pos), std::move(d)), pos.partition());
}

//-------------------------
DataIterator DataSet::erase(const DataIterator& first, const DataIterator& last)
{
    if (first.partition() != last.partition())
        throw exceptions::Exception("Iterators' partitions don't match", "DataSet::erase");
    if (first.partition() == this)
        return dataIterator(DataPoints_.erase(rawIterator(first), rawIterator(last)), first.partition());
    auto it = first;
    while (it != last)
        erase(it++);
    return last;
}

}

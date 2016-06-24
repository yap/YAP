#include "DataSet.h"

#include "DataPoint.h"
#include "Exceptions.h"
#include "FourMomenta.h"
#include "Model.h"

namespace yap {

//-------------------------
DataSet::DataSet(const Model& m) :
    DataPartitionBlock(m.dataAccessors()),
    Model_(&m)
{
}

//-------------------------
DataSet::DataSet(const DataSet& other) :
    DataPartitionBlock(other),
    DataPoints_(other.DataPoints_),
    Model_(other.Model_)
{
}

//-------------------------
DataSet::DataSet(DataSet&& other) :
    DataPartitionBlock(std::move(other)),
    DataPoints_(std::move(other.DataPoints_)),
    Model_(std::move(other.Model_))
{
}

//-------------------------
DataSet& DataSet::operator=(const DataSet& other)
{
    DataPartitionBlock::operator=(other);
    Model_ = other.Model_;
    DataPoints_ = other.DataPoints_;
    return *this;
}

//-------------------------
DataSet& DataSet::operator=(DataSet&& other)
{
    DataPartitionBlock::operator=(std::move(other));
    Model_ = std::move(other.Model_);
    DataPoints_ = std::move(other.DataPoints_);
    return *this;
}

//-------------------------
void DataSet::swap(DataSet& other)
{
    using std::swap;
    swap(static_cast<DataPartitionBlock&>(*this), static_cast<DataPartitionBlock&>(other));
    swap(Model_, other.Model_);
    swap(DataPoints_, other.DataPoints_);
}

//-------------------------
bool DataSet::consistent(const DataPoint& d) const
{
    return points().empty() or equalStructure(points().front(), d);
}

//-------------------------
void DataSet::setFinalStateMomenta(DataPoint& d, const std::vector<FourVector<double> >& P, StatusManager& sm)
{
    if (!model())
        throw exceptions::Exception("Model unset", "DataPoint::setFinalStateMomenta");

    model()->fourMomenta()->setFinalStateMomenta(d, P, sm);

    // call calculate on all static data accessors in model
    for (auto& sda : model()->staticDataAccessors())
        sda->calculate(d, sm);
}

//-------------------------
void DataSet::setFinalStateMomenta(DataPoint& d, const std::vector<FourVector<double> >& P)
{
    setFinalStateMomenta(d, P, *this);
}

//-------------------------
void DataSet::addEmptyDataPoints(size_t n)
{
    if (!model())
        throw exceptions::Exception("Model unset or deleted", "DataSet::add");

	for (size_t i = 0; i < n; ++i) {
		DataPoints_.emplace_back(model()->dataAccessors());

		// check if the created DataPoint is consistent
		if (!consistent(DataPoints_.back()))
			throw exceptions::Exception("produced inconsistent data point", "Model::addDataPoint");
	}
}

//-------------------------
const DataPoint DataSet::createDataPoint(const std::vector<FourVector<double> >& P)
{
	DataPoint d = ((points().empty())
			? DataPoint(model()->dataAccessors())
			: DataPoints_.back());

	setFinalStateMomenta(d, P);
	return d;
}

//-------------------------
void DataSet::push_back(const std::vector<FourVector<double> >& P)
{
	DataPoints_.push_back(createDataPoint(P));
}

}

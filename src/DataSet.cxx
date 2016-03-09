#include "DataSet.h"

#include "DataPoint.h"
#include "Exceptions.h"
#include "Model.h"

namespace yap {

//-------------------------
DataSet::DataSet(const Model& m) :
    DataPartitionBlock(m.dataAccessors()),
    ReportsModel(),
    Model_(&m)
{
}

//-------------------------
bool DataSet::consistent(const DataPoint& d) const
{
    return points().empty() or equalStructure(points().front(), d);
}

//-------------------------
void DataSet::addEmptyPoint()
{
    if (!model())
        throw exceptions::Exception("Model unset or deleted", "DataSet::add");

    DataPoints_.emplace_back(this);
    auto& d = DataPoints_.back();

    if (!consistent(d))
        throw exceptions::Exception("produced inconsistent data point", "Model::addDataPoint");
}

//-------------------------
void DataSet::addEmptyPoints(size_t n)
{
    for (size_t i = 0; i < n; ++i)
        addEmptyPoint();
}

//-------------------------
void DataSet::add(const std::vector<FourVector<double> >& P)
{
    addEmptyPoint();
    DataPoints_.back().setFinalStateMomenta(P);
}

//-------------------------
bool operator==(const DataSet& lhs, const DataSet& rhs)
{ return lhs.Model_ == rhs.Model_ and lhs.DataPoints_ == rhs.DataPoints_; }

}

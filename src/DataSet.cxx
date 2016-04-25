#include "DataSet.h"

#include "DataPoint.h"
#include "Exceptions.h"
#include "Model.h"

namespace yap {

//-------------------------
DataSet::DataSet(const Model& m) :
    DataPartitionBlock(m.dataAccessors()),
    ReportsModel(),
    Model_(&m),
    GlobalStatusManager_(m.dataAccessors())
{
}

//-------------------------
DataSet::DataSet(const DataSet& other) :
    DataPartitionBlock(other),
    ReportsModel(),
    DataPoints_(other.DataPoints_),
    Model_(other.Model_),
    GlobalStatusManager_(other.GlobalStatusManager_)
{
    assertDataPointOwnership();
}

//-------------------------
DataSet::DataSet(DataSet&& other) :
    DataPartitionBlock(std::move(other)),
    ReportsModel(),
    DataPoints_(std::move(other.DataPoints_)),
    Model_(std::move(other.Model_)),
    GlobalStatusManager_(std::move(other.GlobalStatusManager_))
{
    assertDataPointOwnership();
}

//-------------------------
DataSet& DataSet::operator=(const DataSet& other)
{
    DataPartitionBlock::operator=(other);
    Model_ = other.Model_;
    DataPoints_ = other.DataPoints_;
    GlobalStatusManager_ = other.GlobalStatusManager_;
    assertDataPointOwnership();
    return *this;
}

//-------------------------
DataSet& DataSet::operator=(DataSet&& other)
{
    DataPartitionBlock::operator=(std::move(other));
    Model_ = std::move(other.Model_);
    DataPoints_ = std::move(other.DataPoints_);
    GlobalStatusManager_ = std::move(other.GlobalStatusManager_);
    assertDataPointOwnership();
    return *this;
}

//-------------------------
void swap(DataSet& A, DataSet& B)
{
    std::swap(static_cast<DataPartitionBlock&>(A), static_cast<DataPartitionBlock&>(B));
    std::swap(A.Model_, B.Model_);
    std::swap(A.DataPoints_, B.DataPoints_);
    A.assertDataPointOwnership();
    B.assertDataPointOwnership();
}

//-------------------------
void DataSet::assertDataPointOwnership()
{
    for (auto& d : DataPoints_)
        d.DataSet_ = this;
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

    DataPoints_.emplace_back(*this);
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

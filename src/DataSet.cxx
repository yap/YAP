#include "DataSet.h"

#include "DataPoint.h"

#include "logging.h"

namespace yap {
    
//-------------------------
bool DataSet::addDataPoint(DataPoint&& d)
{
    if (!consisent(d))
        return false;
    DataPoints_.push_back(d);
    return true;
}

//-------------------------
bool DataSet::addDataPoint(const DataPoint& d)
{
    return addDataPoint(DataPoint(d));
}

//-------------------------
bool DataSet::consisent(const DataPoint& d) const
{
    bool result = true;

    if (!DataPoints_.empty() && d.fourMomenta().size() != DataPoints_[0].fourMomenta().size()) {
        LOG(ERROR) << "DataSet::consistent(DataPoint) - DataPoint has wrong number of four momenta (" << d.fourMomenta().size() << " != " << DataPoints_[0].fourMomenta().size() << ")";
        result = false;
    }

    return result;
}

}

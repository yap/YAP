#include "DataSet.h"

#include "DataPoint.h"
#include "FourMomenta.h"
#include "logging.h"

namespace yap {

//-------------------------
bool DataSet::addDataPoint(DataPoint&& d)
{
    // consistency check is done in InitialStateParticle
    DataPoints_.push_back(d);
    return true;
}

//-------------------------
bool DataSet::addDataPoint(const DataPoint& d)
{
    return addDataPoint(DataPoint(d));
}

//-------------------------
bool DataSet::consistent(const DataPoint& d) const
{
    bool result = true;

    if (!DataPoints_.empty() && d.FourMomenta_.size() != DataPoints_[0].FourMomenta_.size()) {
        LOG(ERROR) << "DataSet::consistent(DataPoint) - DataPoint has wrong number of four momenta (" << d.FourMomenta_.size() << " != " << DataPoints_[0].FourMomenta_.size() << ")";
        result = false;
    }

    if (d.CalculationStatuses_.size() != d.Data_.size()) {
        LOG(ERROR) << "DataSet::consistent(DataPoint) - number of CalculationStatuses (" << d.CalculationStatuses_.size() << ") != Data size (" << DataPoints_[0].FourMomenta_.size() << ")";
        result = false;
    }

    return result;
}

}

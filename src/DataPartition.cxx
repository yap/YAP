#include "DataPartition.h"

namespace yap {

//-------------------------
DataPartition::DataPartition(const DataPoint& dataPoint, DataIterator begin, DataIterator end, unsigned spacing) :
    CurrentPosition_(begin),
    Begin_(begin),
    End_(end),
    Spacing_(spacing)
{
    // initialize CalculationStatusesDataSet_ with kUncalculated flags
    CalculationStatusesDataSet_.resize(dataPoint.CachedAmplitudes_.size());
    unsigned i(0);
    for (auto& v : dataPoint.CachedAmplitudes_) {
        CalculationStatusesDataSet_[i++].assign(v.size(), kUncalculated);
    }

    CalculationStatusesDataPoint_ = CalculationStatusesDataSet_;
}

}

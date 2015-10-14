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

//-------------------------
bool DataPartition::increment()
{
    // advance
    CurrentPosition_ += Spacing_;

    // reset CalculationStatusesDataPoint_ to CalculationStatusesDataSet_
    CalculationStatusesDataPoint_ = CalculationStatusesDataSet_;

    bool retVal = CurrentPosition_  + Spacing_ < End_;

    if (!retVal) {
        // reset
        CurrentPosition_ = Begin_;

        // set flags to calculated
        for (auto& v : CalculationStatusesDataSet_)
            for (auto& s : v)
                s = kCalculated;

        CalculationStatusesDataPoint_ = CalculationStatusesDataSet_;
    }


    return retVal;
}

}

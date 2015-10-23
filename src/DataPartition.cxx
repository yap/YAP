#include "DataPartition.h"

namespace yap {

//-------------------------
DataPartition::DataPartition(const DataPoint& dataPoint, DataIterator begin, DataIterator end, unsigned spacing) :
    CurrentPosition_(begin),
    Begin_(begin),
    End_(end),
    Spacing_(spacing)
{
    // set DataPartitionIndex_
    /// \todo make this smarter?
    static unsigned index(0);
    DataPartitionIndex_ = (index++);
}

//-------------------------
bool DataPartition::increment()
{
    // advance
    CurrentPosition_ += Spacing_;

    bool retVal = CurrentPosition_  + Spacing_ <= End_;

    if (!retVal) {
        // reset
        CurrentPosition_ = Begin_;
    }

    return retVal;
}

}

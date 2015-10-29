#include "DataPartition.h"

#include "logging.h"

namespace yap {

//-------------------------
DataIterator& DataIterator::operator++()
{
    Owner_->increment(*this);
    return *this;
}

//-------------------------
DataPartitionBase::DataPartitionBase(std::vector<DataPoint>::iterator begin, std::vector<DataPoint>::iterator end) :
    Begin_(this, begin),
    End_(this, end)
{
    // set DataPartitionIndex_
    /// \todo make this smarter?
    static unsigned index(0);
    DataPartitionIndex_ = (index++);
}

//-------------------------
void DataPartitionWeave::increment(DataIterator& it)
{
    // check that iterator belongs to this partition
    if (!it.ownedBy(this)) {
        LOG(FATAL) << "DataPartition::increment - called with DataIterator belonging to different partition!";
        it = End_;
        return;
    }

    for (unsigned i = 0; i < Spacing_ && it != End_; ++i)
        ++rawIterator(it);
}


}

#include "DataAccessor.h"

namespace yap {

//-------------------------
DataAccessor::DataAccessor() :
    Recalculate_(true),
    Index_(0)
{
    // assign a running index to this DataAccessor
    // TODO: come up with something smarter
    static unsigned int index(0);
    Index_ = index++;
}
/*
//-------------------------
DataAccessor::DataAccessor(const DataAccessor& other) :
    Recalculate_(other.Recalculate_),
    Index_(0)
{
    // increments index
    static unsigned int index(0);
    Index_ = index++;
}

//-------------------------
DataAccessor::DataAccessor(DataAccessor&& other) :
    Recalculate_(std::move(other.Recalculate_)),
    Index_(std::move(other.Index_))
{
    // steals other's index
}

//-------------------------
DataAccessor& DataAccessor::operator=(DataAccessor&& rhs)
{
    Recalculate_ = std::move(rhs.Recalculate_);
    Index_ = std::move(rhs.Index_);
    return *this;
}
*/
}

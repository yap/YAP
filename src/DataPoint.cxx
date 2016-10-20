#include "DataPoint.h"

#include "DataAccessor.h"

namespace yap {

//-------------------------
DataPoint::DataPoint(const DataAccessorSet& dataAccessorSet)
    : Data_(dataAccessorSet.size())
{
    for (auto da : dataAccessorSet)
        Data_[da->index()].resize(da->nSymmetrizationIndices() * da->size(), 0.);
}

//-------------------------
bool equalStructure(const DataPoint& A, const DataPoint& B)
{
    if (A.Data_.size() != B.Data_.size())
        return false;
    for (size_t i = 0; i < A.Data_.size(); ++i) {
        if (A.Data_[i].size() != B.Data_[i].size())
            return false;
    }
    return true;
}

//-------------------------
unsigned DataPoint::bytes() const
{
    unsigned size = sizeof(Data_);
    for (auto& v : Data_) {
        size += sizeof(v);
        for (auto& vv : v)
            size += sizeof(vv);
    }
    return size;
}

}

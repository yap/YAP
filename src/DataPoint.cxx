#include "DataPoint.h"

#include "DataAccessor.h"

namespace yap {

//-------------------------
DataPoint::DataPoint(const DataAccessorSet& dataAccessorSet)
    : Data_(dataAccessorSet.size())
{
    for (auto da : dataAccessorSet)
        Data_[da->index()].assign(da->nSymmetrizationIndices(), std::vector<double>(da->size()));
}

//-------------------------
bool equalStructure(const DataPoint& A, const DataPoint& B)
{
    if (A.Data_.size() != B.Data_.size())
        return false;
    for (size_t i = 0; i < A.Data_.size(); ++i) {
        if (A.Data_[i].size() != B.Data_[i].size())
            return false;
        for (size_t j = 0; j < A.Data_[i].size(); ++j) {
            if (A.Data_[i][j].size() != B.Data_[i][j].size())
                return false;
        }
    }
    return true;
}

//-------------------------
unsigned DataPoint::bytes() const
{
    unsigned size = sizeof(Data_);
    for (auto& v : Data_) {
        size += sizeof(v);
        for (auto& vv : v) {
            size += sizeof(vv);
            for (auto vvv : vv)
                size += sizeof(vvv);
        }
    }
    return size;
}

}

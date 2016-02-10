#include "DataPoint.h"

#include "FourMomenta.h"
#include "HelicityAngles.h"
#include "logging.h"
#include "MeasuredBreakupMomenta.h"

#include <assert.h>
#include <iomanip>
#include <iostream>

namespace yap {

//-------------------------
DataPoint::DataPoint(const DataAccessorSet& S) :
    Data_(S.size())
{
    for (auto da : S)
        Data_[da->index()].assign(da->maxSymmetrizationIndex() + 1, std::vector<double>(da->size(), 0));
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
unsigned DataPoint::dataSize() const
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

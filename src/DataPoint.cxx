#include "DataPoint.h"

#include "DataSet.h"
#include "Exceptions.h"
#include "FourMomenta.h"
#include "HelicityAngles.h"
#include "logging.h"
#include "Model.h"
#include "StaticDataAccessor.h"

namespace yap {

//-------------------------
DataPoint::DataPoint(DataSet& dataSet) :
    DataSet_(&dataSet)
{
    if (!model())
        throw exceptions::Exception("Model unset", "DataPoint::DataPoint");

    Data_.resize(model()->dataAccessors().size());
    for (auto da : model()->dataAccessors())
        Data_[da->index()].assign(da->nSymmetrizationIndices(), std::vector<double>(da->size(), 0));
}

//-------------------------
const Model* DataPoint::model() const
{
    return DataSet_->model();
}

//-------------------------
void DataPoint::setFinalStateMomenta(const std::vector<FourVector<double> >& P, StatusManager& sm)
{
    if (!model())
        throw exceptions::Exception("Model unset", "DataPoint::setFinalStateMomenta");

    model()->fourMomenta()->setFinalStateMomenta(*this, P, sm);

    // call calculate on all static data accessors in model
    for (auto& sda : model()->staticDataAccessors())
        sda->calculate(*this, sm);

    FDEBUG("set final state momenta and calculated StaticDataAccessors")
}

//-------------------------
void DataPoint::setFinalStateMomenta(const std::vector<FourVector<double> >& P)
{
    setFinalStateMomenta(P, *DataSet_);
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

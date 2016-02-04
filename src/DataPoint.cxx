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
bool DataPoint::setFinalStateFourMomenta(const std::vector<FourVector<double> >& fourMomenta, bool check)
{
    if (fourMomenta.size() != FSPFourMomenta_.size())
        throw exceptions::Exception("fourMomenta is wrong size (" + std::to_string(fourMomenta.size())
                                    + "!= " + std::to_string(FSPFourMomenta_.size()) + ")",
                                    "DataPoint::setFinalStateFourMomenta");

    if (check and FSPFourMomenta_ == fourMomenta)
        return true;

    FSPFourMomenta_ = fourMomenta;
    return false;
}

//-------------------------
void DataPoint::allocateStorage(std::shared_ptr<FourMomenta> fourMom, const DataAccessorSet& dataAccessors)
{
    FourMomenta_.resize(fourMom->maxSymmetrizationIndex() + 1);

    // allocate space in vectors
    Data_.resize(dataAccessors.size());

    for (auto d : dataAccessors) {
        Data_[d->index()].assign(d->maxSymmetrizationIndex() + 1, std::vector<double>(d->size(), 0));
        FLOG(INFO) << "assigned  " << d->data_accessor_type() << " at index" << d->index() << " a vector of size "
                   << Data_[d->index()][0].size();
    }
}

//-------------------------
void DataPoint::printDataSize()
{
    unsigned totSize(0);

    unsigned size = sizeof(FSPFourMomenta_);
    size += FSPFourMomenta_.size() * sizeof(FourVector<double>);
    std::cout << "  Size of FSPFourMomenta_:   " << std::right << std::setw(5) << size << " byte  \tNumber of Indices: " << FSPFourMomenta_.size() << "\n";
    totSize += size;

    size = sizeof(FourMomenta_);
    size += FourMomenta_.size() * sizeof(FourVector<double>);
    std::cout << "  Size of FourMomenta_:      " << std::right << std::setw(5) << size << " byte  \tNumber of Indices: " << FourMomenta_.size() << "\n";
    totSize += size;

    size = sizeof(Data_);
    for (std::vector<std::vector<double> >& v : Data_) {
        size += sizeof(v);
        for (std::vector<double>& vv : v) {
            size += sizeof(vv);
            size += vv.size() * sizeof(double);
        }
    }
    std::cout << "+ Size of Data_:             " << std::right << std::setw(5) << size << " byte  \tNumber of Indices: " << Data_.size() << "\n";
    totSize += size;

    std::cout << "= Size of DataPoint:         " << std::right << std::setw(5) << totSize << " byte\n";
}

}

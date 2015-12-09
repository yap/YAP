#include "DataPoint.h"

#include "Constants.h"
#include "FourMomenta.h"
#include "HelicityAngles.h"
#include "logging.h"
#include "MeasuredBreakupMomenta.h"

#include <assert.h>
#include <iomanip>
#include <iostream>

namespace yap {

//-------------------------
DataPoint::DataPoint(const std::vector<FourVector<double> >& P) :
    FSPFourMomenta_(P)
{}

//-------------------------
bool DataPoint::setFinalStateFourMomenta(const std::vector<FourVector<double> >& fourMomenta)
{
    if (FSPFourMomenta_.size() != fourMomenta.size()) {
        LOG(ERROR) << "DataPoint::setFourMomenta - fourMomenta have wrong size "
                   << fourMomenta.size() << " != " << FSPFourMomenta_.size();
        return false;
    }

    FSPFourMomenta_ = fourMomenta;
    return true;
}

//-------------------------
void DataPoint::allocateStorage(const FourMomenta& fourMom, const std::set<DataAccessor*> dataAccessors)
{
    FourMomenta_.resize(fourMom.maxSymmetrizationIndex() + 1);

    // allocate space in vectors
    Data_.resize(dataAccessors.size());

    for (DataAccessor* d : dataAccessors) {
        Data_[d->index()].assign(d->maxSymmetrizationIndex() + 1, std::vector<double>(d->size()));
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

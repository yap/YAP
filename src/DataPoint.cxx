#include "DataPoint.h"

#include "Constants.h"
#include "FourMomenta.h"
#include "HelicityAngles.h"
#include "logging.h"

#include <assert.h>
#include <iostream>

namespace yap {

//-------------------------
DataPoint::DataPoint(const std::vector<TLorentzVector>& P) :
    FourMomenta_(P)
{}

//-------------------------
bool DataPoint::setFourMomenta(const std::vector<TLorentzVector>& fourMomenta)
{
    if (FourMomenta_.size() < fourMomenta.size()) {
        LOG(ERROR) << "DataPoint::setFourMomenta - fourMomenta have wrong size "
                   << fourMomenta.size() << " > " << FourMomenta_.size();
        return false;
    }

    std::copy(fourMomenta.begin(), fourMomenta.end(), FourMomenta_.begin());
    return true;
}

//-------------------------
void DataPoint::allocateStorage(const FourMomenta& fourMom, const HelicityAngles& helAngles, const std::set<DataAccessor*> dataAccessors)
{
    FourMomenta_.resize(fourMom.maxSymmetrizationIndex() + 1);
    HelicityAngles_.resize(helAngles.maxSymmetrizationIndex() + 1);

    // initialize helicity angles
    for (auto& helAngles : HelicityAngles_)
        helAngles = {0, 0};

    // allocate space in vectors
    Data_.resize(dataAccessors.size());
    CachedAmplitudes_.resize(dataAccessors.size());

    for (DataAccessor* d : dataAccessors) {
        Data_[d->index()].assign(d->maxSymmetrizationIndex() + 1, std::vector<double>(d->size()));
        CachedAmplitudes_[d->index()].assign(d->maxSymmetrizationIndex() + 1, Complex_0);
    }
}

//-------------------------
void DataPoint::printDataSize()
{
    unsigned totSize(0);

    unsigned size = sizeof(FourMomenta_);
    size += FourMomenta_.size() * sizeof(TLorentzVector);
    std::cout << "  Size of FourMomenta_:         " << size << " byte  \tNumber of Indices: " << FourMomenta_.size() << "\n";

    size = sizeof(HelicityAngles_);
    for (std::vector<double>& v : HelicityAngles_) {
        size += sizeof(v);
        size += v.size() * sizeof(double);
    }
    std::cout << "+ Size of HelicityAngles_:      " << size << " byte  \tNumber of Indices: " << HelicityAngles_.size() << "\n";
    totSize += size;

    size = sizeof(Data_);
    for (std::vector<std::vector<double> >& v : Data_) {
        size += sizeof(v);
        for (std::vector<double>& vv : v) {
            size += sizeof(vv);
            size += vv.size() * sizeof(double);
        }
    }
    std::cout << "+ Size of Data_:                " << size << " byte  \tNumber of Indices: " << Data_.size() << "\n";
    totSize += size;

    size = sizeof(CachedAmplitudes_);
    for (std::vector<std::complex<double> >& v : CachedAmplitudes_) {
        size += sizeof(v);
        size += v.size() * sizeof(std::complex<double>);
    }
    std::cout << "+ Size of CachedAmplitudes_:    " << size << " byte  \tNumber of Indices: " << CachedAmplitudes_.size() << "\n";
    totSize += size;

    std::cout << "= Size of DataPoint:            " << totSize << " byte\n";
}

}

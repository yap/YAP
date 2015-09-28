#include "DataPoint.h"

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

    CachedAmplitudes_.resize(dataAccessors.size());

    for (DataAccessor* d : dataAccessors)
        CachedAmplitudes_.at(d->index()) = std::vector<Amp>(d->maxSymmetrizationIndex() + 1, 0);
}

//-------------------------
void DataPoint::printDataSize()
{
    unsigned totSize(0);

    unsigned size = sizeof(FourMomenta_);
    size += FourMomenta_.size() * sizeof(TLorentzVector);
    std::cout << "  Size of FourMomenta_:         " << size << " byte\n";

    size = sizeof(HelicityAngles_);
    for (std::vector<double>& v : HelicityAngles_) {
        size += sizeof(v);
        size += v.size() * sizeof(double);
    }
    std::cout << "+ Size of HelicityAngles_:      " << size << " byte\n";
    totSize += size;

    size = sizeof(Data_);
    for (std::vector<std::vector<double> >& v : Data_) {
        size += sizeof(v);
        for (std::vector<double>& vv : v) {
            size += sizeof(vv);
            size += vv.size() * sizeof(double);
        }
    }
    std::cout << "+ Size of Data_:                " << size << " byte\n";
    totSize += size;

    size = sizeof(CachedAmplitudes_);
    for (std::vector<Amp>& v : CachedAmplitudes_) {
        size += sizeof(v);
        size += v.size() * sizeof(Amp);
    }
    std::cout << "+ Size of CachedAmplitudes_:    " << size << " byte\n";
    totSize += size;

    std::cout << "= Size of DataPoint:            " << totSize << " byte\n";
}

}

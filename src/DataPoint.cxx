#include "DataPoint.h"

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

    Data_.resize(dataAccessors.size());
    CachedAmplitudes_.resize(dataAccessors.size());
    CalculationStatuses_.resize(dataAccessors.size());

    for (DataAccessor* d : dataAccessors) {
        Data_.at(d->index()).resize(d->maxSymmetrizationIndex() + 1);
        CachedAmplitudes_.at(d->index()) = std::vector<Amp>(d->maxSymmetrizationIndex() + 1, 0);
        CalculationStatuses_.at(d->index()) = std::vector<CalculationStatus>(d->maxSymmetrizationIndex() + 1, kUncalculated);
        for (unsigned int symInd = 0; symInd < d->maxSymmetrizationIndex() + 1; ++symInd) {
            /// \todo size 1 ok?
            Data_.at(d->index()).at(symInd) = {0.};
        }
    }
}


}

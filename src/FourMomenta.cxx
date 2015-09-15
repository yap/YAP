#include "FourMomenta.h"

#include "DataPoint.h"
#include "logging.h"
#include "ParticleCombination.h"

namespace yap {

//-------------------------
FourMomenta::FourMomenta(InitialStateParticle* isp) :
    DataAccessor(isp, &ParticleCombination::equivByOrderlessContent),
    InitialStateIndex_(-1)
{
}

//-------------------------
int FourMomenta::findInitialStateParticle()
{
    if (InitialStateIndex_ < 0) {

        // count FSP
        unsigned fsp = 0;
        for (auto& kv : SymmetrizationIndices_)
            if (kv.first->indices().size() > fsp)
                fsp = kv.first->indices().size();

        // look for ISP
        for (auto& kv : SymmetrizationIndices_)
            if (kv.first->indices().size() == fsp) {
                InitialStateIndex_ = kv.second;
                break;
            }
    }

    if (InitialStateIndex_ < 0)
        LOG(ERROR) << "FourMomenta::findInitialStateParticle() - could not find InitialStateParticle index.";

    return InitialStateIndex_;
}

//-------------------------
bool FourMomenta::consistent() const
{
    bool result = true;

    // check that the first indices in the SymmetrizationIndices_ are the final state particles in order
    for (auto& kv : SymmetrizationIndices_)
        if (kv.first->isFinalStateParticle() and kv.first->indices()[0] != kv.second) {
            LOG(ERROR) << "FourMomenta::consistent - final-state particle id does not match index ("
                       << kv.first->indices()[0] << " != " << kv.second << ")";
            result = false;
        }

    if (InitialStateIndex_ < 0) {
        LOG(ERROR) << "FourMomenta::consistent - does not contain initial-state particle combination.";
        result = false;
    }

    result &= DataAccessor::consistent();

    return result;
}

//-------------------------
void FourMomenta::calculate(DataPoint& d)
{
    std::vector<CalculationStatus> calculationStatuses(SymmetrizationIndices_.size(), kUncalculated);

    for (auto& kv : SymmetrizationIndices_) {

        // check if calculation necessary
        if (calculationStatuses.at(kv.second) == kCalculated)
            continue;

        // if final state particle, 4-momentum already set; else
        if (!kv.first->isFinalStateParticle()) {
            // reset 4-momentum
            d.FourMomenta_.at(kv.second).SetXYZT(0, 0, 0, 0);

            // add in final-state particle momenta
            for (unsigned i : kv.first->indices())
                d.FourMomenta_.at(kv.second) += d.FourMomenta_.at(i);
        }

        double m2 = d.FourMomenta_.at(kv.second).M2();
        double m = sqrt(m2);

        std::vector<double>& D = data(d, kv.second);
        D = {m2, m};
        calculationStatuses.at(kv.second) = kCalculated;
    }
}

//-------------------------
const TLorentzVector& FourMomenta::p(const DataPoint& d, unsigned i) const
{
    return d.FourMomenta_.at(i);
}

}

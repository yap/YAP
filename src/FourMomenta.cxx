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
            if (kv.first->isFinalStateParticle())
                fsp += 1;

        // look for ISP
        for (auto& kv : SymmetrizationIndices_)
            if (kv.first->indices().size() == fsp)
                InitialStateIndex_ = kv.second;
    }

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
    // reserve space for four momenta
    d.FourMomenta_.resize(CalculationStatuses_.size());

    CalculationStatuses_.assign(CalculationStatuses_.size(), kUncalculated);

    for (auto& kv : SymmetrizationIndices_) {

        // check if calculation necessary
        if (CalculationStatuses_[kv.second] == kCalculated)
            continue;

        // if final state particle, 4-momentum already set; else
        if (!kv.first->isFinalStateParticle()) {
            // reset 4-momentum
            d.FourMomenta_[kv.second].SetXYZT(0, 0, 0, 0);

            // add in final-state particle momenta
            for (unsigned i = 0 : kv.first->indices())
                d.FourMomenta_[kv.second] += d.FourMomenta_[i];
        }

        double m2 = d.FourMomenta_[kv.second].M2();
        double m = sqrt(m2);

        std::vector<double>& D = data(d, kv.second);
        D = {m2, m};
        CalculationStatuses_[kv.second] = kCalculated;
    }
}

//-------------------------
const TLorentzVector& FourMomenta::p(const DataPoint& d, unsigned i) const
{
    return d.FourMomenta_[i];
}

}

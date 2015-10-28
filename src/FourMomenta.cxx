#include "FourMomenta.h"

#include "DataPoint.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "ParticleCombination.h"

namespace yap {

//-------------------------
FourMomenta::FourMomenta() :
    StaticDataAccessor(&ParticleCombination::equivByOrderlessContent),
    InitialStateIndex_(-1),
    M2_(this),
    M_(this)
{
}

//-------------------------
int FourMomenta::findInitialStateParticle()
{
    if (InitialStateIndex_ < 0) {

        // count FSP
        unsigned fsp = 0;
        for (auto& pc : particleCombinations())
            if (pc->indices().size() > fsp)
                fsp = pc->indices().size();

        // look for ISP
        for (auto& kv : symmetrizationIndices())
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
    for (auto& kv : symmetrizationIndices())
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
    // use a default dataPartitionIndex of 0

    M2_.setCalculationStatus(kUncalculated, 0);
    M_.setCalculationStatus(kUncalculated, 0);

    for (auto& kv : symmetrizationIndices()) {

        // check if calculation necessary
        if (M2_.calculationStatus(kv.first, kv.second, 0) == kCalculated)
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
        M2_.setValue(m2, d, kv.second, 0);
        M_.setValue(sqrt(m2), d, kv.second, 0);
    }
}


}

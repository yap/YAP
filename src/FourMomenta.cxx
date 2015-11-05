#include "FourMomenta.h"

#include "DataPoint.h"
#include "FinalStateParticle.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "ParticleCombination.h"

namespace yap {

//-------------------------
FourMomenta::FourMomenta() :
    StaticDataAccessor(&ParticleCombination::equivByOrderlessContent),
    InitialStatePC_(nullptr),
    M_(this)
{
}

//-------------------------
void FourMomenta::prepare()
{
    // count FSP particles
    unsigned fsp = 0;
    for (auto& pc : particleCombinations())
        if (pc->indices().size() > fsp)
            fsp = pc->indices().size();

    // look for ISP
    for (auto& kv : symmetrizationIndices())
        if (kv.first->indices().size() == fsp) {
            InitialStatePC_ = kv.first;

            // set fsp masses
            FinalStateParticleM_.resize(fsp);

            for (ParticleIndex i : kv.first->indices()) {
                for (auto& fsp : initialStateParticle()->finalStateParticles()) {
                    bool found(false);
                    for (auto& pc : fsp->particleCombinations()) {
                        ParticleIndex ii = pc->indices().at(0);

                        if (i == ii) {
                            FinalStateParticleM_[i] = fsp->mass();

                            //DEBUG("set mass for fsp " << unsigned(i) << " to " << m);

                            found = true;
                            break;
                        }
                    }

                    if (found)
                        break;
                }
            }
            // done setting fsp masses

            break;
        }


    if (!InitialStatePC_)
        LOG(ERROR) << "FourMomenta::findInitialStateParticle() - could not find InitialStateParticle.";
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

    if (!InitialStatePC_) {
        LOG(ERROR) << "FourMomenta::consistent - does not contain initial-state particle combination.";
        result = false;
    }

    if (FinalStateParticleM_.empty()) {
        LOG(ERROR) << "FourMomenta::consistent - FinalStateParticleM_ and FinalStateParticleM2_ have not been filled.";
        result = false;
    }

    result &= DataAccessor::consistent();

    return result;
}

//-------------------------
void FourMomenta::calculate(DataPoint& d)
{
    // use a default dataPartitionIndex of 0
    M_.setCalculationStatus(kUncalculated, 0);

    for (auto& kv : symmetrizationIndices()) {

        // check if calculation necessary
        if (M_.calculationStatus(kv.first, kv.second, 0) == kCalculated)
            continue;

        // reset 4-momentum
        d.FourMomenta_.at(kv.second).SetXYZT(0, 0, 0, 0);

        // add in final-state particle momenta
        for (unsigned i : kv.first->indices())
            d.FourMomenta_.at(kv.second) += d.FSPFourMomenta_.at(i);

        M_.setValue(d.FourMomenta_.at(kv.second).M(), d, kv.second, 0);
    }
}

//-------------------------
// std::vector<TLorentzVector> FourMomenta::calculateFourMomenta(const DataPoint& d) const
// {
//     // // calculate |p|
//     // std::vector<double> p(FinalStateParticleM_.size(), -1);
//     // for (unsigned i = 0; i < p.size(); ++i)
//     //     p[i] = ;
// }


}

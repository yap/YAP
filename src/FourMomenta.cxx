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
    InitialStatePC_ = nullptr;
    for (auto& kv : symmetrizationIndices())
        if (kv.first->indices().size() == fsp) {
            InitialStatePC_ = kv.first;
            break;
        }

    if (!InitialStatePC_)
        LOG(ERROR) << "FourMomenta::findInitialStateParticle() - could not find InitialStateParticle.";

    // set FSP masses, and FSP PCs
    FinalStateParticleM_.assign(fsp, nullptr);
    FinalStatePC_.assign(fsp, nullptr);
    for (ParticleIndex i = 0; i < fsp; ++i) {
        FinalStatePC_[i] = ParticleCombination::uniqueSharedPtr(i);
        // set FSP mass
        for (auto& fsp : initialStateParticle()->finalStateParticles()) {
            for (auto& pc : fsp->particleCombinations()) {
                if (i == pc->indices()[0]) {
                    FinalStateParticleM_[i] = fsp->mass();
                    //DEBUG("set mass for fsp " << unsigned(i) << " to " << m);
                    break;
                }
            }
            if (FinalStateParticleM_[i])
                break;
        }
    }

    // Set recoil PC's & pair PC's
    RecoilPC_.assign(fsp, nullptr);
    PairPC_.assign(fsp, ParticleCombinationVector(fsp, nullptr));
    for (ParticleIndex i = 0; (unsigned)i < fsp; ++i) {
        std::shared_ptr<const ParticleCombination> i_pc = ParticleCombination::uniqueSharedPtr(i);
        // build vector of other final state particles
        ParticleCombinationVector rec_pcs;
        for (ParticleIndex j = 0; (unsigned)j < fsp; ++j) {
            if (j == i)
                continue;
            std::shared_ptr<const ParticleCombination> j_pc = ParticleCombination::uniqueSharedPtr(j);
            rec_pcs.push_back(j_pc);
            // set pair pc and add to object
            PairPC_[i][j] = ParticleCombination::uniqueSharedPtr({i_pc, j_pc});
            addSymmetrizationIndex(PairPC_[i][j]);
        }
        // set recoil pc and add to object
        RecoilPC_[i] = ParticleCombination::uniqueSharedPtr(rec_pcs);
        addSymmetrizationIndex(RecoilPC_[i]);
    }
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
std::vector<TLorentzVector> FourMomenta::calculateFourMomenta(const DataPoint& d) const
{
    double M2 = m2(d, InitialStatePC_);

    std::vector<TLorentzVector> P(FinalStateParticleM_.size());
    // calculate |p|, setting each vector in z direction initially
    for (unsigned i = 0; i < P.size(); ++i) {
        double fsp_m = m(d, FinalStatePC_[i]);
        double recoil_mass = m(d, RecoilPC_[i]);
        double p = sqrt((M2 - pow(fsp_m + recoil_mass, 2)) * (M2 - pow(fsp_m - recoil_mass, 2)) / 4 / M2);
        P[i].SetXYZM(0, 0, p, fsp_m);
    }

    // calculate angles
    std::vector<std::vector<double> > angle(P.size(), std::vector<double>(P.size(), 0));
    for (unsigned i = 0; i < P.size(); ++i) {
        for (unsigned j = i + 1; j < P.size(); ++j) {
            angle[i][j] = acos(P[i].E() * P[j].E() + 0.5 * (pow(m(d, FinalStatePC_[i]), 2) + pow(m(d, FinalStatePC_[j]), 2) - m2(d, PairPC_[i][j])));
            angle[j][i] = angle[i][j];
        }
    }

    // leave P[0] aligned with z axis
    // rotate remaining by relative angles around y axis
    /// \todo Make work for 4 or more particles
    for (unsigned i = 1; i < P.size(); ++i) {
        P[i].RotateY(angle[i][i - 1]);
    }

    return P;
}


}

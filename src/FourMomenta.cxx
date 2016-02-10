#include "FourMomenta.h"

#include "container_utils.h"
#include "DataPoint.h"
#include "DecayingParticle.h"
#include "Exceptions.h"
#include "FinalStateParticle.h"
#include "FourVector.h"
#include "Model.h"
#include "logging.h"
#include "MathUtilities.h"
#include "ParticleCombination.h"
#include "ParticleCombinationCache.h"
#include "ThreeVector.h"

#include <iomanip>

namespace yap {

//-------------------------
FourMomenta::FourMomenta(Model* m) :
    StaticDataAccessor(m, &ParticleCombination::equivByOrderlessContent),
    M_(RealCachedDataValue::create(this))
{
    if (!model())
        throw exceptions::Exception("Model unset", "FourMomenta::FourMomenta");
}

//-------------------------
void FourMomenta::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    StaticDataAccessor::addParticleCombination(pc);
    // check for ISP
    if (!InitialStatePC_ and pc->indices().size() == model()->finalStateParticles().size())
        InitialStatePC_ = pc;
}

//-------------------------
bool FourMomenta::consistent() const
{
    bool C = StaticDataAccessor::consistent();

    // check that the first indices in the SymmetrizationIndices_ are the final state particles in order
    for (auto& kv : symmetrizationIndices())
        if (kv.first->isFinalStateParticle() and kv.first->indices()[0] != kv.second) {
            FLOG(ERROR) << "final-state particle id does not match index (" << kv.first->indices()[0] << " != " << kv.second << ")";
            C &= false;
        }

    if (!InitialStatePC_) {
        FLOG(ERROR) << "ISP ParticleCombination has not been recorded.";
        C &= false;
    }

    return C;
}

//-------------------------
const FourVector<double>& FourMomenta::p(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
{
    if (pc->isFinalStateParticle())
        return d.finalStateFourMomenta()[pc->indices()[0]];
    return d.fourMomenta()[symmetrizationIndex(pc)];
}

//-------------------------
void FourMomenta::calculate(DataPoint& d, unsigned dataPartitionIndex)
{
    M_->setCalculationStatus(kUncalculated, dataPartitionIndex);

    for (auto& kv : symmetrizationIndices()) {

        // check if calculation necessary
        if (M_->calculationStatus(kv.first, kv.second, dataPartitionIndex) == kCalculated)
            continue;

        const auto P = std::accumulate(kv.first->indices().begin(), kv.first->indices().end(), FourVector_0,
        [&](const FourVector<double>& V, unsigned i) {return V + d.FSPFourMomenta_.at(i);});

        d.FourMomenta_.at(kv.second) = P;

        FDEBUG(*kv.first << " = " << d.FourMomenta_.at(kv.second));

        M_->setValue(abs(d.FourMomenta_.at(kv.second)), d, kv.second, dataPartitionIndex);

        // DEBUG("FourMomenta::calculate - 4-momentum " << * (kv.first) << ": " << to_string(d.FourMomenta_.at(kv.second)) );
        // DEBUG("FourMomenta::calculate - Set mass for " << * (kv.first) << " to " << M_->value(d, kv.second));
    }
}

//-------------------------
double FourMomenta::m(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
{
    if (pc->isFinalStateParticle())
        return model()->finalStateParticles()[pc->indices()[0]]->mass()->value();
    return M_->value(d, symmetrizationIndex(pc));
}

//-------------------------
std::ostream& FourMomenta::printMasses(const DataPoint& d, std::ostream& os) const
{
    os << "Invariant masses:" << std::endl;
    unsigned n_fsp = model()->finalStateParticles().size();
    unsigned n = n_fsp + 2;
    unsigned m_p = 6;

    std::set<unsigned> used;

    // print the ISP
    for (auto& kv : symmetrizationIndices())
        // if ISP and not yet printed (should only print once)
        if (kv.first->indices().size() == n_fsp and used.find(kv.second) == used.end()) {
            os << "    ISP : " << std::setw(n) << indices_string(*kv.first)
               << " = " << std::setprecision(m_p) << m(d, kv.first) << " GeV/c^2"
               << " (nominally " << std::setprecision(m_p) << model()->initialStateParticle()->mass()->value() << " GeV/c^2)" << std::endl;
            used.insert(kv.second);
        }

    // print the FSP's
    for (size_t i = 0; i < d.finalStateFourMomenta().size(); ++i)
        os << "    FSP : " << std::setw(n) << (std::string("(") + std::to_string(i) + ")")
           << " = " << std::setprecision(m_p) << abs(d.finalStateFourMomenta()[i]) << " GeV/c^2"
           << " (nominally " << std::setprecision(m_p) << model()->finalStateParticles()[i]->mass()->value() << " GeV/c^2)" << std::endl;

    // print the rest in increasing number of particle content
    for (unsigned i = 2; i < n_fsp; ++i)
        for (auto& kv : symmetrizationIndices())
            // if i-particle mass and not yet printed
            if (kv.first->indices().size() == i and used.find(kv.second) == used.end()) {
                os << "    " << i << "-p : " << std::setw(n) << indices_string(*kv.first)
                   << " = " << std::setprecision(m_p) << m(d, kv.first) << " GeV/c^2" << std::endl;
                used.insert(kv.second);
            }

    // print unused
    for (auto& kv : symmetrizationIndices())
        if (used.find(kv.second) == used.end()) {
            os << "        : " << std::setw(n) << indices_string(*kv.first) << " = " << m(d, kv.first) << " GeV/c^2" << std::endl;
            used.insert(kv.second);
        }

    return os;
}

}

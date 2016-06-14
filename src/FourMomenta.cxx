#include "FourMomenta.h"

#include "CachedDataValue.h"
#include "CalculationStatus.h"
#include "container_utils.h"
#include "DecayingParticle.h"
#include "Exceptions.h"
#include "FinalStateParticle.h"
#include "FourVector.h"
#include "logging.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleCombination.h"
#include "StatusManager.h"

#include <iomanip>

namespace yap {

//-------------------------
FourMomenta::FourMomenta(Model* m) :
    StaticDataAccessor(m, ParticleCombination::equalByOrderlessContent),
    ISPIndex_(-1),
    P_(FourVectorCachedDataValue::create(this)),
    M_(RealCachedDataValue::create(this, {}, {P_}))
{
    if (!model())
        throw exceptions::Exception("Model unset", "FourMomenta::FourMomenta");
}

//-------------------------
void FourMomenta::setFinalStateMomenta(DataPoint& d, const std::vector<FourVector<double> >& P, StatusManager& sm) const
{
    if (P.size() != FSPIndices_.size())
        throw exceptions::Exception("Wrong number of momenta provided (" + std::to_string(P.size()) + " != " + std::to_string(FSPIndices_.size()) + ")",
                                    "FourMomenta::setFinalStateMomenta");

    for (size_t i = 0; i < P.size(); ++i)
        P_->setValue(P[i], d, FSPIndices_[i], sm);
}

//-------------------------
void FourMomenta::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    StaticDataAccessor::addParticleCombination(pc);
    auto index = symmetrizationIndex(pc);

    // check for ISP
    if (ISPIndex_ < 0 and pc->indices().size() == model()->finalStateParticles().size())
        ISPIndex_ = index;

    /// check for FSP
    if (pc->isFinalStateParticle()) {
        if (pc->indices()[0] + 1 > FSPIndices_.size())
            FSPIndices_.resize(pc->indices()[0] + 1, -1);
        if (FSPIndices_[pc->indices()[0]] < 0)
            FSPIndices_[pc->indices()[0]] = index;
    }
}

//-------------------------
bool FourMomenta::consistent() const
{
    bool C = StaticDataAccessor::consistent();

    if (ISPIndex_ < 0) {
        FLOG(ERROR) << "ISP symmetrization index has not been recorded.";
        C &= false;
    }

    if (FSPIndices_.size() != model()->finalStateParticles().size() or std::any_of(FSPIndices_.begin(), FSPIndices_.end(), [](int i) {return i < 0;})) {
        FLOG(ERROR) << "FSP symmetrization indices not all recorded";
        C &= false;
    }

    return C;
}

//-------------------------
const FourVector<double> FourMomenta::totalMomentum(const DataPoint& d) const
{
    if (ISPIndex_ < 0)
        throw exceptions::Exception("Initial-state particle unknown", "FourMomenta::totalMomentum");
    return P_->value(d, ISPIndex_);
}

//-------------------------
const std::vector<FourVector<double> > FourMomenta::finalStateMomenta(const DataPoint& d) const
{
    std::vector<FourVector<double> > P;
    P.reserve(FSPIndices_.size());

    for (size_t i = 0; i  < FSPIndices_.size(); ++i) {
        if (FSPIndices_[i] < 0)
            throw exceptions::Exception("Final-state particle " + std::to_string(i) + "unknown", "FourMomenta::finalStateFourMomenta");
        P.push_back(P_->value(d, FSPIndices_[i]));
    }
    return P;
}

//-------------------------
FourVector<double> FourMomenta::p(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
{
    return P_->value(d, symmetrizationIndex(pc));
}

double FourMomenta::m(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
{
    return M_->value(d, symmetrizationIndex(pc));
}

//-------------------------
void FourMomenta::calculate(DataPoint& d, StatusManager& sm) const
{
    // set all masses as uncalculated
    sm.set(*M_, CalculationStatus::uncalculated);

    // get fsp four momenta:
    auto fsp = finalStateMomenta(d);

    for (auto& kv : symmetrizationIndices()) {

        // check if calculation unnecessary
        if (sm.status(*M_, kv.second) == CalculationStatus::calculated)
            continue;

        const auto P = std::accumulate(kv.first->indices().begin(), kv.first->indices().end(), FourVector_0,
        [&](const FourVector<double>& V, unsigned i) {return V + fsp[i];});

        P_->setValue(P, d, kv.second, sm);
        M_->setValue(abs(P), d, kv.second, sm);
    }
}

//-------------------------
std::ostream& print_mp_string(std::ostream& os, unsigned n, unsigned m_p, std::shared_ptr<ParticleCombination> pc, double m, FourVector<double> p, std::shared_ptr<RealParameter> M = nullptr)
{
    os << std::setw(n) << indices_string(*pc) << " : "
       << "m = " << std::setprecision(m_p) << m << " GeV/c^2";
    if (M)
        os << " (nominally " << std::setprecision(m_p) << M->value() << " GeV/c^2)";
    else
        os << "            " << std::setw(m_p) << " "                << "         ";
    os << "\tp = " << p << " GeV";
    return os;
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
            os << "    ISP : ";
            print_mp_string(os, n, m_p, kv.first, m(d, kv.first), p(d, kv.first));
            os << std::endl;
            used.insert(kv.second);
        }

    // print the FSP's
    for (size_t i = 0; i < FSPIndices_.size(); ++i)
        for (auto& kv : symmetrizationIndices())
            if (kv.first->isFinalStateParticle() and kv.first->indices()[0] == i and used.find(kv.second) == used.end()) {
                os << "    FSP : ";
                print_mp_string(os, n, m_p, kv.first, m(d, kv.first), p(d, kv.first), model()->finalStateParticles()[i]->mass());
                os << std::endl;
                used.insert(kv.second);
            }


    // print the rest in increasing number of particle content
    for (unsigned i = 2; i < n_fsp; ++i)
        for (auto& kv : symmetrizationIndices())
            // if i-particle mass and not yet printed
            if (kv.first->indices().size() == i and used.find(kv.second) == used.end()) {
                os << "    " << i << "-p : ";
                print_mp_string(os, n, m_p, kv.first, m(d, kv.first), p(d, kv.first));
                os << std::endl;
                used.insert(kv.second);
            }

    // print unused
    for (auto& kv : symmetrizationIndices())
        if (used.find(kv.second) == used.end()) {
            os << "        : ";
            print_mp_string(os, n, m_p, kv.first, m(d, kv.first), p(d, kv.first));
            os << std::endl;
            used.insert(kv.second);
        }

    return os;
}

}

#include "FourMomenta.h"

#include "container_utils.h"
#include "DataPoint.h"
#include "Exceptions.h"
#include "FinalStateParticle.h"
#include "FourVector.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "MathUtilities.h"
#include "ParticleCombination.h"
#include "ParticleCombinationCache.h"
#include "ThreeVector.h"

namespace yap {

//-------------------------
FourMomenta::FourMomenta(InitialStateParticle* isp) :
    StaticDataAccessor(isp, &ParticleCombination::equivByOrderlessContent),
    M_(std::make_shared<RealCachedDataValue>(this))
{
    if (!initialStateParticle())
        throw exceptions::Exception("InitialStateParticle unset", "FourMomenta::FourMomenta");
}

//-------------------------
void FourMomenta::prepare()
{
    ParticleCombinationVector PCs = particleCombinations();

    // count FSP particles
    unsigned fsp = 0;
    for (auto& pc : PCs)
        if (pc->indices().size() > fsp)
            fsp = pc->indices().size();

    // look for ISP
    auto it = std::find_if(PCs.begin(), PCs.end(),
    [&](const ParticleCombinationVector::value_type & a) {return a->indices().size() == fsp;});

    if (it == PCs.end())
        LOG(ERROR) << "FourMomenta::findInitialStateParticle() - could not find InitialStateParticle.";

    InitialStatePC_ = *it;

    // set FSP masses, and FSP PCs
    FinalStateParticleM_.assign(fsp, nullptr);
    FinalStatePC_.assign(fsp, nullptr);
    for (ParticleIndex i = 0; i < fsp; ++i) {
        FinalStatePC_[i] = initialStateParticle()->particleCombinationCache().fsp(i);
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
    RecoilPC_.assign(FinalStatePC_.size(), nullptr);
    PairPC_.assign(FinalStatePC_.size(), ParticleCombinationVector(FinalStatePC_.size(), nullptr));
    for (ParticleIndex i = 0; (unsigned)i < fsp; ++i) {
        auto i_pc = FinalStatePC_[i];
        // build vector of other final state particles
        ParticleCombinationVector rec_pcs;
        for (ParticleIndex j = 0; (unsigned)j < fsp; ++j) {
            if (j == i)
                continue;
            auto j_pc = FinalStatePC_[j];
            rec_pcs.push_back(j_pc);
            // set pair pc and add to object
            PairPC_[i][j] = initialStateParticle()->particleCombinationCache().composite({i_pc, j_pc});
            addParticleCombination(PairPC_[i][j]);
        }
        // set recoil pc and add to object
        RecoilPC_[i] = initialStateParticle()->particleCombinationCache().composite(rec_pcs);
        addParticleCombination(RecoilPC_[i]);
    }
}

//-------------------------
bool FourMomenta::consistent() const
{
    bool C = StaticDataAccessor::consistent();

    // check that the first indices in the SymmetrizationIndices_ are the final state particles in order
    for (auto& kv : symmetrizationIndices())
        if (kv.first->isFinalStateParticle() and kv.first->indices()[0] != kv.second) {
            FLOG(ERROR) << "final-state particle id does not match index ("
                        << kv.first->indices()[0] << " != " << kv.second << ")";
            C &= false;
        }

    if (!InitialStatePC_) {
        FLOG(ERROR) << "does not contain initial-state particle combination.";
        C &= false;
    }

    if (FinalStateParticleM_.empty()) {
        FLOG(ERROR) << "FinalStateParticleM_ and FinalStateParticleM2_ have not been filled.";
        C &= false;
    }

    return C;
}

//-------------------------
void FourMomenta::calculate(DataPoint& d)
{
    // use a default dataPartitionIndex of 0
    M_->setCalculationStatus(kUncalculated, 0);

    for (auto& kv : symmetrizationIndices()) {

        // check if calculation necessary
        if (M_->calculationStatus(kv.first, kv.second, 0) == kCalculated)
            continue;

        // reset 4-momentum
        d.FourMomenta_.at(kv.second) = {0, 0, 0, 0};

        // add in final-state particle momenta
        for (unsigned i : kv.first->indices())
            d.FourMomenta_.at(kv.second) += d.FSPFourMomenta_.at(i);

        M_->setValue(abs(d.FourMomenta_.at(kv.second)), d, kv.second, 0);

        DEBUG("FourMomenta::calculate - 4-momentum " << * (kv.first) << ": " << to_string(d.FourMomenta_.at(kv.second)) );
        DEBUG("FourMomenta::calculate - Set mass for " << * (kv.first) << " to " << M_->value(d, kv.second));
    }
}

//-------------------------
std::vector<FourVector<double> > FourMomenta::calculateFourMomenta(const DataPoint& d) const
{
    double M2 = m2(d, InitialStatePC_);

    std::vector<FourVector<double> > P(FinalStateParticleM_.size());
    assert(FinalStateParticleM_.size() == FinalStatePC_.size());
    std::vector<double> fsp_m2(P.size(), 0);
    // calculate E and |p| for each initial state particle, storing in Z-aligned FourVector
    for (unsigned i = 0; i < P.size(); ++i) {
        fsp_m2[i] = m2(d, FinalStatePC_[i]);
        double E = (M2 - m2(d, RecoilPC_[i]) + fsp_m2[i]) / 2 / sqrt(M2);
        P[i] = {E, 0, 0, sqrt(pow(E, 2) - fsp_m2[i])};
    }

    // debug
    DEBUG("E:");
    for (auto p : P) {
        DEBUG(to_string(p));
    }

    // if only particle (should not happen)
    if (P.size() < 2)
        return P;

    // if only two particles (should not happen)
    if (P.size() < 3) {
        P[1][3] *= -1;
        return P;
    }

    // store angles between particle i and particles 0, 1, 2
    std::vector<std::array<double, 3>> cosAngle(P.size(), {0, 0, 0});
    for (size_t i = 0; i < P.size(); ++i)
        for (size_t j = 0; j < 3; ++j)
            if (j != i)
                // cos(theta_ij) = (E_i * E_j - 1/2 (m^2_ij - m2_i - m2_j)) / |p_i| / |p_j|
                cosAngle[i][j] = (P[i][0] * P[j][0] - 0.5 * (m2(d, PairPC_[i][j]) - fsp_m2[i] - fsp_m2[j])) / P[i][3] / P[j][3];

    double sin01 = sqrt(1 - pow(cosAngle[0][1], 2));

    // P[0] defined parallel to z

    // P[1] defined to be in the x-z plane
    P[1] = {P[1][0], P[1][3]* sin01, 0, P[1][3]* cosAngle[0][1]};

    // define P[2] to be in +Y direction
    auto v2 = ThreeVector<double>({(cosAngle[2][1] - cosAngle[2][0] * cosAngle[0][1]) / sin01, 0, cosAngle[2][0]});
    v2[2] = sqrt(1 - yap::norm(v2));
    P[2] = FourVector<double>(P[2][0], v2 * P[2][3]);

    // define remaining 4-momenta
    for (unsigned i = 3; i < P.size(); ++i) {
        auto vi = ThreeVector<double>({(cosAngle[i][1] - cosAngle[i][0] * cosAngle[0][1]) / sin01, 0, cosAngle[i][0]});
        vi[2] = (yap::abs(v2) * cosAngle[i][2] - vi * v2) / v2[0];
        P[i] = FourVector<double>(P[i][0], P[i][3] * vi);
    }


    // debug
    DEBUG("final:");
    for (auto p : P) {
        DEBUG(to_string(p));
    }

    return P;
}

//-------------------------
bool FourMomenta::calculateMissingMasses(DataPoint& d)
{
    /// for n finalStateParticles
    /// m_{1..n}^2 = (2-n) \sum_{k=0}^n m_k^2 + \sum_{k=0}^n \sum_{l=k+1}^n m_{kl}^2

    // set initial State particle m
    M_->setValue(initialStateParticle()->mass()->value(), d, symmetrizationIndex(InitialStatePC_), 0u);
    double M2 = m2(d, InitialStatePC_);
    unsigned n = InitialStatePC_->indices().size();

    ParticleCombinationVector pairPCs = pairParticleCombinations();

    if (pairPCs.size() != n * (n - 1) / 2) {
        LOG(ERROR) << "Wrong number of pair PCs.";
        return false;
    }

    // check if we have enough masses to calculate the rest
    unsigned nUnset(0);
    unsigned iUnset(0);
    for (unsigned i = 0; i < pairPCs.size(); ++i) {
        if (m(d, pairPCs[i]) < 0.) {
            nUnset++;
            iUnset = i;
        }
    }

    if (nUnset > 1) {
        LOG(ERROR) << "Cannot calculate masses, not enough two-particle invariant masses set.";
        return false;
    }

    ///
    /// calculate unset pair mass
    ///

    double m2_ab(0);
    for (auto& pc : FinalStatePC_) {
        m2_ab += m2(d, pc);
    }
    m2_ab *= 2 - int(n);
    m2_ab = M2 - m2_ab;

    for (unsigned i = 0; i < pairPCs.size(); ++i) {
        if (i == iUnset)
            continue;
        m2_ab -= m2(d, pairPCs[i]);
    }

    if (m2_ab < 0) {
        LOG(ERROR) << "FourMomenta::calculateMissingMasses : Resulting two-particle m2 is < 0.";
        return false;
    }

    M_->setValue(sqrt(m2_ab), d, symmetrizationIndex(pairPCs[iUnset]), 0u);

    DEBUG("calculated missing mass for " << *pairPCs[iUnset] << " = " << M_->value(d, symmetrizationIndex(pairPCs[iUnset])));


    /// for a 3 particle system, we are done
    if (n < 4) {
        //debug
        for (auto& pc : RecoilPC_) {
            DEBUG("recoil mass for " << *pc << " = " << M_->value(d, symmetrizationIndex(pc)));
        }


        return true;
    }


    /// calculate recoil masses
    for (auto& pc : RecoilPC_) {
        double m2_recoil(0);
        for (auto& i : pc->indices()) {
            m2_recoil += pow(FinalStateParticleM_.at(i)->value(), 2);
        }

        for (auto& pcPair : pairPCs) {
            if (contains(pc->indices(), pcPair->indices()))
                m2_recoil += m2(d, pcPair);
        }

        if (m2_recoil < 0) {
            LOG(ERROR) << "Resulting recoil m2 is < 0.";
            return false;
        }

        M_->setValue(sqrt(m2_recoil), d, symmetrizationIndex(pc), 0u);

        DEBUG("calculated missing recoil mass for " << *pc << " = " << M_->value(d, symmetrizationIndex(pc)));
    }


    return true;
}

//-------------------------
ParticleCombinationVector FourMomenta::getDalitzAxes(std::vector<std::vector<ParticleIndex> > pcs) const
{
    ParticleCombinationVector M;

    ParticleCombinationVector PCV = particleCombinations();

    for (auto& v : pcs) {
        auto A = initialStateParticle()->particleCombinationCache().findByUnorderedContent(v);
        if (A.expired())
            throw exceptions::Exception("ParticleCombination not found", "FourMomenta::getDalitzAxes");
        for (auto& B : PCV)
            if (ParticleCombination::equivByOrderlessContent(A.lock(), B)) {
                M.push_back(B);
                break;
            }
    }

    // Check if all combinations found
    if (M.size() != pcs.size())
        throw exceptions::Exception("ParticleCombination not found", "FourMomenta::getDalitzAxes");

    return M;
}

//-------------------------
ParticleCombinationMap<double> FourMomenta::pairMasses(const DataPoint& d) const
{
    ParticleCombinationMap<double> result;

    for (auto& pc : pairParticleCombinations())
        result[pc] = m(d, pc);

    return result;
}

//-------------------------
ParticleCombinationMap<double> FourMomenta::pairMassSquares(const DataPoint& d) const
{
    ParticleCombinationMap<double> result;

    for (auto& pc : pairParticleCombinations())
        result[pc] = m2(d, pc);

    return result;
}

//-------------------------
void FourMomenta::setMasses(DataPoint& d, const ParticleCombinationVector& axes, const std::vector<double>& masses)
{
    if (axes.size() != masses.size()) {
        FLOG(ERROR) << "axes and masses vectors do not match in size.";
        throw exceptions::Exception("Masses do not match Dalitz axes", "FourMomenta::setMasses");
    }

    // reset all masses to -1
    resetMasses(d);

    // set independent masses according to axes
    for (unsigned i = 0; i < axes.size(); ++i)
        M_->setValue(masses[i], d, symmetrizationIndex(axes[i]), 0u);

    // recalculate dependent masses
    if (!calculateMissingMasses(d)) {
        // not enough masses set or not in phasespace
        resetMasses(d);
        /// \todo replace with specific class for outside-of-phase-space exception if that's the case
        throw exceptions::Exception("Failed to calculate missing masses", "FourMomenta::setMasses");
    }

    // recalculate final state masses
    d.setFinalStateFourMomenta(calculateFourMomenta(d));
}

//-------------------------
void FourMomenta::setSquaredMasses(DataPoint& d, const ParticleCombinationVector& axes, const std::vector<double>& squaredMasses)
{
    std::vector<double> masses(squaredMasses.size(), 0);

    std::cout << "mass squares \n";
    for (auto v : squaredMasses)
        std::cout << "  " << v << "\n";

    std::transform(squaredMasses.begin(), squaredMasses.end(), masses.begin(), sqrt);

    std::cout << "masses \n";
    for (auto v : masses)
        std::cout << "  " << v << "\n";

    setMasses(d, axes, masses);
}

//-------------------------
void FourMomenta::setMasses(DataPoint& d, ParticleCombinationMap<double> m)
{
    resetMasses(d);

    for (auto& kv : m)
        M_->setValue(kv.second, d, symmetrizationIndex(kv.first), 0u);

    // recalculate stuff
    if (!calculateMissingMasses(d)) {
        // not enough masses set or not in phasespace
        resetMasses(d);
        /// \todo replace with specific class for outside-of-phase-space exception if that's the case
        throw exceptions::Exception("Failed to calculate missing masses", "FourMomenta::setMasses");
    }

    d.setFinalStateFourMomenta(calculateFourMomenta(d));
}

//-------------------------
void FourMomenta::setMassSquares(DataPoint& d, ParticleCombinationMap<double> m2)
{
    for (auto& kv : m2)
        kv.second = sqrt(kv.second);

    setMasses(d, m2);
}

//-------------------------
void FourMomenta::printMasses(const DataPoint& d) const
{
    std::cout << "Invariant masses:\n";
    std::cout << "  " << *InitialStatePC_ << ": \t" << m(d, InitialStatePC_) << " GeV\n";

    if (InitialStatePC_->indices().size() > 3)
        for (auto& pc : RecoilPC_)
            std::cout << "  " << *pc << ": \t" << m(d, pc) << " GeV\n";

    for (auto& pc : pairParticleCombinations())
        std::cout << "  " << *pc << ": \t" << m(d, pc) << " GeV\n";

    for (auto& pc : FinalStatePC_)
        std::cout << "  " << *pc << ": \t" << m(d, pc) << " GeV\n";

    for (int i = 0; i <= maxSymmetrizationIndex(); ++i) {
        std::cout << "  symIndex " << i << ": " << M_->value(d, i) << " GeV\n";
    }
}

//-------------------------
void FourMomenta::resetMasses(DataPoint& d)
{
    for (int i = 0; i <= maxSymmetrizationIndex(); ++i) {
        if (i == int(symmetrizationIndex(InitialStatePC_)))
            continue;

        M_->setValue(-1, d, unsigned(i), 0u);
    }
}

//-------------------------
ParticleCombinationVector FourMomenta::pairParticleCombinations() const
{
    ParticleCombinationVector pairPCs;
    unsigned i(0);
    for (auto& v : PairPC_) {
        for (unsigned j = i++; j < v.size(); ++j)
            if (v[j])
                pairPCs.push_back(v[j]);
    }

    return pairPCs;
}

}

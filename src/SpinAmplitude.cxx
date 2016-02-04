#include "SpinAmplitude.h"

#include "ClebschGordan.h"
#include "Exceptions.h"
#include "logging.h"

namespace yap {

//-------------------------
SpinAmplitude::SpinAmplitude(const QuantumNumbers& initial,
                             const QuantumNumbers& final1,
                             const QuantumNumbers& final2,
                             unsigned l, unsigned two_s,
                             InitialStateParticle* isp) :
    StaticDataAccessor(),
    InitialTwoJ_(initial.twoJ()),
    FinalTwoJ_( {final1.twoJ(), final2.twoJ()}),
            L_(l),
            TwoS_(two_s)
{
    // check JLS triangle
    if (!triangle(InitialTwoJ_, 2 * L_, TwoS_))
        throw exceptions::AngularMomentumNotConserved("SpinAmplitude::SpinAmplitude");

    // check j1j2S triangle
    if (!triangle(FinalTwoJ_[0], FinalTwoJ_[1], TwoS_))
        throw exceptions::AngularMomentumNotConserved("SpinAmplitude::SpinAmplitude");

    // if (!conserves(InitialTwoJ_, FinalTwoJ_[0], FinalTwoJ_[1], l))
    //     throw exceptions::AngularMomentumNotConserved();

    // check charge conservation
    // \todo check somewhere else
    /*if (InitialQuantumNumbers_.Q() != FinalQuantumNumbers_[0].Q() + FinalQuantumNumbers_[1].Q())
        throw exceptions::Exception(std::string("charge conservation violated: ")
                                    + "(" + std::to_string(InitialQuantumNumbers_.Q())  + ") -> "
                                    + "(" + std::to_string(FinalQuantumNumbers_[0].Q()) + ") + "
                                    + "(" + std::to_string(FinalQuantumNumbers_[1].Q()) + ")",
                                    "SpinAmplitude::SpinAmplitude");*/
}

//-------------------------
void SpinAmplitude::calculate(DataPoint& d, unsigned dataPartitionIndex)
{
    // set calculation statuses uncalc'ed
    for (auto& a : amplitudeSet())
        a->setCalculationStatus(kUncalculated, dataPartitionIndex);

    // loop over particle combinations
    for (auto& pc : particleCombinations()) {

        unsigned symIndex = symmetrizationIndex(pc);

        // loop over mapping of parent spin projection to AmplitudeSubmap
        for (auto& aM_kv : Amplitudes_)
            // loop over mappin of daughter spin projection pairs to amplitudes
            for (auto& aSM_kv : aM_kv.second)
                // if yet uncalculated
                if (aSM_kv.second->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated) {

                    const auto& two_M = aM_kv.first; // parent spin projection
                    const auto& spp = aSM_kv.first; // SpinProjectionPair of daughters

                    aSM_kv.second->setValue(calc(two_M, spp[0], spp[1], d, pc), d, symIndex, dataPartitionIndex);

                }
    }
}

//-------------------------
SpinAmplitude::operator std::string() const
{
    std::string s = spin_to_string(InitialTwoJ_) + " -> ";
    for (auto& j : FinalTwoJ_)
        s += spin_to_string(j) + " + ";
    s.erase(s.size() - 2, 2);
    s += "with L = " + std::to_string(L_);
    s += " and S = " + spin_to_string(TwoS_);
    return s;
}

//-------------------------
std::set<int> SpinAmplitude::twoM() const
{
    std::set<int> S;
    // loop over amplitudes, key = two_M
    for (auto& kv : Amplitudes_)
        S.insert(kv.first);  // first entry is (twice) parent spin projection
    return S;
}

//-------------------------
void SpinAmplitude::addAmplitude(int two_M, int two_m1, int two_m2)
{
    // retrieve (or create) AmplitudeSubmap for two_M
    auto& ASM = Amplitudes_[two_M];

    SpinProjectionPair m1m2 = {two_m1, two_m2};
    // look for m1m2 in AmplitudeSubmap
    if (ASM.find(m1m2) != ASM.end())
        throw exceptions::Exception("Amplitude already stored for " + spin_to_string(two_M)
                                    + " -> " + spin_to_string(two_m1) + " + " + spin_to_string(two_m2),
                                    "SpinAmplitude::addAmplitude");

    // FDEBUG("adding CachedDataValue for")
    ASM[m1m2] = ComplexCachedDataValue::create(this);
}

//-------------------------
bool SpinAmplitude::equiv(const SpinAmplitude& B) const
{
    return // compare only spin of QuantumNumbers
        InitialTwoJ_ == B.InitialTwoJ_
        and FinalTwoJ_[0] == B.FinalTwoJ_[0]
        and FinalTwoJ_[1] == B.FinalTwoJ_[1]
        and L_ == B.L_
        and TwoS_ == B.TwoS_;
}

//-------------------------
bool SpinAmplitude::equals(const SpinAmplitude& B) const
{
    return symmetrizationIndices() == B.symmetrizationIndices()
           and equiv(B);
}

//-------------------------
CachedDataValueSet SpinAmplitude::amplitudeSet()
{
    CachedDataValueSet V;
    // loop over mapping of parent spin projection to AmplitudeSubmap
    for (auto& aM_kv : Amplitudes_)
        // loop over mappin of daughter spin projection pairs to amplitudes
        for (auto& aSM_kv : aM_kv.second)
            V.insert(aSM_kv.second);
    return V;
}

//-------------------------
std::string to_string(const SpinAmplitudeVector& saV)
{
    if (saV.empty())
        return std::string();

    std::string s = spin_to_string(saV[0]->initialTwoJ()) + " -> ";
    for (auto& j : saV[0]->finalTwoJ())
        s += spin_to_string(j) + " + ";
    s.erase(s.size() - 2, 2);
    s += "with LS =";
    for (auto& sa : saV)
        s += " (" + std::to_string(sa->L()) + ", " + spin_to_string(sa->twoS()) + "),";
    s.erase(s.size() - 1, 1);
    if (saV[0]->formalism().empty())
        return s;
    s += " in " + saV[0]->formalism();
    return s;
}

}

#include "SpinAmplitude.h"

#include "CachedValue.h"
#include "CalculationStatus.h"
#include "Exceptions.h"
#include "Spin.h"
#include "StatusManager.h"

namespace yap {

//-------------------------
SpinAmplitude::SpinAmplitude(Model& m, unsigned two_J, const SpinVector& two_j, unsigned l, unsigned two_s, const ParticleCombinationEqualTo& equal) :
    StaticDataAccessor(m, equal),
    InitialTwoJ_(two_J),
    FinalTwoJ_(two_j),
    L_(l),
    TwoS_(two_s)
{
    if (FinalTwoJ_.size() < 2)
        throw exceptions::Exception("Too few daughter spins specified (" + std::to_string(FinalTwoJ_.size()) + " < 2)",
                                    "SpinAmplitude:SpinAmplitude");
    // check JLS triangle
    if (!triangle(InitialTwoJ_, 2 * L_, TwoS_))
        throw exceptions::AngularMomentumNotConserved("SpinAmplitude::SpinAmplitude");
}

//-------------------------
void SpinAmplitude::calculate(DataPoint& d, StatusManager& sm) const
{
    // set all amplitudes Uncalculated
    sm.set(*this, CalculationStatus::uncalculated);

    // loop over particle combinations -> indices
    for (auto& pc_i : symmetrizationIndices()) {

        // loop over mapping of parent spin projection to AmplitudeSubmap
        for (auto& aM_kv : Amplitudes_) {

            const auto& two_M = aM_kv.first; // parent spin projection

            // loop over mappin of daughter spin projection pairs to amplitudes
            for (auto& aSM_kv : aM_kv.second)

                // if yet uncalculated
                if (sm.status(*aSM_kv.second, pc_i.second) == CalculationStatus::uncalculated)
                    aSM_kv.second->setValue(calc(two_M, aSM_kv.first, d, sm, pc_i.first), d, pc_i.second, sm);
        }
    }
}

//-------------------------
SpinAmplitude::operator std::string() const
{
    return spin_to_string(InitialTwoJ_) + " -> " + to_string(FinalTwoJ_)
           + " with L = " + std::to_string(L_)
           + " and S = " + spin_to_string(TwoS_);
}

//-------------------------
const std::set<int> SpinAmplitude::twoM() const
{
    std::set<int> S;
    // loop over amplitudes, key = two_M
    for (auto& kv : Amplitudes_)
        S.insert(kv.first);  // first entry is (twice) parent spin projection
    return S;
}

//-------------------------
const std::set<SpinProjectionVector> SpinAmplitude::twoM(int two_M) const
{
    auto it = Amplitudes_.find(two_M);

    // if two_M not found, return empty set
    if (it == Amplitudes_.end())
        return std::set<SpinProjectionVector>();

    std::set<SpinProjectionVector> two_m;
    std::transform(it->second.begin(), it->second.end(), std::inserter(two_m, two_m.end()), [](const auto& kv) {return kv.first;});
    return two_m;
}

//-------------------------
void SpinAmplitude::addAmplitude(int two_M, const SpinProjectionVector& two_m, bool store_null)
{
    // retrieve (or create) AmplitudeSubmap for two_M
    auto& ASM = Amplitudes_[two_M];

    // look for two_m in AmplitudeSubmap
    if (ASM.find(two_m) != ASM.end())
        throw exceptions::Exception("Amplitude already stored for " + spin_to_string(two_M) + " -> " + to_string(two_m),
                                    "SpinAmplitude::addAmplitude");

    if (store_null)
        ASM[two_m] = std::shared_ptr<ComplexCachedValue>(nullptr);
    else
        ASM[two_m] = ComplexCachedValue::create(*this);
}

//-------------------------
bool SpinAmplitude::equalTo(const SpinAmplitude& B) const
{
    return InitialTwoJ_ == B.InitialTwoJ_
        and FinalTwoJ_ == B.FinalTwoJ_
        and L_ == B.L_
        and TwoS_ == B.TwoS_;
}

//-------------------------
bool SpinAmplitude::equals(const SpinAmplitude& B) const
{
    return symmetrizationIndices() == B.symmetrizationIndices()
           and equalTo(B);
}

//-------------------------
std::string to_string(const SpinAmplitudeVector& saV)
{
    if (saV.empty())
        return std::string();

    std::string s = spin_to_string(saV[0]->initialTwoJ()) + " -> " + to_string(saV[0]->finalTwoJ());
    s += " with LS =";
    for (auto& sa : saV)
        s += " (" + std::to_string(sa->L()) + ", " + spin_to_string(sa->twoS()) + "),";
    s.erase(s.size() - 1, 1);
    if (saV[0]->formalism().empty())
        return s;
    s += " in " + saV[0]->formalism();
    return s;
}

}

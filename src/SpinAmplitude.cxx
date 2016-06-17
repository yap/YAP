#include "SpinAmplitude.h"

#include "CachedDataValue.h"
#include "CalculationStatus.h"
#include "Constants.h"
#include "Exceptions.h"
#include "logging.h"
#include "Spin.h"
#include "StatusManager.h"

namespace yap {

//-------------------------
SpinAmplitude::SpinAmplitude(unsigned two_J, const SpinVector& two_j, unsigned l, unsigned two_s, ParticleCombination::Equal& equal) :
    StaticDataAccessor(equal),
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

    // overall normalization constant
    double NJ = 1;              // no normalization
    // double NJ = sqrt((initialTwoJ() + 1) / 4. / pi<double>()); // sqrt((2J+1)/4pi)

    // loop over particle combinations -> indices
    for (auto& pc_i : symmetrizationIndices()) {

        // loop over mapping of parent spin projection to AmplitudeSubmap
        for (auto& aM_kv : Amplitudes_) {

            const auto& two_M = aM_kv.first; // parent spin projection

            // loop over mappin of daughter spin projection pairs to amplitudes
            for (auto& aSM_kv : aM_kv.second)

                // if yet uncalculated
                if (sm.status(*aSM_kv.second, pc_i.second) == CalculationStatus::uncalculated) {
                    aSM_kv.second->setValue(calc(two_M, aSM_kv.first, d, pc_i.first) * NJ,
                                            d, pc_i.second, sm);
                }
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
void SpinAmplitude::addAmplitude(int two_M, const SpinProjectionVector& two_m)
{
    // retrieve (or create) AmplitudeSubmap for two_M
    auto& ASM = Amplitudes_[two_M];

    // look for two_m in AmplitudeSubmap
    if (ASM.find(two_m) != ASM.end())
        throw exceptions::Exception("Amplitude already stored for " + spin_to_string(two_M) + " -> " + to_string(two_m),
                                    "SpinAmplitude::addAmplitude");

    FDEBUG("adding CachedDataValue for " << spin_to_string(two_M) << " -> " << to_string(two_m) << " in " << *this);
    ASM[two_m] = ComplexCachedDataValue::create(this);
}

//-------------------------
bool SpinAmplitude::equalTo(const SpinAmplitude& B) const
{
    return // compare only spin of QuantumNumbers
        InitialTwoJ_ == B.InitialTwoJ_
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

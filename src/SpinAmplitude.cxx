#include "SpinAmplitude.h"

#include "CachedDataValue.h"
#include "CalculationStatus.h"
#include "Constants.h"
#include "Exceptions.h"
#include "logging.h"
#include "spin.h"
#include "StatusManager.h"

namespace yap {

//-------------------------
SpinAmplitude::SpinAmplitude(unsigned two_J, unsigned two_j1, unsigned two_j2, unsigned l, unsigned two_s, ParticleCombination::Equal& equal) :
    StaticDataAccessor(equal),
    InitialTwoJ_(two_J),
    FinalTwoJ_( {two_j1, two_j2}),
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

                    const auto& spp = aSM_kv.first; // SpinProjectionPair of daughters

                    auto val = calc(two_M, spp[0], spp[1], d, pc_i.first) * NJ;

                    FDEBUG(*this << " := " << val << ", for " << *pc_i.first << " for " << two_M << " -> " << spp[0] << " + " << spp[1]);

                    aSM_kv.second->setValue(val, d, pc_i.second, sm);

                }
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

    FDEBUG("adding CachedDataValue for " << spin_to_string(two_M) << " -> " << m1m2[0] << " + " << m1m2[1]
           << " in " << *this)
    ASM[m1m2] = ComplexCachedDataValue::create(this);
}

//-------------------------
bool SpinAmplitude::equalTo(const SpinAmplitude& B) const
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
           and equalTo(B);
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

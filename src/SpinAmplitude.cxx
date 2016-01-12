#include "SpinAmplitude.h"

#include "ClebschGordan.h"
#include "Exceptions.h"
#include "logging.h"

namespace yap {

//-------------------------
SpinAmplitude::SpinAmplitude(const QuantumNumbers& initial,
                             const QuantumNumbers& final1,
                             const QuantumNumbers& final2,
                             unsigned L, unsigned two_S,
                             InitialStateParticle* isp) :
    StaticDataAccessor(isp),
    InitialQuantumNumbers_(initial),
    FinalQuantumNumbers_({final1, final2}),
    L_(L),
    TwoS_(two_S),
    Amplitude_(std::make_shared<ComplexCachedDataValue>(this))
{
    // check JLS triangle
    if (!triangle(InitialStateParticle_.twoJ(), 2 * l, two_S))
        throw exceptions::AngularMomentumNotConserved();
            
    // check j1j2S triangle
    if (!triangle(FinalQuantumNumbers_[0].twoJ(), FinalQuantumNumbers_[1].twoJ(), two_S))
        throw exceptions::AngularMomentumNotConserved();
      
    // if (!conserves(InitialQuantumNumbers_.twoJ(), FinalQuantumNumbers_[0].twoJ(), FinalQuantumNumbers_[1].twoJ(), l))
    //     throw exceptions::AngularMomentumNotConserved();
    
    // check charge conservation
    if (InitialQuantumNumbers_.Q() != FinalQuantumNumbers_[0].Q() + FinalQuantumNumbers_[1].Q())
        throw exceptions::Exception(std::string("charge conservation violated: ")
                                    + "(" + std::to_string(InitialQuantumNumbers_.Q())  + ") -> "
                                    + "(" + std::to_string(FinalQuantumNumbers_[0].Q()) + ") + "
                                    + "(" + std::to_string(FinalQuantumNumbers_[1].Q()) + ")",
                                    "SpinAmplitude::SpinAmplitude");
}

//-------------------------
SpinAmplitude::operator std::string() const
{
    std::string s = to_string(InitialQuantumNumbers_) + " -> ";
    for (auto& d : FinalQuantumNumbers_)
        s += to_string(d) + " + ";
    s.erase(s.size() - 2, 2);
    s += "with L = " + std::to_string(L_);
    s += " and S = " + spin_to_string(TwoS_);
    return s;
}

//-------------------------
std::set<int> twoM() const
{
    std::set<int> S;
    // loop over amplitudes, key = 3-array
    for (auto& kv : Amplitudes_)
        S.insert(kv.first[0]);  // first entry is (twice) parent spin projection
    return S;
}
        
//-------------------------
bool SpinAmplitude::equals(const SpinAmplitude& B) const
{
    return symmetrizationIndices() == B.symmetrizationIndices()
        and InitialQuantumNumbers_ == B.InitialQuantumNumbers_
        and FinalQuantumNumbers_ == B.FinalQuantumNumbers_
        and L_ == B.L_;
}

//-------------------------
CachedDataValueSet SpinAmplitude::amplitudeSet()
{
    CachedDataValueSet V;
    for (auto& kv : Amplitudes_)
        V.insert(kv.second);
    return V;
}

}

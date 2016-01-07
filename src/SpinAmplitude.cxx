#include "SpinAmplitude.h"

#include "ClebschGordan.h"
#include "Exceptions.h"
#include "logging.h"

namespace yap {

//-------------------------
SpinAmplitude::SpinAmplitude(const QuantumNumbers& initial, const QuantumNumbers& final1, const QuantumNumbers& final2, unsigned l) :
    InitialStateParticle_(nullptr),
    InitialQuantumNumbers_(initial),
    FinalQuantumNumbers_( {final1, final2}),
                      L_(l)
{
    if (!conserves(InitialQuantumNumbers_.twoJ(), FinalQuantumNumbers_[0].twoJ(), FinalQuantumNumbers_[1].twoJ(), l))
        throw exceptions::AngularMomentumNotConserved();
}

//-------------------------
bool SpinAmplitude::consistent() const
{
    bool C = true;

    // check charge conservation
    if (InitialQuantumNumbers_.Q() != FinalQuantumNumbers_[0].Q() + FinalQuantumNumbers_[1].Q()) {
        FLOG(ERROR) << "charge conservation violated: "
                    << "(" << InitialQuantumNumbers_.Q() << ") -> "
                    << "(" << FinalQuantumNumbers_[0].Q() << ") + "
                    << "(" << FinalQuantumNumbers_[1].Q() << ")";
        C &= false;
    }


    // check angular momentum conservation
    if (!conserves(InitialQuantumNumbers_.twoJ(), FinalQuantumNumbers_[0].twoJ(), FinalQuantumNumbers_[1].twoJ(), L_)) {
        FLOG(ERROR) << "angular momentum conservation violated: "
                    << "(" << spin_to_string(InitialQuantumNumbers_.twoJ()) << ") -> "
                    << "(" << spin_to_string(FinalQuantumNumbers_[0].twoJ()) << ") + "
                    << "(" << spin_to_string(FinalQuantumNumbers_[1].twoJ()) << ") "
                    << "with l = " << L_;
        C &= false;
    }

    return C;
}

//-------------------------
SpinAmplitude::operator std::string() const
{
    std::string s = to_string(InitialQuantumNumbers_) + " -> ";
    for (auto& d : FinalQuantumNumbers_)
        s += to_string(d) + " + ";
    s.erase(s.size() - 2, 2);
    s += "with l = " + std::to_string(L_);
    return s;
}

//-------------------------
bool SpinAmplitude::equals(const SpinAmplitude& B) const
{
    return symmetrizationIndices() == B.symmetrizationIndices()
           and InitialQuantumNumbers_ == B.InitialQuantumNumbers_
           and FinalQuantumNumbers_ == B.FinalQuantumNumbers_
           and L_ == B.L_;
}

}

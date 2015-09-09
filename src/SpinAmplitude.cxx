#include "SpinAmplitude.h"

#include "logging.h"

namespace yap {

//-------------------------
SpinAmplitude::SpinAmplitude(InitialStateParticle* isp, const QuantumNumbers& initial, const QuantumNumbers& final1, const QuantumNumbers& final2)
    : DataAccessor(isp),
      InitialQuantumNumbers_(initial),
      FinalQuantumNumbers_( {{final1, final2}})
{}

//-------------------------
bool SpinAmplitude::consistent() const
{
    bool consistent = true;

    // check charge conservation
    if (InitialQuantumNumbers_.Q() != FinalQuantumNumbers_[0].Q() + FinalQuantumNumbers_[1].Q()) {
        LOG(ERROR) << "SpinAmplitude::consistent() - charge conservation violated. " <<
                   "Q(parent) = " << (int)InitialQuantumNumbers_.Q() <<
                   "; Q(daughterA) = " << (int)FinalQuantumNumbers_[0].Q() <<
                   "; Q(daughterB) = " << (int)FinalQuantumNumbers_[1].Q();
        consistent =  false;
    }

    return consistent;
}

//-------------------------
bool SpinAmplitude::equals(const SpinAmplitude& rhs) const
{
    /*LOG(DEBUG) << "compare " << std::string(*this) << " and " << std::string(rhs);

    std::cout << "SymmetrizationIndices_ == rhs.SymmetrizationIndices_: " << (SymmetrizationIndices_ == rhs.SymmetrizationIndices_) << "\n";
    std::cout << "InitialQuantumNumbers_ == rhs.InitialQuantumNumbers_: " << (InitialQuantumNumbers_ == rhs.InitialQuantumNumbers_) << "\n";
    std::cout << "FinalQuantumNumbers_[0] == rhs.FinalQuantumNumbers_[0]: " << (FinalQuantumNumbers_[0] == rhs.FinalQuantumNumbers_[0]) << "\n";
    std::cout << "FinalQuantumNumbers_[1] == rhs.FinalQuantumNumbers_[1]: " << (FinalQuantumNumbers_[1] == rhs.FinalQuantumNumbers_[1]) << "\n";

    if (!(SymmetrizationIndices_ == rhs.SymmetrizationIndices_)) {
      for (auto& kv : SymmetrizationIndices_) {
        std::cout << std::string(*kv.first) << "=>" << kv.second << "\n";
      }
      std::cout << "--\n";
      for (auto& kv : rhs.SymmetrizationIndices_) {
        std::cout << std::string(*kv.first) << "=>" << kv.second << "\n";
      }
    }*/

    return (SymmetrizationIndices_ == rhs.SymmetrizationIndices_
            and InitialQuantumNumbers_ == rhs.InitialQuantumNumbers_
            and FinalQuantumNumbers_[0] == rhs.FinalQuantumNumbers_[0]
            and FinalQuantumNumbers_[1] == rhs.FinalQuantumNumbers_[1]);
}

}

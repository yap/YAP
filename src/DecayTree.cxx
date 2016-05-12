#include "DecayTree.h"

#include "Exceptions.h"

namespace yap {

//-------------------------
DecayTree::DecayTree(int two_M, const std::array<int, 2>& two_m, const std::shared_ptr<ComplexParameter>& free_amp) :
    TwoM_(two_M),
    DaughtersTwoM_(two_m),
    FreeAmplitude_(free_amp)
{}

//-------------------------
void DecayTree::setDaughterDecayTree(unsigned i, std::shared_ptr<DecayTree> dt)
{
    if (i >= DaughtersTwoM_.size())
        throw exceptions::Exception("index exceeds number of daughters", "DecayTree::setDaughterDecayTree");

    if (dt->TwoM_ != DaughtersTwoM_[i])
        throw exceptions::Exception("Spin projection mismatch", "DecayTree::setDaughterDecayTree");

    if (DaughterDecayTrees_.find(i) != DaughterDecayTrees_.end())
        throw exceptions::Exception("DecayTree for this daughter already set", "DecayTree::setDaughterDecayTree");

    DaughterDecayTrees_[i] = dt;
}

//-------------------------
// std::complex<double> operator()(DataPoint& d) const
// {
//     auto A = Complex_1;

//     // multiply by all free amplitudes
//     for (const auto& a : FreeAmplitudes_)
//         A *= a->value();

//     // multiply by all fixed amplitudes
//     for (const auto& a : FixedAmplitudes_)
//         A *= ;

//     return A;
// }

}

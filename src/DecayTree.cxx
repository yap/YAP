#include "DecayTree.h"

#include "Exceptions.h"
#include "FreeAmplitude.h"

namespace yap {

//-------------------------
DecayTree::DecayTree(std::shared_ptr<FreeAmplitude> free_amp) :
    FreeAmplitude_(free_amp)
{}

//-------------------------
void DecayTree::setDaughterDecayTree(unsigned i, std::shared_ptr<DecayTree> dt)
{
    if (i >= DaughtersTwoM_.size())
        throw exceptions::Exception("index exceeds number of daughters", "DecayTree::setDaughterDecayTree");

    if (!dt->freeAmplitude())
        throw exceptions::Exception("DecayTree's free amplitude is nullptr", "DecayTree:setDaughterDecayTree");

    if (DaughterDecayTrees_.find(i) != DaughterDecayTrees_.end())
        throw exceptions::Exception("DecayTree for this daughter already set", "DecayTree::setDaughterDecayTree");

    DaughterDecayTrees_[i] = dt;
    setDaughterSpinProjection(i, dt->freeAmplitude()->twoM());
}

//-------------------------
std::string DecayTree::asString(std::string offset) const
{
    auto s = to_string(*FreeAmplitude_);
    offset += "    ";
    for (const auto& d_dt : DaughterDecayTrees_)
        s += "\n" + offset + "d[" + std::to_string(d_dt.first) + "] --> " + d_dt.second->asString(offset);
    return s;
}

//-------------------------
unsigned depth(const DecayTree& DT)
{
    unsigned d = 0;
    for (const auto& i_dt : DT.daughterDecayTrees())
        d = std::max(d, depth(*i_dt.second));
    return d + 1;
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

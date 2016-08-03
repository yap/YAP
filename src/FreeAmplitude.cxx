#include "FreeAmplitude.h"

#include "DecayChannel.h"
#include "container_utils.h"
#include "Exceptions.h"
#include "ParticleCombination.h"
#include "Spin.h"
#include "SpinAmplitude.h"

namespace yap {

//-------------------------
FreeAmplitude::FreeAmplitude(std::shared_ptr<DecayChannel> dc, std::shared_ptr<SpinAmplitude> sa,
                             int two_m, std::complex<double> a) :
    ComplexParameter(a),
    DecayChannel_(dc),
    SpinAmplitude_(sa),
    TwoM_(two_m)
{
    if (!DecayChannel_)
        throw exceptions::Exception("DecayChannel is nullptr", "FreeAmplitude::FreeAmplitude");

    if (!SpinAmplitude_)
        throw exceptions::Exception("SpinAmpliutde is nullptr", "FreeAmplitude::FreeAmplitude");

    // if (!checkParticleCombinations(*SpinAmplitude_))
    //     throw exceptions::Exception("SpinAmplitude does not contain all of DecayChannel's ParticleCombinations",
    //                                 "FreeAmplitude::FreeAmplitude");

    // check M
    auto s_M = SpinAmplitude_->twoM();
    if (s_M.find(TwoM_) == s_M.end())
        throw exceptions::Exception("SpinAmplitude not valid for spin projection " + spin_to_string(TwoM_),
                                    "FreeAmplitude::FreeAmplitude");
}

//-------------------------
const Model* FreeAmplitude::model() const
{
    return (DecayChannel_) ? DecayChannel_->model() : nullptr;
}

//-------------------------
bool FreeAmplitude::checkParticleCombinations(const DataAccessor& da) const
{
    if (!DecayChannel_)
        throw exceptions::Exception("DecayChannel_ is nullptr", "FreeAmplitude::checkDataAccessor");

    // check DataAccessor contains all particle combination of DecayChannel
    return contains(da.symmetrizationIndices().begin(), da.symmetrizationIndices().end(),
                    DecayChannel_->particleCombinations().begin(), DecayChannel_->particleCombinations().end(),
                    [&](const ParticleCombinationMap<unsigned>::value_type & pc_i, const ParticleCombinationVector::value_type & pc)
                    {return pc_i.first == pc;});
}

//-------------------------
std::string to_string(const FreeAmplitude& fa)
{
    return to_string(*fa.decayChannel())
           + ", M = " + spin_to_string(fa.twoM())
           + ", " + to_string(*fa.spinAmplitude())
           + (fa.variableStatus() == VariableStatus::fixed ? " [fixed]" : "");
}

}

#include "FreeAmplitude.h"

#include "DecayChannel.h"
#include "Exceptions.h"
#include "spin.h"
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

    if (!checkParticleCombinations(*SpinAmplitude_))
        throw exceptions::Exception("SpinAmplitude does not contain all of DecayChannel's ParticleCombinations",
                                    "FreeAmplitude::FreeAmplitude");

    // check M
    auto s_M = SpinAmplitude_->twoM();
    if (s_M.find(TwoM_) == s_M.end())
        throw exceptions::Exception("SpinAmplitude not valid for spin projection " + spin_to_string(TwoM_),
                                    "FreeAmplitude::FreeAmplitude");
}

//-------------------------
bool FreeAmplitude::checkParticleCombinations(const ReportsParticleCombinations& PC) const
{
    if (!DecayChannel_)
        throw exceptions::Exception("DecayChannel_ is nullptr", "FreeAmplitude::checkDataAccessor");

    return std::all_of(DecayChannel_->particleCombinations().begin(),
                       DecayChannel_->particleCombinations().end(),
                       [&PC](const ParticleCombinationVector::value_type & pc)
    {return PC.hasParticleCombination(pc);});
}

//-------------------------
std::string to_string(const FreeAmplitude& fa)
{
    return to_string(*fa.decayChannel())
           + ", M = " + spin_to_string(fa.twoM())
           + ", " + to_string(*fa.spinAmplitude())
           + (fa.variableStatus() == VariableStatus::fixed ? " [fixed]" : "");
}


//-------------------------
FreeAmplitudeSet find(const FreeAmplitudeSet& fas, int two_m)
{
    FreeAmplitudeSet S;
    std::copy_if(fas.begin(), fas.end(), std::inserter(S, S.end()), [&two_m](const FreeAmplitudeSet::value_type & fa) {return fa->twoM() == two_m;});
    return S;
}

//-------------------------
FreeAmplitudeSet find(const FreeAmplitudeSet& fas, const DecayChannel* const dc)
{
    FreeAmplitudeSet S;
    std::copy_if(fas.begin(), fas.end(), std::inserter(S, S.end()), [&dc](const FreeAmplitudeSet::value_type & fa) {return fa->decayChannel().get() == dc;});
    return S;
}

//-------------------------
FreeAmplitudeSet find(const FreeAmplitudeSet& fas, const SpinAmplitude* const sa)
{
    FreeAmplitudeSet S;
    std::copy_if(fas.begin(), fas.end(), std::inserter(S, S.end()), [&sa](const FreeAmplitudeSet::value_type & fa) {return fa->spinAmplitude().get() == sa;});
    return S;
}

}

#include "AmplitudeComponentDataAccessor.h"

namespace yap {

//-------------------------
AmplitudeComponentDataAccessor::AmplitudeComponentDataAccessor(InitialStateParticle* isp, ParticleCombination::Equiv* equiv) :
    AmplitudeComponent(),
    DataAccessor(isp, equiv)
{
}


//-------------------------
const Amp& AmplitudeComponentDataAccessor::amplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc)
{
    // find symmetrization index
    unsigned sym_index = SymmetrizationIndices_[pc];

    Amp& a = cachedAmplitude(d, sym_index);

    // check whether data-dependent calculation needs to be made
    CalculationStatus& calcStat = CalculationStatuses(d, sym_index);
    if (calcStat == kUncalculated) {
        // calculate amplitude
        a = calcAmplitude(d, pc);

        // set calculation status
        calcStat = kCalculated;
    }

    return a;
}

}


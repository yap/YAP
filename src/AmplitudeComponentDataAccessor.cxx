#include "AmplitudeComponentDataAccessor.h"

#include "InitialStateParticle.h"

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
    CalculationStatus& calcStat = CalculationStatuses_[sym_index];
    if (calcStat == kUncalculated) {

        // calculate amplitude for ALL dataPoints
        for (DataPoint& dataPt : initialStateParticle()->dataSet()) {
            cachedAmplitude(dataPt, sym_index) = calcAmplitude(dataPt, pc);
        }

        // set calculation status
        calcStat = kCalculated;
    }

    return a;
}

}


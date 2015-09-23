#include "AmplitudeComponentDataAccessor.h"

#include "Constants.h"
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
    // has to made shure in caller!
    //if (! hasSymmetrizationIndex(pc))
    //    return Complex_0;

    // find symmetrization index
    unsigned sym_index = SymmetrizationIndices_.at(pc);

    Amp& a = cachedAmplitude(d, sym_index);

    // check whether data-dependent calculation needs to be made
    CalculationStatus& calcStat = CalculationStatuses_.at(sym_index);
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


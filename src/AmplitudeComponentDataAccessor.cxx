#include "AmplitudeComponentDataAccessor.h"

#include "Constants.h"
#include "InitialStateParticle.h"
#include "logging.h"

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
    // has to made sure in caller!
    /*if (! hasSymmetrizationIndex(pc)) {
        LOG(ERROR) << "AmplitudeComponentDataAccessor::amplitude - called with wrong symmetrization index!";
        static Amp zero(Complex_0);
        return zero;
    }*/

    // find symmetrization index
    unsigned sym_index = SymmetrizationIndices_.at(pc);

    Amp& a = cachedAmplitude(d, sym_index);

    // check whether data-dependent calculation needs to be made
    if (calculationStatus(sym_index) == kUncalculated) {

        // calculate amplitude for ALL dataPoints
        for (DataPoint& dataPt : initialStateParticle()->dataSet()) {
            cachedAmplitude(dataPt, sym_index) = calcAmplitude(dataPt, pc);
        }

        // set calculation status
        setCalculationStatus(sym_index, kCalculated);
    }

    return a;
}

}


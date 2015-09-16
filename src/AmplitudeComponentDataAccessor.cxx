#include "AmplitudeComponentDataAccessor.h"

namespace yap {

//-------------------------
AmplitudeComponentDataAccessor::AmplitudeComponentDataAccessor(InitialStateParticle* isp, ParticleCombination::Equiv* equiv) :
    AmplitudeComponent(),
    DataAccessor(isp, equiv)
{
}


//-------------------------
Amp AmplitudeComponentDataAccessor::amplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc)
{
    // find symmetrization index
    unsigned sym_index = SymmetrizationIndices_[pc];

    // check whether data-dependent calculation needs to be made
    CalculationStatus& calcStat = CalculationStatuses(d, sym_index);
    if (calcStat == kUncalculated) {

        // calculate amplitude
        Amp a = calcAmplitude(d, pc);

        // store into data point
        data(d, sym_index) = {real(a), imag(a)};

        // set calculation status
        calcStat = kCalculated;

        return a;
    }

    // else return from cached data
    const std::vector<double>& D = data(d, sym_index);
    return Amp(D[0], D[1]);
}

}


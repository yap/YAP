#include "AmplitudeComponentDataAccessor.h"

#include "Constants.h"
#include "InitialStateParticle.h"
#include "logging.h"

namespace yap {

//-------------------------
AmplitudeComponentDataAccessor::AmplitudeComponentDataAccessor(InitialStateParticle* isp, ParticleCombination::Equiv* equiv) :
    DataAccessor(isp, equiv)
{
}


//-------------------------
const Amp& AmplitudeComponentDataAccessor::amplitude(DataPartition& d, std::shared_ptr<const ParticleCombination> pc) const
{
#ifdef ELPP_DISABLE_DEBUG_LOGS
    // find symmetrization index
    unsigned sym_index = SymmetrizationIndices_[pc];
#else
    // has to made sure in caller! This is just for safety while debugging!
    if (! hasSymmetrizationIndex(pc)) {
        LOG(ERROR) << "AmplitudeComponentDataAccessor::amplitude - called with wrong symmetrization index!";
        static Amp zero(Complex_0);
        return zero;
    }

    // find symmetrization index
    unsigned sym_index = SymmetrizationIndices_.at(pc);
#endif

    Amp& a = cachedAmplitude(d.dataPoint(), sym_index);

    CalculationStatus& calcStat = d.CalculationStatusesDataPoint(index(), sym_index);

    if (calcStat == kCalculated and d.CalculationStatusesDataSet(index(), sym_index) == kCalculated)
        return a;

    // calculate amplitude and store in cache
    a = calcAmplitude(d, pc);

    // set calculation status for dataPoint
    calcStat = kCalculated;

    return a;
}

}


#include "BreitWigner.h"

#include "Constants.h"
#include "FourMomenta.h"
#include "logging.h"
#include "InitialStateParticle.h"
#include "ParticleCombination.h"

namespace yap {

//-------------------------
BreitWigner::BreitWigner(InitialStateParticle* isp, double mass, double width) :
    MassShape(isp)
{
    ParameterSet::operator=({mass, width});
}

//-------------------------
Amp BreitWigner::amplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc)
{
    // find symmetrization index
    unsigned sym_index = SymmetrizationIndices_[pc];

    // check whether data-dependent calculation needs to be made
    CalculationStatus& calcStat = CalculationStatuses(d, sym_index);
    if (calcStat == kUncalculated) {

        // check whether data-independent calculation needs to be made
        if (CalcStatus_ == kUncalculated) {
            M2iMG_ = Amp(mass() * mass(), -mass() * width());
            CalcStatus_ = kCalculated;
        }

        // calculate amplitude
        Amp a = 1. / (M2iMG_ - Amp(initialStateParticle()->fourMomenta().m2(d, pc), 0));

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

//-------------------------
Amp BreitWigner::amplitude(double s)
{
    if (CalcStatus_ == kUncalculated)
        M2iMG_ = Amp(mass() * mass(), -mass() * width());
    return 1. / (M2iMG_ - Amp(s, 0));
}

//-------------------------
bool BreitWigner::consistent() const
{
    bool consistent = true;
    if (mass() <= 0) {
        LOG(ERROR) << "BreitWigner::consistent() - mass <= 0";
        consistent = false;
    }
    if (width() <= 0) {
        LOG(ERROR) << "BreitWigner::consistent() - width <= 0";
        consistent = false;
    }

    return consistent;
}

}





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
Amp BreitWigner::calcAmplitudeS(double s)
{
    if (CalcStatus_ == kUncalculated)
        M2iMG_ = Amp(mass() * mass(), -mass() * width());

    LOG(DEBUG) << "BreitWigner amplitude (s = " << s << ") = " << 1. / (M2iMG_ - Amp(s, 0));

    return 1. / (M2iMG_ - Amp(s, 0));
}

//-------------------------
bool BreitWigner::consistent() const
{
    bool consistent = MassShape::consistent();

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





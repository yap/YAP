#include "BreitWigner.h"

#include "Constants.h"
#include "FourMomenta.h"
#include "logging.h"
#include "InitialStateParticle.h"
#include "ParticleCombination.h"

namespace yap {

//-------------------------
BreitWigner::BreitWigner(double mass, double width) :
    MassShape( {mass, width})
{
}

//-------------------------
void BreitWigner::calcPrecalculate()
{
    M2iMG_ = std::complex<double>(mass() * mass(), -mass() * width());
}

//-------------------------
std::complex<double> BreitWigner::calcAmplitudeS(double s) const
{
    std::complex<double> a = 1. / (M2iMG_ - Complex_1 * s);

    DEBUG("BreitWigner amplitude (s = " << s << ") = " << a);

    return a;
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





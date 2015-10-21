#include "BreitWigner.h"

#include "Constants.h"
#include "FourMomenta.h"
#include "logging.h"
#include "InitialStateParticle.h"
#include "ParticleCombination.h"

namespace yap {

//-------------------------
BreitWigner::BreitWigner(double mass, double width) :
    MassShape( {mass, width}),
           M2iMG_(new ComplexCachedValue({Parameters_[0], Parameters_[1]}))
{
}

std::complex<double> BreitWigner::amplitude(DataPartition& d, std::shared_ptr<const ParticleCombination> pc) const
{
    /// \todo implement
    return Complex_1;
}

//-------------------------
/*void BreitWigner::precalculate()
{
    if (M2iMG_->calculationStatus() == kUncalculated) {
        // m*m -i*m*w
        M2iMG_->setValue(Parameters_[0]->value() * Parameters_[0]->value() - Complex_i * Parameters_[0]->value() * Parameters_[1]->value());
    }
}*/

//-------------------------
std::complex<double> BreitWigner::calcAmplitudeS(double s) const
{
    std::complex<double> a = 1. / (M2iMG_->value() - Complex_1 * s);

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





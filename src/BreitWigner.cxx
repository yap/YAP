#include "BreitWigner.h"

#include "Constants.h"
#include "FourMomenta.h"
#include "logging.h"
#include "InitialStateParticle.h"
#include "ParticleCombination.h"
#include "Resonance.h"

namespace yap {

//-------------------------
BreitWigner::BreitWigner(double mass, double width) :
    MassShape( {std::make_shared<RealParameter>(mass), std::make_shared<RealParameter>(width)}),
           M2iMG_(new ComplexCachedValue({Parameters_[0], Parameters_[1]}))
{
}

//-------------------------
void BreitWigner::borrowParametersFromResonance(Resonance* R)
{
    // Remove existing mass parameter from M2iMG_
    M2iMG_->removeDependency(mass());
    // borrow mass from Owner_
    Parameters_[0] = R->mass();
    // add new mass parameter into M2iMG_
    M2iMG_->addDependency(mass());
}

//-------------------------
std::complex<double> BreitWigner::amplitude(DataPartition& d, const std::shared_ptr<const ParticleCombination>& pc) const
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

    if (mass()->value() <= 0) {
        LOG(ERROR) << "BreitWigner::consistent() - mass <= 0";
        consistent = false;
    }
    if (width()->value() <= 0) {
        LOG(ERROR) << "BreitWigner::consistent() - width <= 0";
        consistent = false;
    }

    return consistent;
}

}





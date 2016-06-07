#include "PoleMass.h"

#include "CalculationStatus.h"
#include "DataPartition.h"
#include "FourMomenta.h"
#include "logging.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleFactory.h"
#include "Resonance.h"
#include "StatusManager.h"

namespace yap {

//-------------------------
PoleMass::PoleMass(std::complex<double> mass) :
    MassShape(),
    Mass_(std::make_shared<ComplexParameter>(mass))
{
}

//-------------------------
void PoleMass::setParameters(const ParticleTableEntry& entry)
{
    auto m = mass()->value();
    if (real(m) < 0)
        m.real(entry.Mass);
    if (imag(m) < 0 and !entry.MassShapeParameters.empty())
        m.imag(entry.MassShapeParameters[0] / 2.);
    Mass_->setValue(m);
}

//-------------------------
void PoleMass::setResonance(Resonance* r)
{
    MassShape::setResonance(r);

    // set real component of mass from resonance, if yet unset
    if (real(Mass_->value()) < 0)
        Mass_->setValue(std::complex<double>(resonance()->mass()->value(), imag(Mass_->value())));

    // replace resonance's mass with real component of pole mass
    replaceResonanceMass(std::make_shared<RealComponentParameter>(Mass_));
}

//-------------------------
std::complex<double> PoleMass::amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc, StatusManager& sm) const
{
    unsigned symIndex = symmetrizationIndex(pc);

    // recalculate, cache, & return, if necessary
    if (sm.status(*T(), symIndex) == CalculationStatus::uncalculated) {

        // T = 1 / (M^2 - m^2)
        std::complex<double> t = 1. / (pow(Mass_->value(), 2) - model()->fourMomenta()->m2(d, pc));

        T()->setValue(t, d, symIndex, sm);

        FDEBUG("calculated T = " << t << " and stored it in the cache");
        return t;
    }

    FDEBUG("using cached T = " << T()->value(d, symIndex));

    // else return cached value
    return T()->value(d, symIndex);
}

//-------------------------
void PoleMass::calculateT(DataPartition& D, const std::shared_ptr<ParticleCombination>& pc, unsigned si) const
{
    auto M2 = pow(Mass_->value(), 2);

    // T := 1 / (M^2 - m^2)
    for (auto& d : D)
        T()->setValue(1. / (M2 - model()->fourMomenta()->m2(d, pc)), d, si, D);
}

//-------------------------
bool PoleMass::consistent() const
{
    bool C = MassShape::consistent();

    if (real(Mass_->value()) <= 0) {
        FLOG(ERROR) << "real(mass) <= 0";
        C &= false;
    }
    if (imag(Mass_->value()) <= 0) {
        FLOG(ERROR) << "imag(mass) <= 0";
        C &= false;
    }

    return C;
}

}





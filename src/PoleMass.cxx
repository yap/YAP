#include "PoleMass.h"

#include "FourMomenta.h"
#include "logging.h"
#include "Model.h"
#include "Resonance.h"

namespace yap {

//-------------------------
PoleMass::PoleMass(std::complex<double> mass) :
    MassShape(),
    Mass_(std::make_shared<ComplexParameter>(mass)),
    T_(ComplexCachedDataValue::create(this, ParameterSet{Mass_}))
{
}

//-------------------------
void PoleMass::setParameters(const ParticleTableEntry& entry)
{
    auto mass = Complex_1 * entry.Mass;
    if (!entry.MassShapeParameters.empty())
        mass += Complex_i * entry.MassShapeParameters[0] / 2.;
    Mass_->setValue(mass);
}

//-------------------------
void PoleMass::borrowParametersFromResonance()
{
    // set real component of mass from resonance
    Mass_->setValue(std::complex<double>(resonance()->mass()->value(), imag(Mass_->value())));

    // replace resonance's mass with real component of pole mass
    replaceResonanceMass(std::make_shared<RealComponentParameter>(Mass_));
}

//-------------------------
void PoleMass::setDependenciesFromModel()
{
    if (!model())
        throw exceptions::Exception("Model unset", "PoleMass::setDependenciesFromResonance");
    if (!model()->fourMomenta())
        throw exceptions::Exception("Model's FourMomenta unset", "PoleMass::setDependenciesFromResonance");

    T_->addDependency(model()->fourMomenta()->mass());
}

//-------------------------
std::complex<double> PoleMass::amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc, unsigned dataPartitionIndex) const
{
    unsigned symIndex = symmetrizationIndex(pc);

    // recalculate, cache, & return, if necessary
    if (T_->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated) {

        // T = 1 / (M^2 - m^2)
        std::complex<double> T = 1. / (pow(Mass_->value(), 2) - model()->fourMomenta()->m2(d, pc));

        T_->setValue(T, d, symIndex, dataPartitionIndex);

        FDEBUG("calculated T = " << T << " and stored it in the cache");
        return T;
    }

    FDEBUG("using cached T = " << T_->value(d, symIndex));

    // else return cached value
    return T_->value(d, symIndex);
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





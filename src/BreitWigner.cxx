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
    MassShape(),
    Mass_(new RealParameter(mass)),
    Width_(new RealParameter(width)),
    T_(new ComplexCachedDataValue(this, {Mass_, Width_}))
{
}

//-------------------------
bool BreitWigner::setParameters(const ParticleTableEntry& entry)
{
    Mass_->setValue(entry.Mass_);

    if (entry.MassShapeParameters_.size() < 1)
        return false;

    Width_->setValue(entry.MassShapeParameters_[0]);
    return true;
}

//-------------------------
void BreitWigner::borrowParametersFromResonance(Resonance* R)
{
    // Remove existing mass parameter from M2iMG_
    T_->removeDependency(Mass_);
    // borrow mass from Owner_
    Mass_ = R->mass();
    // add new mass parameter into M2iMG_
    T_->addDependency(Mass_);
}

//-------------------------
std::complex<double> BreitWigner::amplitude(DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc, unsigned dataPartitionIndex) const
{
    unsigned symIndex = symmetrizationIndex(pc);

    // recalculate, cache, & return, if necessary
    if (T_->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated) {

        // T = 1 / (M^2 - m^2 - iMG)
        std::complex<double> T = 1. / (pow(Mass_->value(), 2) - initialStateParticle()->fourMomenta().m2(d, pc) - Complex_i * Mass_->value() * Width_->value());

        T_->setValue(T, d, symIndex, dataPartitionIndex);

        DEBUG("BreitWigner::amplitude - calculated T = " << T << " and stored it in the cache");
        return T;
    }

    DEBUG("BreitWigner::amplitude - using cached T = " << T_->value(d, symIndex));

    // else return cached value
    return T_->value(d, symIndex);
}

//-------------------------
bool BreitWigner::consistent() const
{
    bool consistent = MassShape::consistent();

    if (Mass_->value() <= 0) {
        LOG(ERROR) << "BreitWigner::consistent() - mass <= 0";
        consistent = false;
    }
    if (Width_->value() <= 0) {
        LOG(ERROR) << "BreitWigner::consistent() - width <= 0";
        consistent = false;
    }

    return consistent;
}

}





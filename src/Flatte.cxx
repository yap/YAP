#include "Flatte.h"

#include "Exceptions.h"
#include "FourMomenta.h"
#include "logging.h"
#include "Model.h"
#include "Resonance.h"

namespace yap {

//-------------------------
Flatte::Flatte(double mass) :
    MassShape(),
    Mass_(std::make_shared<RealParameter>(mass)),
    WidthTerm_(ComplexCachedDataValue::create(this)),
    T_(ComplexCachedDataValue::create(this, ParameterSet{Mass_}, CachedDataValueSet{WidthTerm_}))
{
}

//-------------------------
void Flatte::addChannel(std::shared_ptr<RealParameter> coupling, std::shared_ptr<RealParameter> mass)
{
    if (!coupling)
        throw exceptions::Exception("Coupling is unset", "Flatte::addChannel");
    if (!mass)
        throw exceptions::Exception("Mass is unset", "Flatte::addChannel");
    FlatteChannels_.push_back(FlatteChannel(coupling, mass));
    WidthTerm_->addDependency(FlatteChannels_.back().Coupling);
    WidthTerm_->addDependency(FlatteChannels_.back().Mass);
}

//-------------------------
void Flatte::addChannel(double coupling, double mass)
{
    addChannel(std::make_shared<RealParameter>(coupling), std::make_shared<RealParameter>(mass));
}

//-------------------------
void Flatte::setParameters(const ParticleTableEntry& entry)
{
    Mass_->setValue(entry.Mass);
}

//-------------------------
void Flatte::setDependenciesFromModel()
{
    WidthTerm_->addDependency(model()->fourMomenta()->mass());
    T_->addDependency(model()->fourMomenta()->mass());
}

//-------------------------
void Flatte::borrowParametersFromResonance()
{
    // Remove existing mass parameter from M2iMG_
    T_->removeDependency(Mass_);

    // borrow mass from Owner_
    Mass_ = resonance()->mass();

    // add new mass parameter into M2iMG_
    T_->addDependency(Mass_);
}

//-------------------------
std::complex<double> Flatte::amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc, unsigned dataPartitionIndex) const
{
    unsigned symIndex = symmetrizationIndex(pc);

    if (WidthTerm_->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated) {
        auto w = Complex_0;
        // sum of coupling * complex-breakup-momentum
        for (const auto& fc : FlatteChannels_)
            w += fc.Coupling->value() * std::sqrt(std::complex<double>(model()->fourMomenta()->m2(d, pc) / 4. - pow(fc.Mass->value(), 2), 0));
        // sum * i * 2 / mass
        w *= Complex_i * 2. / model()->fourMomenta()->m(d, pc);
        WidthTerm_->setValue(w, d, symIndex, dataPartitionIndex);
    }

    // recalculate, cache, & return, if necessary
    if (T_->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated) {

        // T = 1 / (M^2 - m^2 - width-term)
        std::complex<double> T = 1. / (pow(Mass_->value(), 2) - model()->fourMomenta()->m2(d, pc) - WidthTerm_->value(d, symIndex));

        T_->setValue(T, d, symIndex, dataPartitionIndex);

        FDEBUG("calculated T = " << T << " and stored it in the cache");
        return T;
    }

    FDEBUG("using cached T = " << T_->value(d, symIndex));

    // else return cached value
    return T_->value(d, symIndex);
}

//-------------------------
bool Flatte::consistent() const
{
    bool C = MassShape::consistent();

    if (Mass_->value() <= 0) {
        FLOG(ERROR) << "mass <= 0";
        C &= false;
    }
    for (const auto& fc : FlatteChannels_) {
        if (fc.Coupling->value() <= 0) {
            FLOG(ERROR) << "coupling constant <= 0";
            C &= false;
        }
        if (fc.Mass->value() <= 0) {
            FLOG(ERROR) << "mass <= 0";
            C &= false;
        }
    }

    return C;
}

}





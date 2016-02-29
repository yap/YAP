#include "Flatte.h"

#include "Exceptions.h"
#include "FourMomenta.h"
#include "logging.h"
#include "Model.h"
#include "Resonance.h"

namespace yap {

//-------------------------
Flatte::Flatte() :
    MassShapeWithNominalMass(),
    WidthTerm_(ComplexCachedDataValue::create(this))
{
    T()->addDependency(WidthTerm_);
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
void Flatte::setDependenciesFromModel()
{
    MassShapeWithNominalMass::setDependenciesFromModel();
    WidthTerm_->addDependency(model()->fourMomenta()->mass());
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
    if (T()->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated) {

        // T = 1 / (M^2 - m^2 - width-term)
        std::complex<double> t = 1. / (pow(mass()->value(), 2) - model()->fourMomenta()->m2(d, pc) - WidthTerm_->value(d, symIndex));

        T()->setValue(t, d, symIndex, dataPartitionIndex);

        FDEBUG("calculated T = " << t << " and stored it in the cache");
        return t;
    }

    FDEBUG("using cached T = " << T()->value(d, symIndex));

    // else return cached value
    return T()->value(d, symIndex);
}

//-------------------------
bool Flatte::consistent() const
{
    bool C = MassShapeWithNominalMass::consistent();

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





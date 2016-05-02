#include "Flatte.h"

#include "CalculationStatus.h"
#include "Constants.h"
#include "DataPoint.h"
#include "DataPartition.h"
#include "Exceptions.h"
#include "FourMomenta.h"
#include "logging.h"
#include "Model.h"
#include "Parameter.h"
#include "Resonance.h"
#include "StatusManager.h"

namespace yap {

//-------------------------
void Flatte::addChannel(std::shared_ptr<RealParameter> coupling, std::shared_ptr<RealParameter> mass)
{
    if (!coupling)
        throw exceptions::Exception("Coupling is unset", "Flatte::addChannel");
    if (!mass)
        throw exceptions::Exception("Mass is unset", "Flatte::addChannel");
    FlatteChannels_.push_back(FlatteChannel(coupling, mass));
    T()->addDependency(FlatteChannels_.back().Coupling);
    T()->addDependency(FlatteChannels_.back().Mass);
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
    T()->addDependency(model()->fourMomenta()->mass());
}

//-------------------------
std::complex<double> Flatte::amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc, StatusManager& sm) const
{
    unsigned symIndex = symmetrizationIndex(pc);

    if (sm.status(*T(), symIndex) == CalculationStatus::uncalculated) {

        // calculate width term
        auto w = Complex_0;
        // sum of coupling * complex-breakup-momentum
        for (const auto& fc : FlatteChannels_)
            w += fc.Coupling->value() * std::sqrt(std::complex<double>(model()->fourMomenta()->m2(d, pc) / 4. - pow(fc.Mass->value(), 2), 0));
        // sum * i * 2 / mass
        w *= Complex_i * 2. / model()->fourMomenta()->m(d, pc);

        // T = 1 / (M^2 - m^2 - width-term)
        std::complex<double> t = 1. / (pow(mass()->value(), 2) - model()->fourMomenta()->m2(d, pc) - w);

        T()->setValue(t, d, symIndex, sm);

        FDEBUG("calculated T = " << t << " and stored it in the cache");
        return t;
    }

    FDEBUG("using cached T = " << T()->value(d, symIndex));

    // else return cached value
    return T()->value(d, symIndex);
}

//-------------------------
void Flatte::calculate(DataPartition& D, const std::shared_ptr<ParticleCombination>& pc) const
{
    unsigned symIndex = symmetrizationIndex(pc);

    if (D.status(*T(), symIndex) == CalculationStatus::uncalculated) {

        for (auto& d : D) {

            // calculate width term
            auto w = Complex_0;
            // sum of coupling * complex-breakup-momentum
            for (const auto& fc : FlatteChannels_)
                w += fc.Coupling->value() * std::sqrt(std::complex<double>(model()->fourMomenta()->m2(d, pc) / 4. - pow(fc.Mass->value(), 2), 0));
            // sum * i * 2 / mass
            w *= Complex_i * 2. / model()->fourMomenta()->m(d, pc);

            // T = 1 / (M^2 - m^2 - width-term)
            std::complex<double> t = 1. / (pow(mass()->value(), 2) - model()->fourMomenta()->m2(d, pc) - w);

            T()->setValue(t, d, symIndex, D);
        }

        D.status(*T(), symIndex) = CalculationStatus::calculated;
    }
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





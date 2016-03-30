#include "HelicityAngles.h"

#include "CalculationStatus.h"
#include "CoordinateSystem.h"
#include "FourMomenta.h"
#include "FourVector.h"
#include "LorentzTransformation.h"
#include "Model.h"
#include "ParticleCombination.h"
#include "ThreeVector.h"

namespace yap {

//-------------------------
HelicityAngles::HelicityAngles(Model* m) :
    StaticDataAccessor(m, &ParticleCombination::equivUpAndDown),
    Phi_(RealCachedDataValue::create(this)),
    Theta_(RealCachedDataValue::create(this))
{
    /// \todo add check that FourMomenta exists, after changing to return shared_ptr
    Phi_->addDependency(model()->fourMomenta()->momentum());
    Theta_->addDependency(model()->fourMomenta()->momentum());
}

//-------------------------
void HelicityAngles::calculate(DataPoint& d, StatusManager& sm) const
{
    // set angles uncalculated
    sm.set(*this, CalculationStatus::uncalculated);

    // call on ISP PC's
    // \todo allow for designating the boost that takes from the data frame to the lab frame (possibly event dependent)
    for (auto& kv : symmetrizationIndices())
        if (kv.first->indices().size() == model()->finalStateParticles().size())
            calculateAngles(d, kv.first, model()->coordinateSystem(), unitMatrix<double, 4>(), sm);
}

//-------------------------
void HelicityAngles::calculateAngles(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc,
                                     const CoordinateSystem<double, 3>& C, const FourMatrix<double>& boosts,
                                     StatusManager& sm) const
{
    // terminate recursion
    if (pc->isFinalStateParticle())
        return;

    // get pc's 4-mom in data frame
    const FourVector<double> P = model()->fourMomenta()->p(d, pc);

    // calculate reference frame for P from parent's RF
    const CoordinateSystem<double, 3> cP = helicityFrame(boosts * P, C);

    // calculate boost from data frame into pc rest frame
    const FourMatrix<double> boost = lorentzTransformation(-(boosts * P));

    const unsigned symIndex = symmetrizationIndex(pc);

    for (auto& daughter : pc->daughters()) {

        // if unset, calculate and set angles of parent to first daughter's
        if (sm.status(*Phi_, symIndex) == CalculationStatus::uncalculated or sm.status(*Theta_, symIndex) == CalculationStatus::uncalculated) {

            // boost daughter momentum from data frame into pc rest frame
            const FourVector<double> p = boost * boosts * model()->fourMomenta()->p(d, daughter);

            auto phi_theta = angles<double>(vect<double>(p), cP);

            // set ambiguous phi to theta
            // todo: in this cases, theta should be 0 or pi. In most cases it is, but sometimes not.
            // Not checking if theta == 0 or pi results in tests passing which would otherwise not
            if (std::isnan(phi_theta[0]))
                phi_theta[0] = phi_theta[1];

            Phi_->setValue(phi_theta[0], d, symIndex, sm);
            Theta_->setValue(phi_theta[1], d, symIndex, sm);

            DEBUG("calculated helicity angles for " << to_string(*daughter)
                  << ": (phi, theta) = (" << phi_theta[0] << ", " << phi_theta[1] << ")");
        }

        // recurse down the decay tree
        calculateAngles(d, daughter, cP, boost, sm);
    }
}

//-------------------------
unsigned HelicityAngles::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    if (pc->isFinalStateParticle())
        throw exceptions::FinalStateParticleCombination("cannot calculate helicity angles for fsp", "HelicityAngles::addParticleCombination");
    return StaticDataAccessor::addParticleCombination(pc);
}

}

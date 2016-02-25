#include "HelicityAngles.h"

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
void HelicityAngles::calculate(DataPoint& d, unsigned dataPartitionIndex)
{
    Phi_->setCalculationStatus(kUncalculated, dataPartitionIndex);
    Theta_->setCalculationStatus(kUncalculated, dataPartitionIndex);

    // call an ISP PC's
    // \todo allow for designating the boost that takes from the data frame to the lab frame (possibly event dependent)
    for (auto& kv : symmetrizationIndices())
        if (kv.first->indices().size() == model()->finalStateParticles().size())
            calculateAngles(d, kv.first, model()->coordinateSystem(), unitMatrix<double, 4>(), dataPartitionIndex);
}

//-------------------------
void HelicityAngles::calculateAngles(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc,
                                     const CoordinateSystem<double, 3>& C, const FourMatrix<double>& boosts,
                                     unsigned dataPartitionIndex)
{
    // terminate recursion
    if (pc->isFinalStateParticle())
        return;

    // get pc's 4-mom in data frame
    const auto P = model()->fourMomenta()->p(d, pc);

    // calculate reference frame for P from parent's RF
    const auto cP = helicityFrame(boosts * P, C);

    // calculate boost from data frame into pc rest frame
    const auto b = lorentzTransformation(-P);

    const unsigned symIndex = symmetrizationIndex(pc);

    for (auto& daughter : pc->daughters()) {

        // boost daughter momentum from data frame into pc rest frame
        const auto p = b * model()->fourMomenta()->p(d, daughter);

        // if unset, set angles of parent to first daughter's
        if (Phi_->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated or
                Theta_->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated ) {

            const auto phi_theta = angles<double>(vect<double>(p), cP);

            Phi_->setValue(phi_theta[0], d, symIndex, dataPartitionIndex);
            Theta_->setValue(phi_theta[1], d, symIndex, dataPartitionIndex);

        }

        // continue down the decay tree
        calculateAngles(d, daughter, cP, b, dataPartitionIndex);
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

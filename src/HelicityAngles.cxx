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
    for (auto& kv : symmetrizationIndices())
        if (kv.first->indices().size() == model()->finalStateParticles().size())
            calculateAngles(d, kv.first, model()->coordinateSystem(),
                            unitMatrix<double, 4>(),
                            //lorentzTransformation(-(model()->fourMomenta().finalStateMomenta(d))),
                            dataPartitionIndex);
}

//-------------------------
void HelicityAngles::calculateAngles(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc,
                                     const CoordinateSystem<double, 3>& C, const FourMatrix<double>& boosts,
                                     unsigned dataPartitionIndex)
{
    // terminate recursion
    if (pc->isFinalStateParticle())
        return;

    // calculate pc momentum in parent frame:
    FourVector<double> P = boosts * model()->fourMomenta()->p(d, pc);

    // calculate reference frame for pc
    CoordinateSystem<double, 3> cP = helicityFrame<double>(P, C);

    // calculate boost from lab frame into pc rest frame
    FourMatrix<double> b = lorentzTransformation<double>(-P) * boosts;

    unsigned symIndex = symmetrizationIndex(pc);

    for (auto& daughter : pc->daughters()) {

        // boost daughter momentum from lab frame into pc rest frame
        FourVector<double> p = b * model()->fourMomenta()->p(d, daughter);

        // if unset, set angles of parent to first daughter's
        if (Phi_->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated or
                Theta_->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated ) {

            const auto phi_theta = angles<double>(vect<double>(p), C);

            Phi_->setValue(phi_theta[0], d, symIndex, dataPartitionIndex);
            Theta_->setValue(phi_theta[1], d, symIndex, dataPartitionIndex);

            //DEBUG("calculated helicity angles: phi = " << phi_theta[0] << ", theta = " << phi_theta[1] << " for " << *pc);
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

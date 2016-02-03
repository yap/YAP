#include "HelicityAngles.h"

#include "Constants.h"
#include "CoordinateSystem.h"
#include "FourVector.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "LorentzTransformation.h"
#include "MathUtilities.h"
#include "ParticleCombination.h"
#include "Rotation.h"
#include "ThreeVector.h"

namespace yap {

//-------------------------
HelicityAngles::HelicityAngles(InitialStateParticle* isp) :
    StaticDataAccessor(isp, &ParticleCombination::equivUpAndDown),
    Phi_(RealCachedDataValue::create(this)),
    Theta_(RealCachedDataValue::create(this))
{
}

//-------------------------
// void HelicityAngles::addSymmetrizationIndex(std::shared_ptr<ParticleCombination> c)
// {
//     /// dFunctions for J == 0 are 0, so we don't need to calculate and store helicity angles
//     if (initialStateParticle()->quantumNumbers().twoJ() == 0
//             and c->parent() == nullptr)
//         return;

//     DataAccessor::addSymmetrizationIndex(c);
// }

//-------------------------
void HelicityAngles::calculate(DataPoint& d, unsigned dataPartitionIndex)
{
    Phi_->setCalculationStatus(kUncalculated, dataPartitionIndex);
    Theta_->setCalculationStatus(kUncalculated, dataPartitionIndex);

    for (auto& pc : initialStateParticle()->particleCombinations())
        calculateAngles(d, pc, initialStateParticle()->coordinateSystem(), unitMatrix<double, 4>(), dataPartitionIndex);
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
    FourVector<double> P = boosts * initialStateParticle()->fourMomenta().p(d, pc);

    // calculate reference frame for pc
    CoordinateSystem<double, 3> cP = helicityFrame<double>(P, C);

    // calculate boost from lab frame into pc rest frame
    FourMatrix<double> b = lorentzTransformation<double>(-P) * boosts;

    unsigned symIndex = symmetrizationIndex(pc);

    for (auto& daughter : pc->daughters()) {

        // boost daughter momentum from lab frame into pc rest frame
        FourVector<double> p = b * initialStateParticle()->fourMomenta().p(d, daughter);

        // if unset, set angles of parent to first daughter's
        if (Phi_->calculationStatus(pc, symIndex, 0) == kUncalculated or
                Theta_->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated ) {

            const auto phi_theta = angles<double>(vect<double>(p), C);
            Phi_->setValue(phi_theta[0], d, symIndex, dataPartitionIndex);
            Theta_->setValue(phi_theta[1], d, symIndex, dataPartitionIndex);

            DEBUG("calculated helicity angles: phi = " << phi_theta[0] << ", theta = " << phi_theta[1] << " for " << *pc);
        }

        // continue down the decay tree
        calculateAngles(d, daughter, cP, b, dataPartitionIndex);
    }
}

}

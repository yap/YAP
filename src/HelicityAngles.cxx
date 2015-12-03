#include "HelicityAngles.h"

#include "Constants.h"
#include "CoordinateSystem.h"
#include "FourVector.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "LorentzTransformation.h"
#include "MathUtilities.h"
#include "Rotation.h"
#include "ThreeVector.h"

namespace yap {

//-------------------------
HelicityAngles::HelicityAngles() :
    StaticDataAccessor(&ParticleCombination::equivUpAndDownButLambda),
    Phi_(new RealCachedDataValue(this)),
    Theta_(new RealCachedDataValue(this))
{
}

//-------------------------
void HelicityAngles::addSymmetrizationIndex(std::shared_ptr<const ParticleCombination> c)
{
    /// dFunctions for J == 0 are 0, so we don't need to calculate and store helicity angles
    if (initialStateParticle()->quantumNumbers().twoJ() == 0
            and c->parent() == nullptr)
        return;

    DataAccessor::addSymmetrizationIndex(c);
}

//-------------------------
void HelicityAngles::calculate(DataPoint& d)
{
    // use a default dataPartitionIndex of 0

    Phi_->setCalculationStatus(kUncalculated, 0);
    Theta_->setCalculationStatus(kUncalculated, 0);

    // calculate initial coordinate system
    CoordinateSystem<double, 3> cISP = helicityFrame<double>(initialStateParticle()->fourMomenta().initialStateMomentum(d),
                                                             initialStateParticle()->coordinateSystem());

    FourMatrix<double> boost = lorentzTransformation<double>(-initialStateParticle()->fourMomenta().initialStateMomentum(d));

    for (auto& pc : initialStateParticle()->particleCombinations())
        calculateAngles(d, pc, cISP, boost);
}

//-------------------------
void HelicityAngles::calculateAngles(DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc, const CoordinateSystem<double, 3>& C, const FourMatrix<double>& boosts)
{
    // calculate boost to put all daughters into pc rest frame
    FourMatrix<double> b = lorentzTransformation<double>(-initialStateParticle()->fourMomenta().p(d, pc)) * boosts;
    
    for (auto& daughter : pc->daughters()) {

        unsigned symIndex = symmetrizationIndex(daughter);
        
        // boost daughter momentum into pc rest frame
        FourVector<double> P = b * initialStateParticle()->fourMomenta().p(d, daughter);

        // set angles if unset
        if (Phi_->calculationStatus(daughter, symIndex, 0) == kUncalculated or
            Theta_->calculationStatus(daughter, symIndex, 0) == kUncalculated ) {

            auto phi_theta = angles<double>(vect<double>(P), C);
            Phi_->setValue(phi_theta[0], d, symIndex, 0);
            Theta_->setValue(phi_theta[1], d, symIndex, 0);
        }

        // if particle is not FSP, call recursively
        if (!pc->isFinalStateParticle())
            calculateAngles(d, pc, helicityFrame<double>(P, C), b);
    }
}

}

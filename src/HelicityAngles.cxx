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
                            //unitMatrix<double, 4>(),
                            lorentzTransformation( -(model()->fourMomenta()->initialStateMomentum(d)) ), // boost into ISP rest frame
                            dataPartitionIndex);
}

//-------------------------
void HelicityAngles::calculateAngles(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc,
                                     const CoordinateSystem<double, 3>& C, const FourMatrix<double>& boosts,
                                     unsigned dataPartitionIndex)
{
    yap::plainLogs(el::Level::Debug);

    // terminate recursion
    if (pc->isFinalStateParticle()) {
        DEBUG(" done calculating Helicity angles");
        return;
    }

    // calculate pc momentum in parent frame:
    const FourVector<double> P = boosts * model()->fourMomenta()->p(d, pc);

    DEBUG(to_string(*pc) << " momentum in parent frame: " << to_string(P) << "; norm = " << norm(vect(P)));

    // calculate reference frame for pc
    const CoordinateSystem<double, 3> cP = helicityFrame<double>(P, C);

    DEBUG("coordinate system: " << to_string(C));
    DEBUG("helicity frame:    " << to_string(cP));

    // calculate boost from lab frame into pc rest frame
    const FourMatrix<double> b = lorentzTransformation<double>(-P) * boosts;

    const unsigned symIndex = symmetrizationIndex(pc);

    for (auto& daughter : pc->daughters()) {

        // boost daughter momentum from lab frame into pc rest frame
        const FourVector<double> p = b * model()->fourMomenta()->p(d, daughter);

        DEBUG(to_string(*daughter) << " momentum in pc rest frame: " << to_string(p));

        // if unset, set angles of parent to first daughter's
        if (Phi_->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated or
                Theta_->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated ) {

            const auto phi_theta = angles<double>(vect<double>(p), C);

            Phi_->setValue(phi_theta[0], d, symIndex, dataPartitionIndex);
            Theta_->setValue(phi_theta[1], d, symIndex, dataPartitionIndex);

            DEBUG("calculated helicity angles: phi = " << phi_theta[0] << ", theta = " << phi_theta[1] << " for " << *pc);
        }

        // continue down the decay tree
        calculateAngles(d, daughter, cP, b, dataPartitionIndex);
    }

    yap::disableLogs(el::Level::Debug);
}

//-------------------------
unsigned HelicityAngles::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    if (pc->isFinalStateParticle())
        throw exceptions::FinalStateParticleCombination("cannot calculate helicity angles for fsp", "HelicityAngles::addParticleCombination");
    return StaticDataAccessor::addParticleCombination(pc);
}

}

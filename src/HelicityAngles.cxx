#include "HelicityAngles.h"

#include "Constants.h"
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
    HelicityAngles_(new CachedDataValue(this, 2))
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

    HelicityAngles_->setCalculationStatus(kUncalculated, 0);

    if (initialStateParticle()->quantumNumbers().twoJ() != 0)
        LOG(ERROR) << "Helicity angles are at the moment only implemented for initial state particles with spin 0.";

    // initial helicity frame.
    // const FourMatrix<double> trans = hfTransform(initialStateLab); // ??
    // boost into RF of initialState
    FourMatrix<double> boost = lorentzTransformation<double>(-initialStateParticle()->fourMomenta().initialStateMomentum(d));
    std::vector<FourVector<double> > finalStatesHf = d.finalStateFourMomenta();
    for (FourVector<double>& v : finalStatesHf)
        v = boost * v;

    for (auto& pc : initialStateParticle()->particleCombinations()) {
        transformDaughters(d, pc, finalStatesHf);
    }

}

//-------------------------
FourMatrix<double> HelicityAngles::hfTransform(const FourVector<double>& daughter)
{
    // the Y axis of the helicity frame
    ThreeVector<double> D = vect<double>(daughter);
    ThreeVector<double> hfY = cross(Axis_Z, D);

    // Rotate to put hfY parallel to Y, and hfZ in the X--Z plane
    ThreeMatrix<double> R1 = rotation(Axis_X, yap::theta(hfY) - PI / 2) * rotation(Axis_Y, PI / 2 - yap::phi(hfY));

    // rotate daughter by R1
    ThreeVector<double> rD = R1 * D;

    // Rotate about hfY (now Y) to put daughter momentum along Z
    ThreeMatrix<double> R2 = rotation(Axis_Y, -signum(rD[0]) * yap::theta(rD));

    // use daughter rotated by R2 for boost
    // to fom helicity frame transformation
    return lorentzTransformation(R2 * R1, -(R2 * rD));
}

//-------------------------
void HelicityAngles::transformDaughters(DataPoint& d,
                                        const std::shared_ptr<const ParticleCombination>& pc,
                                        std::vector<FourVector<double> > finalStatesHf)
{
    // loop over daughters
    for (auto& daugh : pc->daughters()) {

        if (daugh->daughters().empty())
            continue;

        // construct 4-vector of daughter
        FourVector<double> daughter = FourVector_0;
        for (ParticleIndex i : daugh->indices())
            daughter += finalStatesHf.at(i);

        ThreeVector<double> d3 = vect<double>(daughter);

        HelicityAngles_->setValue(0, yap::phi(d3),   d, symmetrizationIndex(daugh));
        HelicityAngles_->setValue(1, yap::theta(d3), d, symmetrizationIndex(daugh));

        HelicityAngles_->setCalculationStatus(kCalculated, symmetrizationIndex(daugh), 0);

        // next helicity frame
        const FourMatrix<double> transDaugh = hfTransform(daughter);
        for (ParticleIndex i : daugh->indices())
            finalStatesHf[i] = transDaugh * finalStatesHf[i];

        transformDaughters(d, daugh, finalStatesHf);
    }
}


}

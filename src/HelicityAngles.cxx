#include "HelicityAngles.h"

#include "Constants.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "LorentzTransformation.h"
#include "MathUtilities.h"
#include "Rotation.h"

#include <TLorentzRotation.h>

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
    //const TLorentzRotation trans = hfTransform(initialStateLab); // ??
    // boost into RF of initialState
    TLorentzRotation boost;
    boost.Boost(-initialStateParticle()->fourMomenta().initialStateMomentum(d).BoostVector());
    std::vector<TLorentzVector> finalStatesHf = d.finalStateFourMomenta();
    for (TLorentzVector& lv : finalStatesHf)
        lv.Transform(boost);

    for (auto& pc : initialStateParticle()->particleCombinations()) {
        transformDaughters(d, pc, finalStatesHf);
    }

}

//-------------------------
TLorentzRotation HelicityAngles::hfTransform(const TLorentzVector& daughterLv)
{
    // code copied from rootpwa
    TLorentzVector daughter = daughterLv;
    const TVector3 zAxisParent(0, 0, 1);  // take z-axis as defined in parent frame
    const TVector3 yHfAxis = zAxisParent.Cross(daughter.Vect());  // y-axis of helicity frame
    // rotate so that yHfAxis becomes parallel to y-axis and zHfAxis ends up in (x, z)-plane
    TRotation rot1;
    rot1.RotateZ(0.5 * PI - yHfAxis.Phi());
    rot1.RotateX(yHfAxis.Theta() - 0.5 * PI);
    daughter *= rot1;
    // rotate about yHfAxis so that daughter momentum is along z-axis
    TRotation rot2;
    rot2.RotateY(-signum(daughter.X()) * daughter.Theta());
    daughter *= rot2;
    // boost to daughter RF
    rot1.Transform(rot2);
    TLorentzRotation hfTransform(rot1);
    hfTransform.Boost(-daughter.BoostVector());
    return hfTransform;
}

//-------------------------
FourMatrix<double> HelicityAngles::hfTransform(const FourVector<double>& daughter)
{
    // the Y axis of the helicity frame
    ThreeVector<double> D = vect<double>(daughter);
    ThreeVector<double> hfY = cross(Axis_Z, D);

    // Rotate to put hfY parallel to Y, and hfZ in the X--Z plane
    ThreeMatrix<double> R1 = rotation<double>(Axis_X, atan2(sqrt(pow(hfY[0], 2) + pow(hfY[1], 2)), hfY[2]) - PI / 2)
                             * rotation<double>(Axis_Y, PI / 2 - atan2(hfY[1], hfY[0]));

    // rotate daughter by R1
    ThreeVector<double> rD = R1 * D;

    // Rotate about hfY (now Y) to put daughter momentum along Z
    ThreeMatrix<double> R2 = rotation<double>(Axis_Y, -signum(rD[0]) * atan2(sqrt(pow(rD[0], 2) + pow(rD[1], 2)), rD[2]));

    // use daughter rotated by R2 for boost
    // to fom helicity frame transformation
    return lorentzTransformation<double>(R2 * R1, -(R2 * rD));
}

//-------------------------
void HelicityAngles::transformDaughters(DataPoint& d,
                                        const std::shared_ptr<const ParticleCombination>& pc,
                                        std::vector<TLorentzVector> finalStatesHf)
{
    // loop over daughters
    for (auto& daugh : pc->daughters()) {

        if (daugh->daughters().empty())
            continue;

        // construct 4-vector of daughter
        TLorentzVector daughter;
        for (ParticleIndex i : daugh->indices())
            daughter += finalStatesHf.at(i);

        double phi = daughter.Phi();
        double theta = daughter.Theta();
        HelicityAngles_->setValue(0, phi,   d, symmetrizationIndex(daugh));
        HelicityAngles_->setValue(1, theta, d, symmetrizationIndex(daugh));

        HelicityAngles_->setCalculationStatus(kCalculated, symmetrizationIndex(daugh), 0);

        // next helicity frame
        const TLorentzRotation transDaugh = hfTransform(daughter);
        for (ParticleIndex i : daugh->indices())
            finalStatesHf.at(i).Transform(transDaugh);

        transformDaughters(d, daugh, finalStatesHf);
    }
}


}

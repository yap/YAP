#include "HelicityAngles.h"

#include "Constants.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "MathUtilities.h"

#include <TLorentzRotation.h>

namespace yap {

//-------------------------
HelicityAngles::HelicityAngles(InitialStateParticle* isp) :
    DataAccessor(isp)
{

}

//-------------------------
void HelicityAngles::calculateHelicityAngles(DataPoint& d)
{
/// \todo Use kinematics::fourMomenta

    if (initialStateParticle()->quantumNumbers().twoJ() != 0)
        LOG(ERROR) << "Helicity angles are at the moment only implemented for initial state particles with spin 0.";


    // final state 4-momenta
    std::vector<TLorentzVector> finalStatesLab;
    for (unsigned i = 0; i < initialStateParticle()->particleCombinations()[0]->indices().size(); ++i) {
        finalStatesLab.push_back(initialStateParticle()->fourMomenta().p(d, i));
    }

    // initial helicity frame. \todo Not sure if correct
    //const TLorentzRotation trans = hfTransform(initialStateLab); // ??
    // boost into RF of initialState
    TLorentzRotation boost;
    boost.Boost(-initialStateParticle()->fourMomenta().initialStateMomentum(d).BoostVector());
    std::vector<TLorentzVector> finalStatesHf = finalStatesLab;
    for (TLorentzVector& lv : finalStatesHf)
        lv.Transform(boost);

    for (std::shared_ptr<ParticleCombination>& pc : initialStateParticle()->particleCombinations()) {
        transformDaughters(pc, finalStatesHf);
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
void HelicityAngles::transformDaughters(std::shared_ptr<ParticleCombination> pc,
                                        std::vector<TLorentzVector> finalStatesHf)
{
    // loop over daughters
    for (const std::shared_ptr<ParticleCombination>& daugh : pc->daughters()) {

        if (daugh->daughters().empty())
            continue;

        // testing ...
        /*std::cout << std::string(*daugh) << "\n";
        // find matching decay channels
        for (unsigned i = 0; i < initialStateParticle()->nChannels(); ++i) {
            std::shared_ptr<yap::SpinAmplitude> sa = initialStateParticle()->channel(i)->spinAmplitude();
            for (std::shared_ptr<ParticleCombination> pcSa : sa->particleCombinations()) {
                if (pcSa == daugh) {
                    std::cout << " matches " << sa << std::string(*sa) << ", " << std::string(*pcSa) << "; ";
                }
            }
        }
        std::cout << "\n";*/


        // construct 4-vector of daughter
        TLorentzVector daughter;
        for (ParticleIndex i : daugh->indices())
            daughter += finalStatesHf.at(i);

        double phi = daughter.Phi();
        double theta = daughter.Theta();
        LOG(DEBUG) << std::string(*daugh) << " helicity angles (phi, theta) = (" << phi << ", " << theta << ")\n";

        // next helicity frame
        const TLorentzRotation transDaugh = hfTransform(daughter);
        for (ParticleIndex i : daugh->indices())
            finalStatesHf.at(i).Transform(transDaugh);

        transformDaughters(daugh, finalStatesHf);
    }
}


}

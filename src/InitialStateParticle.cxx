#include "InitialStateParticle.h"

#include "CanonicalSpinAmplitude.h"
#include "logging.h"

#include <TLorentzRotation.h>

namespace yap {

//-------------------------
InitialStateParticle::InitialStateParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize) :
    DecayingParticle(q, mass, name, radialSize)
{
}

//-------------------------
bool InitialStateParticle::consistent() const
{
    return DecayingParticle::consistent();
}

//-------------------------
void InitialStateParticle::calculateHelicityAngles(DataPoint& d) const
{
    /// \todo put this function somewhere else?

    if (quantumNumbers().twoJ() != 0)
        LOG(ERROR) << "Helicity angles are at the moment only implemented for initial particles with spin 0.";


    // final state 4-momenta
    const std::vector<TLorentzVector>& finalStatesLab = d.fourMomenta();

    // construct 4-vector of initial state
    TLorentzVector initialStateLab;
    for (const TLorentzVector& lv : finalStatesLab)
        initialStateLab += lv;

    // initial helicity frame. \todo Not sure if correct
    const TLorentzRotation trans = CanonicalSpinAmplitude::hfTransform(initialStateLab); // ??
    std::vector<TLorentzVector> finalStatesHf = finalStatesLab;
    for (TLorentzVector& lv : finalStatesHf)
        lv.Transform(trans);

    for (auto& kv : SymmetrizationIndices_) {
        // loop over daughters
        CanonicalSpinAmplitude::transformDaughters(kv.first, finalStatesHf);
    }
}


}

#include "CanonicalSpinAmplitude.h"

#include "Constants.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "MathUtilities.h"
#include "SpinUtilities.h"

namespace yap {

//-------------------------
CanonicalSpinAmplitude::CanonicalSpinAmplitude(const QuantumNumbers& initial, const QuantumNumbers& final1, const QuantumNumbers& final2, unsigned char twoL)
    : SpinAmplitude(initial, final1, final2),
      TwoL_(twoL)
{
    /// \todo put this somewhere else?
    calculateClebschGordanCoefficients();
}

//-------------------------
Amp CanonicalSpinAmplitude::amplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc)
{
    /*if (Recalculate_) {
        std::vector<std::vector<double> >& data = d.data(index());

        // loop over symmetrization indices
        for (auto& kv : SymmetrizationIndices_) {
            std::vector<double>& dataVector = data.at(kv.second);


        }

        Recalculate_ = false;
    }*/


    // \todo implement

    return Complex_0;
}

//-------------------------
bool CanonicalSpinAmplitude::consistent() const
{
    bool consistent = SpinAmplitude::consistent();

    // check angular momentum conservation laws
    int twoJ_P = InitialQuantumNumbers_.twoJ();
    int twoJ_A = FinalQuantumNumbers_[0].twoJ();
    int twoJ_B = FinalQuantumNumbers_[1].twoJ();

    // check if
    // \vect{s_P} = \vect{l} + \vect{s_A} + \vect{s_B}
    bool ok = false;
    for (int twoL_AB = abs(twoJ_A - twoJ_B); twoL_AB <= abs(twoJ_A + twoJ_B); twoL_AB += 2) {
        for (int rhs = abs(TwoL_ - twoL_AB); rhs <= abs(TwoL_ + twoL_AB); rhs += 2) {
            if (twoJ_P == rhs) {
                ok = true;
                break;
            }
            if (ok)
                break;
        }
    }

    if (!ok) {
        LOG(ERROR) << "CanonicalSpinAmplitude::consistent() - angular momentum conservation violated. " <<
                   "J(parent) = " << spinToString(twoJ_P) << "; J(daughter1) = " << spinToString(twoJ_A) << "; J(daughter2) = " << spinToString(twoJ_B) << "; l = " << spinToString(TwoL_);
        consistent =  false;
    }

    if (ClebschGordanCoefficients_.empty()) {
        LOG(ERROR) << "CanonicalSpinAmplitude::consistent() - ClebschGordanCoefficients_ are empty. They are probably all 0 and you can remove this channel.";
        consistent =  false;
    }

    return consistent;
}

//-------------------------
CanonicalSpinAmplitude::operator std::string() const
{
    std::string result = "(l=" + spinToString(TwoL_);

    if (not ClebschGordanCoefficients_.empty()) {
        result += "; λ=";
        auto& last = *(--ClebschGordanCoefficients_.end());
        for (auto& kv : ClebschGordanCoefficients_) {
            result += spinToString(kv.first[0]) + "," + spinToString(kv.first[1]);
            if (&kv != &last)
                result += "; ";
        }
    }

    result += ")";

    return result;
}

//-------------------------
void CanonicalSpinAmplitude::printClebschGordanCoefficients() const
{
    std::cout << "Clebsch-Gordan coefficients for decay: (" << InitialQuantumNumbers_ << ") -> ("
              << FinalQuantumNumbers_[0] << ") + ("
              << FinalQuantumNumbers_[1] << "), " << std::string(*this) << "\n";
    for (auto& kv : ClebschGordanCoefficients_) {
        std::cout << "  λ_1, λ_2 = (" << spinToString(kv.first[0]) << "," << spinToString(kv.first[1])
                  << "): \t" << kv.second << "\n";
    }
}

//-------------------------
void CanonicalSpinAmplitude::calculateClebschGordanCoefficients()
{
    /// code is copied in parts from rootpwa

    const int J  = InitialQuantumNumbers_.twoJ();
    const int s1 = FinalQuantumNumbers_[0].twoJ();
    const int s2 = FinalQuantumNumbers_[1].twoJ();

    for (int lambda1 = -s1; lambda1 <= +s1; lambda1 += 2) {
        for (int lambda2 = -s2; lambda2 <= +s2; lambda2 += 2) {

            // \todo: cross check that S is really meant to be s1 +s2
            const int    S         = s1 + s2;
            const int    lambda    = lambda1 - lambda2;

            // calculate Clebsch-Gordan coefficient for L-S coupling
            const double lsClebsch = clebschGordan(TwoL_, 0, S, lambda, J, lambda);
            if (lsClebsch == 0.)
                continue;

            // calculate Clebsch-Gordan coefficient for S-S coupling
            const double ssClebsch = clebschGordan(s1, lambda1, s2, -lambda2, S, lambda);
            if (ssClebsch == 0.)
                continue;

            LOG(DEBUG) << "Clebsch-Gordan coefficient for λ_1, λ_2 = (" << spinToString(lambda1)
                       << "," << spinToString(lambda2) << "): " << ssClebsch << " * " << lsClebsch
                       << " = " << ssClebsch* lsClebsch << "\n";

            ClebschGordanCoefficients_[ {lambda1, lambda2}] = ssClebsch * lsClebsch;
        }
    }

}

//-------------------------
void CanonicalSpinAmplitude::calculateHelicityAngles(DataPoint& d, std::shared_ptr<InitialStateParticle> initialState)
{

    if (initialState->quantumNumbers().twoJ() != 0)
        LOG(ERROR) << "Helicity angles are at the moment only implemented for initial particles with spin 0.";


    // final state 4-momenta
    const std::vector<TLorentzVector>& finalStatesLab = d.fourMomenta();

    // construct 4-vector of initial state
    TLorentzVector initialStateLab;
    for (const TLorentzVector& lv : finalStatesLab)
        initialStateLab += lv;

    // initial helicity frame. \todo Not sure if correct
    const TLorentzRotation trans = hfTransform(initialStateLab); // ??
    std::vector<TLorentzVector> finalStatesHf = finalStatesLab;
    for (TLorentzVector& lv : finalStatesHf)
        lv.Transform(trans);

    for (std::shared_ptr<ParticleCombination>& pc : initialState->particleCombinations()) {
        transformDaughters(pc, finalStatesHf);
    }
}

//-------------------------
TLorentzRotation CanonicalSpinAmplitude::hfTransform(const TLorentzVector& daughterLv)
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
void CanonicalSpinAmplitude::transformDaughters(std::shared_ptr<ParticleCombination> pc, std::vector<TLorentzVector> finalStatesHf)
{
    // loop over daughters
    for (const std::shared_ptr<ParticleCombination>& daugh : pc->daughters()) {

        if (daugh->daughters().empty())
            continue;

        // construct 4-vector of daughter
        TLorentzVector daughter;
        for (ParticleIndex i : daugh->indices())
            daughter += finalStatesHf.at(i);

        double phi = daughter.Phi();
        double theta = daughter.Theta();
        LOG(DEBUG) << std::string(*daugh) << " helicity angles (phi, theta) = (" << phi << ", " << theta << ")\n";

        // next helicity frame
        const TLorentzRotation transDaugh = CanonicalSpinAmplitude::hfTransform(daughter);
        for (ParticleIndex i : daugh->indices())
            finalStatesHf.at(i).Transform(transDaugh);

        transformDaughters(daugh, finalStatesHf);
    }
}

//-------------------------
bool CanonicalSpinAmplitude::equals(const SpinAmplitude& rhs) const
{
    const CanonicalSpinAmplitude* cSA = dynamic_cast<const CanonicalSpinAmplitude*>(&rhs);
    if (!cSA) return false;

    return (SpinAmplitude::equals(rhs)
            && TwoL_ == cSA->TwoL_);
}

}

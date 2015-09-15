#include "HelicitySpinAmplitude.h"

#include "Constants.h"
#include "DecayingParticle.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "SpinUtilities.h"

namespace yap {

//-------------------------
HelicitySpinAmplitude::HelicitySpinAmplitude(InitialStateParticle* isp, const QuantumNumbers& initial, const QuantumNumbers& final1, const QuantumNumbers& final2, unsigned char twoL)
    : SpinAmplitude(isp, initial, final1, final2),
      TwoL_(twoL)
{
    calculateClebschGordanCoefficients();
}

//-------------------------
Amp HelicitySpinAmplitude::amplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc)
{
    // \todo implement
    return Complex_0;
}

//-------------------------
bool HelicitySpinAmplitude::consistent() const
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
        LOG(ERROR) << "HelicitySpinAmplitude::consistent() - angular momentum conservation violated. " <<
                   "J(parent) = " << spinToString(twoJ_P) << "; J(daughter1) = " << spinToString(twoJ_A) << "; J(daughter2) = " << spinToString(twoJ_B) << "; l = " << spinToString(TwoL_);
        consistent =  false;
    }

    if (ClebschGordanCoefficients_.empty()) {
        LOG(ERROR) << "HelicitySpinAmplitude::consistent() - ClebschGordanCoefficients_ are empty. They are probably all 0 and you can remove this channel.";
        consistent =  false;
    }

    return consistent;
}

//-------------------------
HelicitySpinAmplitude::operator std::string() const
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
void HelicitySpinAmplitude::printClebschGordanCoefficients() const
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
void HelicitySpinAmplitude::calculateClebschGordanCoefficients()
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
bool HelicitySpinAmplitude::equals(const SpinAmplitude& rhs) const
{
    //LOG(DEBUG) << "compare " << std::string(*this) << " and " << std::string(rhs);

    const HelicitySpinAmplitude* cSA = dynamic_cast<const HelicitySpinAmplitude*>(&rhs);
    if (!cSA) return false;

    return (TwoL_ == cSA->TwoL_
            and ClebschGordanCoefficients_ == cSA->ClebschGordanCoefficients_
            and SpinAmplitude::equals(rhs) );
}

}

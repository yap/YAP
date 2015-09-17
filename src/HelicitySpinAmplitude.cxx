#include "HelicitySpinAmplitude.h"

#include "Constants.h"
#include "DecayingParticle.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "SpinUtilities.h"

namespace yap {

//-------------------------
HelicitySpinAmplitude::HelicitySpinAmplitude(InitialStateParticle* isp, const QuantumNumbers& initial,
        const QuantumNumbers& final1, const QuantumNumbers& final2, unsigned char twoL,
        std::pair<std::array<int, 2>, double> clebschGordanCoefficient) :
    SpinAmplitude(isp, initial, final1, final2, twoL),
    ClebschGordanCoefficient_(clebschGordanCoefficient)
{
}

//-------------------------
Amp HelicitySpinAmplitude::calcAmplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc)
{
    // \todo implement

    // needs helicity of parent!!!!
    // DecayChannel needs to sum over daughter helicities

    // make multiple channels, each with definite helicities?

    return Complex_0;
}

//-------------------------
bool HelicitySpinAmplitude::consistent() const
{
    bool consistent = SpinAmplitude::consistent();

    if (ClebschGordanCoefficient_.second == 0) {
        LOG(ERROR) << "HelicitySpinAmplitude::consistent() - ClebschGordanCoefficient is 0.";
        consistent =  false;
    }

    return consistent;
}

//-------------------------
HelicitySpinAmplitude::operator std::string() const
{
    std::string result = "(l=" + spinToString(TwoL_);

    result += "; λ=";
    result += spinToString(ClebschGordanCoefficient_.first[0]) + ",";
    result += spinToString(ClebschGordanCoefficient_.first[1]) + ")";

    return result;
}

//-------------------------
std::map<std::array<int, 2>, double> HelicitySpinAmplitude::calculateClebschGordanCoefficients(
    const QuantumNumbers& initial,
    const QuantumNumbers& final1, const QuantumNumbers& final2,
    unsigned char twoL)
{
    std::map<std::array<int, 2>, double> coeffs;

    /// code is copied in parts from rootpwa

    const int J  = initial.twoJ();
    const int s1 = final1.twoJ();
    const int s2 = final2.twoJ();

    for (int lambda1 = -s1; lambda1 <= +s1; lambda1 += 2) {
        for (int lambda2 = -s2; lambda2 <= +s2; lambda2 += 2) {

            // \todo: cross check that S is really meant to be s1 +s2
            const int    S         = s1 + s2;
            const int    lambda    = lambda1 - lambda2;

            // calculate Clebsch-Gordan coefficient for L-S coupling
            const double lsClebsch = clebschGordan(twoL, 0, S, lambda, J, lambda);
            if (lsClebsch == 0.)
                continue;

            // calculate Clebsch-Gordan coefficient for S-S coupling
            const double ssClebsch = clebschGordan(s1, lambda1, s2, -lambda2, S, lambda);
            if (ssClebsch == 0.)
                continue;

            LOG(DEBUG) << "Clebsch-Gordan coefficient for λ_1, λ_2 = (" << spinToString(lambda1)
                       << "," << spinToString(lambda2) << "): " << ssClebsch << " * " << lsClebsch
                       << " = " << ssClebsch* lsClebsch << "\n";

            coeffs[ {lambda1, lambda2}] = ssClebsch * lsClebsch;
        }
    }

    return coeffs;
}

//-------------------------
bool HelicitySpinAmplitude::equals(const SpinAmplitude& rhs) const
{
    //LOG(DEBUG) << "compare " << std::string(*this) << " and " << std::string(rhs);

    const HelicitySpinAmplitude* cSA = dynamic_cast<const HelicitySpinAmplitude*>(&rhs);
    if (!cSA) return false;

    return (ClebschGordanCoefficient_ == cSA->ClebschGordanCoefficient_
            and SpinAmplitude::equals(rhs) );
}

}

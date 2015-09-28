#include "HelicitySpinAmplitude.h"

#include "Constants.h"
#include "DecayingParticle.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "SpinUtilities.h"
#include "WignerD.h"

namespace yap {

//-------------------------
HelicitySpinAmplitude::HelicitySpinAmplitude(InitialStateParticle* isp, const QuantumNumbers& initial,
        const QuantumNumbers& final1, const QuantumNumbers& final2, unsigned char twoL,
        double clebschGordanCoefficient) :
    SpinAmplitude(isp, initial, final1, final2, twoL),
    ClebschGordanCoefficient_(clebschGordanCoefficient)
{
}

//-------------------------
Amp HelicitySpinAmplitude::calcAmplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc)
{

    /// \todo Take a look at momentum-dependent Clebsch-Gordan coefficients by J. Friedrich and S.U. Chung
    /// implemented in rootPWA by C. Bicker

    // \todo angular normalization factor??? sqrt(2*L + 1)
    Amp a = ClebschGordanCoefficient_;

    const int J = InitialQuantumNumbers_.twoJ();
    // DFunction == 1  for  J == 0
    if (J != 0) {
        const int Lambda  = InitialQuantumNumbers_.twoLambda();
        const int P       = InitialQuantumNumbers_.P();

        const int lambda1 = FinalQuantumNumbers_[0].twoLambda();
        const int lambda2 = FinalQuantumNumbers_[1].twoLambda();
        const int lambda  = lambda1 - lambda2;

        const std::vector<double>& helAngles = initialStateParticle()->helicityAngles().helicityAngles(d, pc);
        const double phi   = helAngles.at(0);  // use daughter1 as analyzer
        const double theta = helAngles.at(1);

        a *= DFunctionConj(J, Lambda, lambda, P, phi, theta);
    }

    DEBUG("HelicitySpinAmplitude = " << a);
    return a;
}

//-------------------------
bool HelicitySpinAmplitude::consistent() const
{
    bool consistent = SpinAmplitude::consistent();

    if (ClebschGordanCoefficient_ == 0) {
        LOG(ERROR) << "HelicitySpinAmplitude::consistent() - ClebschGordanCoefficient is 0.";
        consistent =  false;
    }

    return consistent;
}

//-------------------------
HelicitySpinAmplitude::operator std::string() const
{
    std::string result = "(l=" + spinToString(TwoL_);

    result += "; λ_p=" + spinToString(InitialQuantumNumbers_.twoLambda());
    result += "; λ=";
    result += spinToString(FinalQuantumNumbers_[0].twoLambda()) + ",";
    result += spinToString(FinalQuantumNumbers_[1].twoLambda()) + ")";

    return result;
}

//-------------------------
double HelicitySpinAmplitude::calculateClebschGordanCoefficient(
    const QuantumNumbers& initial,
    const QuantumNumbers& final1, const QuantumNumbers& final2,
    unsigned char twoL)
{
    /// code is copied in parts from rootpwa

    const int J  = initial.twoJ();
    const int s1 = final1.twoJ();
    const int s2 = final2.twoJ();

    int lambda1 = final1.twoLambda();
    int lambda2 = final2.twoLambda();

    // \todo: cross check that S is really meant to be s1 +s2
    const int    S         = s1 + s2;
    const int    lambda    = lambda1 - lambda2;

    // calculate Clebsch-Gordan coefficient for L-S coupling
    const double lsClebsch = clebschGordan(twoL, 0, S, lambda, J, lambda);
    if (lsClebsch == 0.)
        return 0;

    // calculate Clebsch-Gordan coefficient for S-S coupling
    const double ssClebsch = clebschGordan(s1, lambda1, s2, -lambda2, S, lambda);
    if (ssClebsch == 0.)
        return 0;

    DEBUG("Clebsch-Gordan coefficient for λ_1, λ_2 = (" << spinToString(lambda1)
          << "," << spinToString(lambda2) << "): " << ssClebsch << " * " << lsClebsch
          << " = " << ssClebsch * lsClebsch);


    return ssClebsch * lsClebsch;
}

//-------------------------
bool HelicitySpinAmplitude::equals(const SpinAmplitude& rhs) const
{
    //DEBUG("compare " << std::string(*this) << " and " << std::string(rhs));

    const HelicitySpinAmplitude* cSA = dynamic_cast<const HelicitySpinAmplitude*>(&rhs);
    if (!cSA) return false;

    return (ClebschGordanCoefficient_ == cSA->ClebschGordanCoefficient_
            and SpinAmplitude::equals(rhs) );
}

}

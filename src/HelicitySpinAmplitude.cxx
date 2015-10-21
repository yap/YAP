#include "HelicitySpinAmplitude.h"

#include "Constants.h"
#include "DecayingParticle.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "SpinUtilities.h"
#include "WignerD.h"

namespace yap {

//-------------------------
HelicitySpinAmplitude::HelicitySpinAmplitude(const QuantumNumbers& initial,
        const QuantumNumbers& final1, const QuantumNumbers& final2, unsigned char twoL) :
    SpinAmplitude(initial, final1, final2, twoL)
{
}

//-------------------------
std::complex<double> HelicitySpinAmplitude::amplitude(DataPartition& d, std::shared_ptr<const ParticleCombination> pc) const
{
    /// \todo implement
    /*
    if (ClebschGordanCoefficients_.empty()) {
        for (auto& pc : particleCombinations())
            ClebschGordanCoefficients_[pc] = calculateClebschGordanCoefficient(pc);
    }

    /// \todo Take a look at momentum-dependent Clebsch-Gordan coefficients by J. Friedrich and S.U. Chung
    /// implemented in rootPWA by C. Bicker

    std::complex<double> a = ClebschGordanCoefficients_.at(pc);

    // \todo angular normalization factor??? sqrt(2*L + 1)

    const int J = InitialQuantumNumbers_.twoJ();
    // DFunction == 1  for  J == 0
    if (J != 0) {
        const int Lambda  = pc->twoLambda();
        const int P       = InitialQuantumNumbers_.P();

        const int lambda1 = pc->daughters()[0]->twoLambda();
        const int lambda2 = pc->daughters()[1]->twoLambda();
        const int lambda  = lambda1 - lambda2;

        const double phi   = initialStateParticle()->helicityAngles().phi(d.dataPoint(), pc);
        const double theta = initialStateParticle()->helicityAngles().theta(d.dataPoint(), pc);

        a *= DFunctionConj(J, Lambda, lambda, P, phi, theta);
    }

    DEBUG("HelicitySpinAmplitude = " << a);
    return a;*/
    return Complex_0;
}

//-------------------------
bool HelicitySpinAmplitude::consistent() const
{
    bool consistent = SpinAmplitude::consistent();

    /*if (ClebschGordanCoefficients_.empty()) {
        LOG(ERROR) << "HelicitySpinAmplitude::consistent() - ClebschGordanCoefficients are empty.";
        consistent =  false;
    }*/

    return consistent;
}

//-------------------------
HelicitySpinAmplitude::operator std::string() const
{
    std::string result = "(l=" + spinToString(TwoL_) + ")";

    return result;
}

//-------------------------
double HelicitySpinAmplitude::calculateClebschGordanCoefficient(std::shared_ptr<const ParticleCombination> c) const
{
    /// code is copied in parts from rootpwa

    const int J  = InitialQuantumNumbers_.twoJ();
    const int s1 = FinalQuantumNumbers_[0].twoJ();
    const int s2 = FinalQuantumNumbers_[1].twoJ();

    int lambda1 = c->daughters()[0]->twoLambda();
    int lambda2 = c->daughters()[1]->twoLambda();

    const int    S         = s1 + s2;
    const int    lambda    = lambda1 - lambda2;

    // calculate Clebsch-Gordan coefficient for L-S coupling
    const double lsClebsch = clebschGordan(TwoL_, 0, S, lambda, J, lambda);
    if (lsClebsch == 0) {
        //DEBUG("lsClebsch == 0");
        return 0;
    }

    // calculate Clebsch-Gordan coefficient for S-S coupling
    const double ssClebsch = clebschGordan(s1, lambda1, s2, -lambda2, S, lambda);
    if (ssClebsch == 0) {
        //DEBUG("ssClebsch == 0");
        return 0;
    }

    /*DEBUG("Clebsch-Gordan coefficient for λ_1, λ_2 = (" << spinToString(lambda1)
          << "," << spinToString(lambda2) << "): " << ssClebsch << " * " << lsClebsch
          << " = " << ssClebsch * lsClebsch);*/


    return ssClebsch * lsClebsch;

}

//-------------------------
bool HelicitySpinAmplitude::equals(const SpinAmplitude& rhs) const
{
    //DEBUG("compare " << std::string(*this) << " and " << std::string(rhs));

    const HelicitySpinAmplitude* cSA = dynamic_cast<const HelicitySpinAmplitude*>(&rhs);
    if (!cSA) return false;

    return (ClebschGordanCoefficients_ == cSA->ClebschGordanCoefficients_
            and SpinAmplitude::equals(rhs) );
}

}

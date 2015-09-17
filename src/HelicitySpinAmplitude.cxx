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
        double clebschGordanCoefficient) :
    SpinAmplitude(isp, initial, final1, final2, twoL),
    ClebschGordanCoefficient_(clebschGordanCoefficient)
{
}

//-------------------------
Amp HelicitySpinAmplitude::calcAmplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc)
{
    /*const int       Lambda = parent->spinProj();
    const int       P      = parent->P();
    const int       refl   = parent->reflectivity();
    const double    phi    = daughter1->lzVec().Phi();  // use daughter1 as analyzer
    const double    theta  = daughter1->lzVec().Theta();

    // \todo normalization factor???
    return clebschGordanValue() * DFunctionConj();*/
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

    result += "; λ_p=" + spinToString(InitialQuantumNumbers_.twoHelicity());
    result += "; λ=";
    result += spinToString(FinalQuantumNumbers_[0].twoHelicity()) + ",";
    result += spinToString(FinalQuantumNumbers_[1].twoHelicity()) + ")";

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

    int lambda1 = final1.twoHelicity();
    int lambda2 = final2.twoHelicity();

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

    LOG(DEBUG) << "Clebsch-Gordan coefficient for λ_1, λ_2 = (" << spinToString(lambda1)
               << "," << spinToString(lambda2) << "): " << ssClebsch << " * " << lsClebsch
               << " = " << ssClebsch* lsClebsch << "\n";


    return ssClebsch * lsClebsch;
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

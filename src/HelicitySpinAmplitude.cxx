#include "HelicitySpinAmplitude.h"

#include "Constants.h"
#include "DecayingParticle.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "WignerD.h"

namespace yap {

//-------------------------
HelicitySpinAmplitude::HelicitySpinAmplitude(const QuantumNumbers& initial,
        const QuantumNumbers& final1, const QuantumNumbers& final2, unsigned char twoL) :
    SpinAmplitude(initial, final1, final2, twoL),
    SpinAmplitude_(new ComplexCachedDataValue(this))
{
    // SpinAmplitude_ dependency on helicity angles is set in setInitialStateParticle
}

//-------------------------
std::complex<double> HelicitySpinAmplitude::amplitude(DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc, unsigned dataPartitionIndex) const
{
    /// \todo check
    unsigned symIndex = symmetrizationIndex(pc);

    if (SpinAmplitude_->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated) {

        /// \todo Take a look at momentum-dependent Clebsch-Gordan coefficients by J. Friedrich and S.U. Chung
        /// implemented in rootPWA by C. Bicker

        std::complex<double> a = ClebschGordanCoefficients_.at(pc);

        // \todo angular normalization factor??? sqrt(2*L + 1)

        unsigned char twoJ = InitialQuantumNumbers_.twoJ();

        char twoM = pc->twoLambda();

        char twoLambda1 = pc->daughters()[0]->twoLambda();
        char twoLambda2 = pc->daughters()[1]->twoLambda();
        char twoLambda  = twoLambda1 - twoLambda2;

        double phi   = initialStateParticle()->helicityAngles().phi(d, pc);
        double theta = initialStateParticle()->helicityAngles().theta(d, pc);

        a *= std::conj(DFunction(twoJ, twoM, twoLambda, phi, theta, 0));

        SpinAmplitude_->setValue(a, d, symIndex, dataPartitionIndex);

        DEBUG("HelicitySpinAmplitude::amplitude - calculated amplitude for symIndex " << symIndex << " = " << a);
        return a;
    }

    DEBUG("HelicitySpinAmplitude::amplitude - using cached amplitude for symIndex " << symIndex << " = " << SpinAmplitude_->value(d, symIndex));
    return SpinAmplitude_->value(d, symIndex);
}

//-------------------------
bool HelicitySpinAmplitude::consistent() const
{
    bool consistent = SpinAmplitude::consistent();

    if (ClebschGordanCoefficients_.size() != particleCombinations().size()) {
        LOG(ERROR) << "HelicitySpinAmplitude::consistent() - ClebschGordanCoefficients have wrong size.";
        consistent =  false;
    }

    return consistent;
}

//-------------------------
void HelicitySpinAmplitude::addSymmetrizationIndex(std::shared_ptr<const ParticleCombination> c)
{
    ClebschGordanCoefficients_[c] = calculateClebschGordanCoefficient(c);
    SpinAmplitude::addSymmetrizationIndex(c);
}

//-------------------------
void HelicitySpinAmplitude::clearSymmetrizationIndices()
{
    ClebschGordanCoefficients_.clear();
    SpinAmplitude::clearSymmetrizationIndices();
}

//-------------------------
HelicitySpinAmplitude::operator std::string() const
{
    std::string result = "(l=" + spin_to_string(TwoL_) + ")";

    return result;
}

//-------------------------
double HelicitySpinAmplitude::calculateClebschGordanCoefficient(std::shared_ptr<const ParticleCombination> c) const
{
    unsigned char two_j1 = FinalQuantumNumbers_[0].twoJ();
    unsigned char two_j2 = FinalQuantumNumbers_[1].twoJ();
    unsigned char two_S = two_j1 + two_j2;

    char two_lambda1 = c->daughters()[0]->twoLambda();
    char two_lambda2 = c->daughters()[1]->twoLambda();
    char two_lambda = two_lambda1 - two_lambda2;
 
    // calculate Clebsch-Gordan coefficient for L-S coupling
    double lsClebsch = clebschGordan(TwoL_, 0, two_S, two_lambda, InitialQuantumNumbers_.twoJ(), two_lambda);
    if (lsClebsch == 0)
        return 0;
    
    // calculate Clebsch-Gordan coefficient for S-S coupling
    double ssClebsch = clebschGordan(two_j1, two_lambda1, two_j2, -two_lambda2, two_S, two_lambda);
    if (ssClebsch == 0)
        return 0;

    /*DEBUG("Clebsch-Gordan coefficient for λ_1, λ_2 = (" << spin_to_string(lambda1)
          << "," << spin_to_string(lambda2) << "): " << ssClebsch << " * " << lsClebsch
          << " = " << ssClebsch * lsClebsch);*/

    return ssClebsch * lsClebsch;
}

//-------------------------
void HelicitySpinAmplitude::setInitialStateParticle(InitialStateParticle* isp)
{
    SpinAmplitude::setInitialStateParticle(isp);

    if (isp) {
        SpinAmplitude_->addDependency(isp->helicityAngles().phi());
        SpinAmplitude_->addDependency(isp->helicityAngles().theta());
    }
}

//-------------------------
bool HelicitySpinAmplitude::equals(const SpinAmplitude& rhs) const
{
    const HelicitySpinAmplitude* cSA = dynamic_cast<const HelicitySpinAmplitude*>(&rhs);
    if (!cSA) return false;

    return (ClebschGordanCoefficients_ == cSA->ClebschGordanCoefficients_
            and SpinAmplitude::equals(rhs) );
}

}

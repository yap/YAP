#include "HelicitySpinAmplitude.h"

#include "ClebschGordan.h"
#include "Constants.h"
#include "DecayingParticle.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "WignerD.h"

namespace yap {

//-------------------------
HelicitySpinAmplitude::HelicitySpinAmplitude(const QuantumNumbers& initial,
        const QuantumNumbers& final1, const QuantumNumbers& final2, unsigned l) :
    SpinAmplitude(initial, final1, final2, l),
    SpinAmplitude_(new ComplexCachedDataValue(this))
{
}

//-------------------------
std::complex<double> HelicitySpinAmplitude::amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc, unsigned dataPartitionIndex) const
{
    /// \todo check
    unsigned symIndex = symmetrizationIndex(pc);

    if (SpinAmplitude_->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated) {

        /// \todo Take a look at momentum-dependent Clebsch-Gordan coefficients by J. Friedrich and S.U. Chung
        /// implemented in rootPWA by C. Bicker
        std::complex<double> a = ClebschGordanCoefficients_.at(pc);

        // \todo angular normalization factor??? sqrt(2*L + 1)

        unsigned twoJ = initialQuantumNumbers().twoJ();

        int twoM = pc->twoLambda();

        int twoLambda1 = pc->daughters()[0]->twoLambda();
        int twoLambda2 = pc->daughters()[1]->twoLambda();
        int twoLambda  = twoLambda1 - twoLambda2;

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
    bool C = SpinAmplitude::consistent();

    if (ClebschGordanCoefficients_.size() != particleCombinations().size()) {
        FLOG(ERROR) << "ClebschGordanCoefficients have wrong size.";
        C &= false;
    }

    return C;
}

//-------------------------
ParticleCombinationVector HelicitySpinAmplitude::addSymmetrizationIndices(std::shared_ptr<ParticleCombination> pc)
{
    if (!initialStateParticle())
        throw exceptions::Exception("InitialStateParticle unset", "HelicitySpinAmplitude::addSymmetrizationIndices");

    ParticleCombinationVector PCs;

    unsigned two_j1 = finalQuantumNumbers()[0].twoJ();
    unsigned two_j2 = finalQuantumNumbers()[1].twoJ();

    auto d1 = pc->daughters()[0];
    auto d2 = pc->daughters()[1];

    // loop over possible spin projections of d1
    for (int two_lambda1 = -two_j1; two_lambda1 <= (int)two_j1; two_lambda1 += 2) {

        // copy d1, which has unset lambda, into one with set lambda
        auto d1_lambda = initialStateParticle()->particleCombinationCache.copy(*d1, two_lambda1);

        // loop over possible spin projections of d2
        for (int two_lambda2 = -two_j2; two_lambda2 <= (int)two_j2; two_lambda2 += 2) {

            // copy d2, which has unset lambda, into one with set lambda
            auto d2_lambda = initialStateParticle()->particleCombinationCache.copy(*d2, two_lambda2);

            try {
                if (!ClebschGordan::nonzeroCoupling(two_j1, two_lambda1, two_j2, two_lambda2, l(), two_s(), initialQuantumNumbers().twoJ()))
                    continue;
                PCs.push_back(initialStateParticle()->particleCombinationCache.composite({d1_lambda, d2_lambda}));
                addSymmetrizationIndex(PCs.back());
            } catch (const exceptions::InconsistentSpinProjection&) {/*ignore*/}
        }
    }
    return PCs;

    // // if C-G coefficient vanishes, return empty vector
    // FLOG(INFO) << spin_to_string(finalQuantumNumbers()[0].twoJ()) << " " << spin_to_string(pc->daughters()[0]->twoLambda())
    //            << ", "
    //            << spin_to_string(finalQuantumNumbers()[1].twoJ()) << " " << spin_to_string(pc->daughters()[1]->twoLambda())
    //            << " with J = " << spin_to_string(initialQuantumNumbers().twoJ())
    //            << " s = " << spin_to_string(two_s())
    //            << " l = " << l();
    // if (!ClebschGordan::nonzeroCoupling(finalQuantumNumbers()[0].twoJ(), pc->daughters()[0]->twoLambda(),
    //                                     finalQuantumNumbers()[1].twoJ(), pc->daughters()[1]->twoLambda(),
    //                                     l(), two_s(), initialQuantumNumbers().twoJ()))
    //     return PCs;

    // // add all spin projections from -J to J
    // for (int two_lambda = -initialQuantumNumbers().twoJ(); two_lambda <= (int)initialQuantumNumbers().twoJ(); two_lambda += 2) {
    //     auto pc_lambda = std::make_shared<ParticleCombination>(*pc);
    //     pc_lambda->setTwoLambda(two_lambda);
    //     PCs.push_back(initialStateParticle()->particleCombinationCache[pc_lambda]);
    //     addSymmetrizationIndex(PCs.back());
    // }

    // return PCs;
}


//-------------------------
void HelicitySpinAmplitude::addSymmetrizationIndex(std::shared_ptr<ParticleCombination> pc)
{
    // may throw (and be caught in DecayChannel::DecayChannel), so is called before actually adding pc
    auto coupling = ClebschGordan::couple(finalQuantumNumbers()[0].twoJ(), pc->daughters()[0]->twoLambda(),
                                          finalQuantumNumbers()[1].twoJ(), pc->daughters()[1]->twoLambda(),
                                          l(), two_s(), initialQuantumNumbers().twoJ());
    SpinAmplitude::addSymmetrizationIndex(pc);
    ClebschGordanCoefficients_[pc] = coupling;
}

//-------------------------
void HelicitySpinAmplitude::clearSymmetrizationIndices()
{
    ClebschGordanCoefficients_.clear();
    SpinAmplitude::clearSymmetrizationIndices();
}

//-------------------------
bool HelicitySpinAmplitude::equals(const SpinAmplitude& other) const
{
    return SpinAmplitude::equals(other)
           and dynamic_cast<const HelicitySpinAmplitude*>(&other)
           and ClebschGordanCoefficients_ == dynamic_cast<const HelicitySpinAmplitude*>(&other)->ClebschGordanCoefficients_;
}


//-------------------------
void HelicitySpinAmplitude::setInitialStateParticle(InitialStateParticle* isp)
{
    SpinAmplitude::setInitialStateParticle(isp);
    if (!initialStateParticle())
        throw exceptions::Exception("InitialStateParticle unset", "HelicitySpinAmplitude::setInitialStateParticle");
    SpinAmplitude_->addDependency(initialStateParticle()->helicityAngles().phi());
    SpinAmplitude_->addDependency(initialStateParticle()->helicityAngles().theta());
}

}

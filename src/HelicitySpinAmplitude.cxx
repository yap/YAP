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
        const QuantumNumbers& final1,
        const QuantumNumbers& final2,
        unsigned l,
        InitialStateParticle* isp) :
    SpinAmplitude(initial, final1, final2, l, isp),
    SpinAmplitude_(new ComplexCachedDataValue(this))
{
    // set SpinAmplitude_'s dependencies
    if (!initialStateParticle())
        throw exceptions::Exception("InitialStateParticle unset", "HelicitySpinAmplitude::setInitialStateParticle");
    SpinAmplitude_->addDependency(initialStateParticle()->helicityAngles().phi());
    SpinAmplitude_->addDependency(initialStateParticle()->helicityAngles().theta());
}

//-------------------------
void HelicitySpinAmplitude::calculate(DataPoint& d)
{
    // use a default dataPartitionIndex of 0

    Amplitude_->setCalculationStatus(kUncalculated, 0);

    unsigned twoJ = initialQuantumNumbers().twoJ();

    // loop over particle combinations
    auto PCs = particleCombinations();
    for (auto& pc : PCs) {

        unsigned symIndex = symmetrizationIndex(pc);

        // if Amplitude is already set, continue
        if (Amplitude_->calculationStatus(pc, symIndex, 0) != kUncalculated)
            continue;

        int twoM = pc->twoLambda();

        int twoLambda1 = pc->daughters()[0]->twoLambda();
        int twoLambda2 = pc->daughters()[1]->twoLambda();
        int twoLambda  = twoLambda1 - twoLambda2;

        double phi   = initialStateParticle()->helicityAngles().phi(d, pc);
        double theta = initialStateParticle()->helicityAngles().theta(d, pc);

        // calculate D*
        auto a = std::conj(DFunction(twoJ, twoM, twoLambda, phi, theta, 0));
        // multiply by Clebsch-Gordan coefficient
        a *= ClebschGordanCoefficients_.at(pc);

        /// \todo angular normalization factor??? sqrt(2*L + 1)

        // and store
        Amplitude_->setValue(a, d, symIndex, 0);
        DEBUG("HelicitySpinAmplitude::amplitude - calculated amplitude for symIndex " << symIndex << " = " << a);
    }

    /// \todo Take a look at momentum-dependent Clebsch-Gordan
    /// coefficients by J. Friedrich and S.U. Chung implemented in
    /// rootPWA by C. Bicker
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

}

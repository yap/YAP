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
        unsigned l, unsigned two_s,
        InitialStateParticle* isp) :
    SpinAmplitude(initial, final1, final2, l, two_s, isp)
{
    // set cached spin amplitudes' dependencies on helicity angles
    for (auto& a : amplitudeSet()) {
        a->addDependency(initialStateParticle()->helicityAngles().phi());
        a->addDependency(initialStateParticle()->helicityAngles().theta());
    }

    // cache coefficients
    double c = sqrt((2 * L() + 1) / 4 / PI);
    for (int two_m1 = -finalQuantumNumbers()[0].twoJ(); two_m1 <= (int)finalQuantumNumbers()[0].twoJ(); two_m1 += 2)
        for (int two_m2 = -finalQuantumNumbers()[1].twoJ(); two_m2 <= (int)finalQuantumNumbers()[1].twoJ(); two_m2 += 2) {
            double CG = 0;
            try {
                CG = c * ClebschGordan::couple(finalQuantumNumbers()[0].twoJ(), two_m1,
                                               finalQuantumNumbers()[1].twoJ(), two_m2,
                                               L(), twoS(), initialQuantumNumbers().twoJ());
            } catch (const exceptions::InconsistentSpinProjection&) { /* ignore */
                FLOG(INFO) << ClebschGordan::to_string(2 * L(), 0 , twoS(), two_m1 - two_m2, initialQuantumNumbers().twoJ(), two_m1 - two_m2)
                           << " "
                           << ClebschGordan::to_string(finalQuantumNumbers()[0].twoJ(), two_m1, finalQuantumNumbers()[1].twoJ(), -two_m2, twoS(), two_m1 - two_m2)
                           << " = " << CG;
            }

            FLOG(INFO) << ClebschGordan::to_string(2 * L(), 0 , twoS(), two_m1 - two_m2, initialQuantumNumbers().twoJ(), two_m1 - two_m2)
                       << " "
                       << ClebschGordan::to_string(finalQuantumNumbers()[0].twoJ(), two_m1, finalQuantumNumbers()[1].twoJ(), -two_m2, twoS(), two_m1 - two_m2)
                       << " = " << CG;

            if (CG != 0)
                Coefficients_[two_m1][two_m2] = c * CG;
        }

    if (Coefficients_.empty()) {
        FLOG(ERROR) << "no valid nonzero Clebsch-Gordan coefficients stored in " << *this;
        throw exceptions::Exception("no valid nonzero Clebsch-Gordan coefficients stored", "HelicitySpinAmplitude::HelicitySpinAmplitude");
    }
}
//-------------------------
std::complex<double> HelicitySpinAmplitude::calc(int two_M, int two_m1, int two_m2,
        const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
{

    // helicity angles
    double phi   = initialStateParticle()->helicityAngles().phi(d, pc);
    double theta = initialStateParticle()->helicityAngles().theta(d, pc);

    return std::conj(DFunction(initialQuantumNumbers().twoJ(), two_M, two_m1 - two_m2, phi, theta, 0))
           * Coefficients_.at(two_m1).at(two_m2);

    /// \todo Take a look at momentum-dependent Clebsch-Gordan
    /// coefficients by J. Friedrich and S.U. Chung implemented in
    /// rootPWA by C. Bicker
}

}

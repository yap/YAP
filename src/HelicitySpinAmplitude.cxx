#include "HelicitySpinAmplitude.h"

#include "ClebschGordan.h"
#include "Constants.h"
#include "DecayingParticle.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "WignerD.h"

namespace yap {

//-------------------------
HelicitySpinAmplitude::HelicitySpinAmplitude(unsigned two_J, unsigned two_j1, unsigned two_j2, unsigned l, unsigned two_s,
        InitialStateParticle* isp) :
    SpinAmplitude(two_J, two_j1, two_j2, l, two_s, isp)
{
    // set cached spin amplitudes' dependencies on helicity angles
    for (auto& a : amplitudeSet()) {
        a->addDependency(initialStateParticle()->helicityAngles().phi());
        a->addDependency(initialStateParticle()->helicityAngles().theta());
    }

    // cache coefficients
    double c = sqrt((2. * L() + 1) / 4. / PI);
    for (int two_m1 = -finalTwoJ()[0]; two_m1 <= (int)finalTwoJ()[0]; two_m1 += 2)
        for (int two_m2 = -finalTwoJ()[1]; two_m2 <= (int)finalTwoJ()[1]; two_m2 += 2)
            try {

                double CG = c * ClebschGordan::couple(finalTwoJ()[0], two_m1,
                                                      finalTwoJ()[1], two_m2,
                                                      L(), twoS(), initialTwoJ());
                if (CG == 0)
                    continue;

                Coefficients_[two_m1][two_m2] = c * CG;

                // add amplitudes for all initial spin projections
                for (int two_M = -initialTwoJ(); two_M <= (int)initialTwoJ(); two_M += 2)
                    addAmplitude(two_M, two_m1, two_m2);

            } catch (const exceptions::InconsistentSpinProjection&) { /* ignore */ }

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

    return std::conj(DFunction(initialTwoJ(), two_M, two_m1 - two_m2, phi, theta, 0))
           * Coefficients_.at(two_m1).at(two_m2);

    /// \todo Take a look at momentum-dependent Clebsch-Gordan
    /// coefficients by J. Friedrich and S.U. Chung implemented in
    /// rootPWA by C. Bicker
}

}

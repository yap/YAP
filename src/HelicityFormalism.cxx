#include "HelicityFormalism.h"

#include "CachedValue.h"
#include "ClebschGordan.h"
#include "HelicityAngles.h"
#include "logging.h"
#include "Model.h"
#include "Spin.h"
#include "WignerD.h"

namespace yap {

//-------------------------
HelicitySpinAmplitude::HelicitySpinAmplitude(unsigned two_J, const SpinVector& two_j, unsigned l, unsigned two_s) :
    SpinAmplitude(two_J, two_j, l, two_s, equal_by_shared_pointer),
    RequiresHelicityAngles(two_J != 0)
{
    if (finalTwoJ().size() != 2)
        throw exceptions::Exception("Wrong number of daughter spins specified (" + std::to_string(finalTwoJ().size()) + " != 2)",
                                    "HelicitySpinAmplitude::HelicitySpinAmplitude");

    // check j1j2S triangle
    if (!triangle(finalTwoJ()[0], finalTwoJ()[1], twoS()))
        throw exceptions::AngularMomentumNotConserved("HelicitySpinAmplitude::HelicitySpinAmplitude");

    // angular momentum normalization factor
    /// \todo check which is the right one
    // double c = sqrt(2. * L() + 1);
    // double c = sqrt((2. * L() + 1) / 4. / pi<double>() );
    double c  = sqrt((2. * L() + 1) / (initialTwoJ() + 1.));

    // cache coefficients for each spin projection state
    for (const auto& two_m : projections(finalTwoJ()))
        try {
            double CG = ClebschGordan::couple(finalTwoJ()[0], two_m[0], finalTwoJ()[1], two_m[1], L(), twoS(), initialTwoJ());

            if (CG == 0)
                continue;

            Coefficients_[two_m] = c * CG;

            // add amplitudes for all initial spin projections
            for (auto two_M : projections(initialTwoJ()))
                // for J==0, only the Clebsch-Gordan coefficients are returned
                // ==> we don't need storage space in the DataPoint
                addAmplitude(two_M, two_m, initialTwoJ() == 0);

        } catch (const exceptions::InconsistentSpinProjection&) { /* ignore */ }

    if (Coefficients_.empty())
        throw exceptions::Exception("no valid nonzero Clebsch-Gordan coefficients stored", "HelicitySpinAmplitude::HelicitySpinAmplitude");
}

//-------------------------
const std::complex<double> HelicitySpinAmplitude::calc(int two_M, const SpinProjectionVector& two_m,
        const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
{
    // helicity angles
    double phi   = model()->helicityAngles()->phi(d, pc);
    double theta = model()->helicityAngles()->theta(d, pc);

    return std::conj(DFunction(initialTwoJ(), two_M, two_m[0] - two_m[1], phi, theta, 0))
           * Coefficients_.at(two_m);

    /// \todo Take a look at momentum-dependent Clebsch-Gordan
    /// coefficients by J. Friedrich and S.U. Chung implemented in
    /// rootPWA by C. Bicker
}

}

#include "HelicityFormalism.h"

#include "ClebschGordan.h"
#include "HelicityAngles.h"
#include "logging.h"
#include "Model.h"
#include "spin.h"
#include "WignerD.h"

namespace yap {

//-------------------------
HelicitySpinAmplitude::HelicitySpinAmplitude(unsigned two_J, unsigned two_j1, unsigned two_j2, unsigned l, unsigned two_s, ParticleCombination::Equiv* equiv) :
    SpinAmplitude(two_J, two_j1, two_j2, l, two_s, equiv), RequiresHelicityAngles()
{
    // angular momentum normalization factor
    /// \todo check which is the right one
    // double c = sqrt(2. * L() + 1);
    //double c = sqrt((2. * L() + 1) / 4. / pi<double>() );
    double c  = sqrt((2. * L() + 1) / (initialTwoJ() + 1.));

    // cache coefficients
    for (int two_m1 = -finalTwoJ()[0]; two_m1 <= (int)finalTwoJ()[0]; two_m1 += 2)
        for (int two_m2 = -finalTwoJ()[1]; two_m2 <= (int)finalTwoJ()[1]; two_m2 += 2)
            try {
                double CG = ClebschGordan::couple(finalTwoJ()[0], two_m1, finalTwoJ()[1], two_m2, L(), twoS(), initialTwoJ());

                if (CG == 0)
                    continue;

                Coefficients_[two_m1][two_m2] = c * CG;

                // add amplitudes for all initial spin projections
                for (int two_M = -initialTwoJ(); two_M <= (int)initialTwoJ(); two_M += 2)
                    addAmplitude(two_M, two_m1, two_m2);

            } catch (const exceptions::InconsistentSpinProjection&) { /* ignore */ }

    if (Coefficients_.empty())
        throw exceptions::Exception("no valid nonzero Clebsch-Gordan coefficients stored", "HelicitySpinAmplitude::HelicitySpinAmplitude");
}

//-------------------------
void HelicitySpinAmplitude::setDependencies(std::shared_ptr<CachedDataValue> a)
{
    a->addDependency(model()->helicityAngles()->phi());
    a->addDependency(model()->helicityAngles()->theta());
}

//-------------------------
std::complex<double> HelicitySpinAmplitude::calc(int two_M, int two_m1, int two_m2,
        const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
{

    // helicity angles
    double phi   = model()->helicityAngles()->phi(d, pc);
    double theta = model()->helicityAngles()->theta(d, pc);

    DEBUG("HelicitySpinAmplitude::calc : "
          << spin_to_string(initialTwoJ()) << ", " << spin_to_string(two_M)
          << " -> " << spin_to_string(two_m1) << " + " << spin_to_string(two_m2)
          << " for pc = " << *pc << " with CG coeff " << Coefficients_.at(two_m1).at(two_m2)
          << " and helicity angles (" << phi << ", " << theta << ")");

    return std::conj(DFunction(initialTwoJ(), two_M, two_m1 - two_m2, phi, theta, 0))
           * Coefficients_.at(two_m1).at(two_m2);;

    /// \todo Take a look at momentum-dependent Clebsch-Gordan
    /// coefficients by J. Friedrich and S.U. Chung implemented in
    /// rootPWA by C. Bicker
}

}

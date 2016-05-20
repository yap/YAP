#include "ZemachFormalism.h"

#include "CachedDataValue.h"
#include "Constants.h"
#include "Exceptions.h"
#include "FourMomenta.h"
#include "logging.h"
#include "LorentzTransformation.h"
#include "Model.h"
#include "spin.h"

namespace yap {

//-------------------------
ZemachSpinAmplitude::ZemachSpinAmplitude(unsigned two_J, unsigned two_j1, unsigned two_j2, unsigned l, unsigned two_s) :
    SpinAmplitude(two_J, two_j1, two_j2, l, two_s, ParticleCombination::equivZemach)
{
    if (is_odd(two_J) or is_odd(two_j1) or is_odd(two_j2) or is_odd(two_s))
        throw exceptions::Exception("only supporting integer spins", "ZemachSpinAmplitude::ZemachSpinAmplitude");

    if (initialTwoJ() != 0) {
        if (finalTwoJ()[0] != 0 or finalTwoJ()[1] != 0)
            throw exceptions::Exception("Zemach not valid for " + spin_to_string(initialTwoJ())
                                        + " -> " + spin_to_string(finalTwoJ()[0])
                                        + " + " + spin_to_string(finalTwoJ()[1]),
                                        "ZemachSpinAmplitude::ZemachSpinAmplitude");
        if (2 * L() != initialTwoJ() and twoS() != 0)
            throw exceptions::Exception("Zemach not valid for " + spin_to_string(initialTwoJ())
                                        + " -> " + spin_to_string(finalTwoJ()[0])
                                        + " + " + spin_to_string(finalTwoJ()[1])
                                        + " with L = " + std::to_string(L()) + " and S = " + spin_to_string(twoS()),
                                        "ZemachSpinAmplitude::ZemachSpinAmplitude");
    } else {

        if (finalTwoJ()[0] != 0 and finalTwoJ()[1] != 0)
            throw exceptions::Exception("one daughter must be a scalar",
                                        "ZemachSpinAmplitude::ZemachSpinAmplitude");

        // get resonance spin
        auto two_J_r = std::max(finalTwoJ()[0], finalTwoJ()[1]);
        if (two_J_r != 2 * L())
            throw exceptions::Exception("orbital angular momentum != resonance spin (" + std::to_string(L()) + " != " + spin_to_string(two_J_r) + ")",
                                        "ZemachSpinAmplitude::ZemachSpinAmplitude");
        if (two_J_r != twoS())
            throw exceptions::Exception("total spin != resonance spin (" + std::to_string(L()) + " != " + spin_to_string(twoS()) + ")",
                                        "ZemachSpinAmplitude::ZemachSpinAmplitude");

    }

    if (twoS() > 2 * 2)
        throw exceptions::Exception("currently only supporting up to spin 2",
                                    "ZemachSpinAmplitude::ZemachSpinAmplitude");

    addAmplitude(0, 0, 0);
    // dependencies?
}

//-------------------------
void ZemachSpinAmplitude::setDependencies(std::shared_ptr<CachedDataValue> a)
{
    a->addDependency(DaughterCachedDataValue(model()->fourMomenta()->momentum(), 0));
    a->addDependency(DaughterCachedDataValue(model()->fourMomenta()->momentum(), 1));
}

//-------------------------
std::complex<double> ZemachSpinAmplitude::calc(int two_M, int two_m1, int two_m2,
        const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
{
    if (pc->indices().size() < 3)
        return Complex_1;

    if (pc->indices().size() > 3)
        throw exceptions::Exception("Zemach not valid for more than 3 particles",
                                    "ZemachSpinAmplitude::calc");

    if (twoS() == 0)
        return Complex_1;

    if (twoS() > 4)
        throw exceptions::Exception("spin greater than 2 not supported", "ZemachSpinAmplitude::calc");

    // find resonance and spectator ParticleCombinations
    auto pcR = pc->daughters()[0];
    auto pcS = pc->daughters()[1];
    if (pcR->indices().size() == 1)
        std::swap(pcR, pcS);

    // get resonance four-momentum
    auto R4 = model()->fourMomenta()->p(d, pcR);
    // get spectator four-momentum
    auto p4 = model()->fourMomenta()->p(d, pcS);
    // get four-momentum of R's daughter
    auto q4 = model()->fourMomenta()->p(d, pcR->daughters()[0]);

    FDEBUG("R4 = " << R4 << " = " << model()->fourMomenta()->p(d, pcR->daughters()[0]) << " + " << model()->fourMomenta()->p(d, pcR->daughters()[1]));
    FDEBUG("p4 = " << p4);
    FDEBUG("q4 = " << q4);

    // boost p and q into R rest frame, and get spacial components
    auto B = lorentzTransformation<double>(-R4);
    auto p = vect(B * p4);
    auto q = vect(B * q4);

    if (twoS() == 2) {
        FDEBUG("p = " << to_string(p));
        FDEBUG("q = " << to_string(q));
        FDEBUG("p * q = " << (p * q));
        return -2. * (p * q) * Complex_1;
    }

    // else twoS() == 4
    return 4 * (pow(p * q, 2) - norm(p) * norm(q) / 3.) * Complex_1;
}

}

#include "CanonicalSpinAmplitude.h"

#include "Constants.h"
#include "logging.h"
#include "MathUtilities.h"
#include "SpinUtilities.h"

namespace yap {

//-------------------------
CanonicalSpinAmplitude::CanonicalSpinAmplitude(const QuantumNumbers& initial, const QuantumNumbers& final1, const QuantumNumbers& final2, unsigned char twoL)
    : SpinAmplitude(initial, final1, final2),
      TwoL_(twoL)
{}

//-------------------------
Amp CanonicalSpinAmplitude::amplitude(DataPoint& d)
{
    // \todo implement
    return Complex_0;
}

//-------------------------
bool CanonicalSpinAmplitude::consistent() const
{
    bool consistent = SpinAmplitude::consistent();

    // check angular momentum conservation laws
    int twoJ_P = InitialQuantumNumbers_.twoJ();
    int twoJ_A = FinalQuantumNumbers_[0].twoJ();
    int twoJ_B = FinalQuantumNumbers_[1].twoJ();

    // check if
    // \vect{s_P} = \vect{l} + \vect{s_A} + \vect{s_B}
    bool ok = false;
    for (int twoL_AB = abs(twoJ_A - twoJ_B); twoL_AB <= abs(twoJ_A + twoJ_B); twoL_AB += 2) {
        for (int rhs = abs(TwoL_ - twoL_AB); rhs <= abs(TwoL_ + twoL_AB); rhs += 2) {
            if (twoJ_P == rhs) {
                ok = true;
                break;
            }
            if (ok)
                break;
        }
    }

    if (!ok) {
        LOG(ERROR) << "CanonicalSpinAmplitude::consistent() - angular momentum conservation violated. " <<
                   "J(parent) = " << spinToString(twoJ_P) << "; J(daughterA) = " << spinToString(twoJ_A) << "; J(daughterB) = " << spinToString(twoJ_B) << "; l = " << spinToString(TwoL_);
        consistent =  false;
    }

    return consistent;
}

//-------------------------
CanonicalSpinAmplitude::operator std::string() const
{
    return "(l=" + spinToString(TwoL_) + ")";
}

//-------------------------
void CanonicalSpinAmplitude::calculateClebschGordanCoefficients()
{
    /// code is copied in parts from rootpwa

    const int J  = InitialQuantumNumbers_.twoJ();
    const int s1 = FinalQuantumNumbers_[0].twoJ();
    const int s2 = FinalQuantumNumbers_[1].twoJ();

    for (int lambda1 = -s1; lambda1 <= +s1; lambda1 += 2) {
        for (int lambda2 = -s2; lambda2 <= +s2; lambda2 += 2) {

            // \todo
            const int    S         = 0;//vertex->S();
            const int    lambda    = lambda1 - lambda2;

            // calculate Clebsch-Gordan coefficient for L-S coupling
            const double lsClebsch = clebschGordan(TwoL_, 0, S, lambda, J, lambda);

            // calculate Clebsch-Gordan coefficient for S-S coupling
            const double ssClebsch = clebschGordan(s1, lambda1, s2, -lambda2, S, lambda);

        }
    }


}

//-------------------------
bool CanonicalSpinAmplitude::equals(const SpinAmplitude& rhs) const
{
    const CanonicalSpinAmplitude* cSA = dynamic_cast<const CanonicalSpinAmplitude*>(&rhs);
    if (!cSA) return false;

    return (SpinAmplitude::equals(rhs)
            && TwoL_ == cSA->TwoL_);
}

}

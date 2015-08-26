#include "CanonicalSpinAmplitude.h"

#include "Constants.h"
#include "logging.h"

namespace yap {

//-------------------------
CanonicalSpinAmplitude::CanonicalSpinAmplitude(const QuantumNumbers& initial, const QuantumNumbers& final1, const QuantumNumbers& final2, unsigned char l)
    : SpinAmplitude(initial, final1, final2),
      L_(l)
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
    int twoL = 2 * L_;
    int twoJ_P = InitialQuantumNumbers_.twoJ();
    int twoJ_A = FinalQuantumNumbers_[0].twoJ();
    int twoJ_B = FinalQuantumNumbers_[1].twoJ();

    // check if
    // \vect{s_P} = \vect{l} + \vect{s_A} + \vect{s_B}
    bool ok = false;
    for (int twoL_AB = abs(twoJ_A - twoJ_B); twoL_AB <= abs(twoJ_A + twoJ_B); twoL_AB += 2) {
        for (int rhs = abs(twoL - twoL_AB); rhs <= abs(twoL + twoL_AB); rhs += 2) {
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
                   "J(parent) = " << .5 * twoJ_P << "; J(daughterA) = " << .5 * twoJ_A << "; J(daughterB) = " << .5 * twoJ_B << "; l = " << .5 * twoL;
        consistent =  false;
    }

    return consistent;
}

//-------------------------
CanonicalSpinAmplitude::operator std::string() const
{
    return "(l=" + std::to_string(static_cast<unsigned>(L_)) + ")";
}

//-------------------------
bool CanonicalSpinAmplitude::equals(const SpinAmplitude& rhs) const
{
    const CanonicalSpinAmplitude* cSA = dynamic_cast<const CanonicalSpinAmplitude*>(&rhs);
    if (!cSA) return false;

    return (SpinAmplitude::equals(rhs)
            && L_ == cSA->L_);
}

}

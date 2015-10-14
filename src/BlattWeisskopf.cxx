#include "BlattWeisskopf.h"

#include "Constants.h"
#include "DecayChannel.h"
#include "DecayingParticle.h"
#include "logging.h"
#include "Resonance.h"
#include "SpinAmplitude.h"
#include "SpinUtilities.h"

namespace yap {

//-------------------------
BlattWeisskopf::BlattWeisskopf(DecayChannel* decayChannel) :
    DecayChannel_(decayChannel),
    Value_(new CachedValue())
{
    Value_->addDependency(DecayChannel_->parent()->radialSize());
    Value_->addDependency(DecayChannel_->breakupMomentum());
}

void BlattWeisskopf::precalculate()
{
    if (Value_->calculationStatus() == kUncalculated) {

        /// \todo What if we want to fit masses?
        double breakupMom = DecayChannel_->breakupMomentum()->value().real();
        unsigned twoL = DecayChannel_->spinAmplitude()->twoL();
        double R = DecayChannel_->parent()->radialSize()->realValue();

        /// \todo ? in denominator, there must be q^2_(ab), the breakup momentum calculated from the invariant masses

        //
        // The following code was copied from rootPWA
        //
        const double z   = (breakupMom * breakupMom) * (R * R);
        double       bf2 = 0;
        switch (twoL) {
            case 0:  // L = 0
                bf2 = 1;
                break;
            case 2:  // L = 1
                bf2 = (2 * z) / (z + 1);
                break;
            case 4:  // L = 2
                bf2 = (13 * z * z) / (z * (z + 3) + 9);
                break;
            case 6:  // L = 3
                bf2 = (277 * z * z * z) / (z * (z * (z + 6) + 45) + 225);
                break;
            case 8: { // L = 4
                const double z2 = z * z;
                bf2 = (12746 * z2 * z2) / (z * (z * (z * (z + 10) + 135) + 1575) + 11025);
            }
            break;
            case 10: { // L = 5
                const double z2 = z * z;
                bf2 = (998881 * z2 * z2 * z)
                      / (z * (z * (z * (z * (z + 15) + 315) + 6300) + 99225) + 893025);
            }
            break;
            case 12: { // L = 6
                const double z3 = z * z * z;
                bf2 = (118394977 * z3 * z3)
                      / (z * (z * (z * (z * (z * (z + 21) + 630) + 18900) + 496125) + 9823275) + 108056025);
            }
            break;
            case 14: { // L = 7
                const double z3 = z * z * z;
                bf2 = (19727003738LL * z3 * z3 * z)
                      / (z * (z * (z * (z * (z * (z * (z + 28) + 1134) + 47250) + 1819125) + 58939650)
                              + 1404728325L) + 18261468225LL);
            }
            break;
            default:
                LOG(ERROR) << "calculation of Blatt-Weisskopf barrier factor is not (yet) implemented for L = "
                           << spinToString(twoL) << ". returning 0." << std::endl;
                Value_->setValue(Complex_0);
        }

        Value_->setValue(sqrt(bf2));

        DEBUG("Blatt-Weisskopf barrier factor (L = " << spinToString(twoL) << ", " << "q = " << breakupMom << " GeV/c; R = " << R << " 1/GeV) = " << Value_->value());
    }
}

//-------------------------
bool BlattWeisskopf::consistent() const
{
    return true;
}

}


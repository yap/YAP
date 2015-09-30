#include "BlattWeisskopf.h"

#include "Constants.h"
#include "DecayChannel.h"
#include "DecayingParticle.h"
#include "logging.h"
#include "SpinAmplitude.h"
#include "SpinUtilities.h"

namespace yap {

//-------------------------
BlattWeisskopf::BlattWeisskopf(DecayChannel* decayChannel) :
    AmplitudeComponent(),
    DecayChannel_(decayChannel),
    CachedAmplitude_(0),
    CalculationStatus_(kUncalculated)
{}

//-------------------------
const Amp& BlattWeisskopf::amplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc)
{
    if (CalculationStatus_ == kUncalculated) {

        /// \todo What if we want to fit masses?
        double breakupMom = DecayChannel_->breakupMomentum();
        unsigned twoL = DecayChannel_->spinAmplitude()->twoL();
        double R = DecayChannel_->parent()->radialSize();

        /// \todo ? in denominator, there must be q^2_(ab), the breakup momentum calculated from the invariant masses

        /**
         * The following code was copied from rootPWA
         */
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
                CachedAmplitude_ = Complex_0;
        }

        CachedAmplitude_ = sqrt(bf2);
        CalculationStatus_ = kCalculated;

        DEBUG("Blatt-Weisskopf barrier factor (L = " << spinToString(twoL) << ", " << "q = " << breakupMom << " GeV/c; R = " << R << " 1/GeV) = " << CachedAmplitude_);

    }

    return CachedAmplitude_;
}

//-------------------------
bool BlattWeisskopf::consistent() const
{
    return true;
}

}


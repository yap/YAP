#include "BlattWeisskopf.h"

#include "Constants.h"
#include "DecayChannel.h"
#include "DecayingParticle.h"
#include "HelicitySpinAmplitude.h"
#include "logging.h"
#include "SpinUtilities.h"

namespace yap {

//-------------------------
BlattWeisskopf::BlattWeisskopf(DecayChannel* decayChannel) :
    AmplitudeComponentDataAccessor(decayChannel->initialStateParticle()),
    DecayChannel_(decayChannel)
{}

//-------------------------
Amp BlattWeisskopf::calcAmplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc)
{
    double breakupMom = DecayChannel_->breakupMomentum();
    if (not std::dynamic_pointer_cast<HelicitySpinAmplitude>(DecayChannel_->spinAmplitude())) {
        LOG(ERROR) << "Blatt-Weisskopf barrier factor can only be calculated with a  HelicitySpinAmplitude.";
        return Complex_0;
    }

    unsigned twoL = std::static_pointer_cast<HelicitySpinAmplitude>(DecayChannel_->spinAmplitude())->twoL();

    /// \todo confirm
    // const double Pr    = 0.1973)  // momentum scale 0.1973 GeV/c corresponds to 1 fm interaction radius
    double Pr = 1. / DecayChannel_->parent()->radialSize();

    /**
     * The following code was copied from rootPWA
     */
    const double z   = (breakupMom * breakupMom) / (Pr * Pr);
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
            return 0;
    }

    LOG(DEBUG) << "Blatt-Weisskopf barrier factor (L = " << spinToString(twoL) << ", " << "q = " << breakupMom << " GeV/c; P_r = " << Pr << " GeV/c) = " << sqrt(bf2) << std::endl;

    return sqrt(bf2);
}

//-------------------------
bool BlattWeisskopf::consistent() const
{
    return AmplitudeComponentDataAccessor::consistent();
}

}


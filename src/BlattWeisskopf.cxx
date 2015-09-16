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
    AmplitudeComponent(),
    DataAccessor(decayChannel->initialStateParticle()),
    DecayChannel_(decayChannel)
{}

//-------------------------
Amp BlattWeisskopf::amplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc)
{
    // find symmetrization index
    unsigned sym_index = SymmetrizationIndices_[pc];

    // check whether data-dependent calculation needs to be made
    CalculationStatus& calcStat = CalculationStatuses(d, sym_index);
    if (calcStat == kUncalculated) {

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
        /*if (debug)
            LOG(ERROR) << "squared Blatt-Weisskopf barrier factor(L = " << spinQn(twoL) << ", "
                       << "q = " << maxPrecision(breakupMom) << " GeV/c; P_r = " << Pr << " GeV/c) = "
                       << maxPrecision(bf2) << std::endl;*/

        double bf = sqrt(bf2);

        // store into data point
        data(d, sym_index) = {bf};

        // set calculation status
        calcStat = kCalculated;

        return bf;
    }

    // else return from cached data
    return Amp(data(d, sym_index)[0]);
}

//-------------------------
bool BlattWeisskopf::consistent() const
{
    /// \todo implement
    return true;
}

}


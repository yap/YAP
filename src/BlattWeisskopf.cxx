#include "BlattWeisskopf.h"

#include "Constants.h"
#include "DecayChannel.h"
#include "DecayingParticle.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "Resonance.h"
#include "SpinAmplitude.h"
#include "SpinUtilities.h"

namespace yap {

//-------------------------
BlattWeisskopf::BlattWeisskopf(DecayChannel* decayChannel) :
    DataAccessor(&ParticleCombination::equivDownByOrderlessContent),
    DecayChannel_(decayChannel),
    Fq_r(new RealCachedDataValue(this)),
    Fq_ab(new RealCachedDataValue(this))
{
    Fq_r->addDependency(DecayChannel_->parent()->mass());
    Fq_r->addDependency(DecayChannel_->parent()->radialSize());

    Fq_ab->addDependency(DecayChannel_->parent()->radialSize());

    /// measured breakup momenta and four momenta dependencies are set in setInitialStateParticle
}

//-------------------------
std::complex<double> BlattWeisskopf::amplitude(DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc, unsigned dataPartitionIndex) const
{
    unsigned symIndex = symmetrizationIndex(pc);
    bool calc(false); // for debugging

    if (Fq_r->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated) {
        // nominal breakup momentum
        double m2_R = pow(DecayChannel_->parent()->mass()->value(), 2);
        double m_a = initialStateParticle()->fourMomenta().m(d, pc->daughters().at(0));
        double m_b = initialStateParticle()->fourMomenta().m(d, pc->daughters().at(1));
        double q2 = MeasuredBreakupMomenta::calcQ2(m2_R, m_a, m_b);

        DEBUG(std::string(*pc->daughters().at(0)) << " " << std::string(*pc->daughters().at(1)));
        //DEBUG(initialStateParticle()->fourMomenta().symmetrizationIndex(pc->daughters().at(0)) <<
        //        " " << initialStateParticle()->fourMomenta().symmetrizationIndex(pc->daughters().at(1)));
        DEBUG(m_a << " " << m_b);

        double R = DecayChannel_->parent()->radialSize()->value();
        double f = sqrt(F2(DecayChannel_->spinAmplitude()->twoL(), R * R, q2));
        Fq_r->setValue(f, d, symIndex, dataPartitionIndex);

        calc = true;
        DEBUG("BlattWeisskopf::amplitude - calculated barrier factor Fq_r (L = " << spinToString(DecayChannel_->spinAmplitude()->twoL()) << ") = " << Fq_r->value(d, symIndex));
    }

    if (Fq_ab->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated) {
        // measured breakup momentum
        double q2 = initialStateParticle()->measuredBreakupMomenta().q2(d, pc);

        double R = DecayChannel_->parent()->radialSize()->value();
        double f = sqrt(F2(DecayChannel_->spinAmplitude()->twoL(), R * R, q2));
        Fq_ab->setValue(f, d, symIndex, dataPartitionIndex);

        calc = true;
        DEBUG("BlattWeisskopf::amplitude - calculated barrier factor Fq_ab (L = " << spinToString(DecayChannel_->spinAmplitude()->twoL()) << ") = " << Fq_ab->value(d, symIndex));
    }

    double Fq_rOFq_ab = Fq_r->value(d, symIndex) / Fq_ab->value(d, symIndex);

    if (calc) {
        DEBUG("BlattWeisskopf::amplitude - using calculated values to calculate Blatt-Weisskopf barrier factor ratio (L = " << spinToString(DecayChannel_->spinAmplitude()->twoL()) << ") = " << Fq_rOFq_ab);
    } else {
        DEBUG("BlattWeisskopf::amplitude - using cached values to calculate Blatt-Weisskopf barrier factor ratio (L = " << spinToString(DecayChannel_->spinAmplitude()->twoL()) << ") = " << Fq_rOFq_ab);
    }

    return std::complex<double>(Fq_rOFq_ab);
}

//-------------------------
bool BlattWeisskopf::consistent() const
{
    return true;
}

//-------------------------
double BlattWeisskopf::F2(int twoL, double R2, double q2)
{
    const double z   = q2 * R2;
    switch (twoL) {
        case 0:  // L = 0
            return 1.;
        case 2:  // L = 1
            return 1. + z;
        case 4:  // L = 2
            return 9. + 3.*z + z * z;
        default:
            /// \todo put in generic formula for L > 2
            LOG(ERROR) << "calculation of Blatt-Weisskopf barrier factor is not (yet) implemented for L = "
                       << spinToString(twoL) << ". returning 0." << std::endl;
            return 0;
    }
}

//-------------------------
void BlattWeisskopf::setInitialStateParticle(InitialStateParticle* isp)
{
    DataAccessor::setInitialStateParticle(isp);

    if (initialStateParticle()) {
        Fq_r->addDependency(initialStateParticle()->fourMomenta().masses());
        Fq_ab->addDependency(initialStateParticle()->measuredBreakupMomenta().breakupMomenta());
    }
}

}


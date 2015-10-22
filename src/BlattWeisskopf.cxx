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
    DataAccessor(&ParticleCombination::equivDownByOrderlessContent),
    DecayChannel_(decayChannel),
    Fq_r(new RealCachedDataValue(this)),
    Fq_ab(new RealCachedDataValue(this))
{
    Fq_r->addDependency(DecayChannel_->parent()->mass());
    Fq_r->addDependency(DecayChannel_->parent()->radialSize());

    Fq_ab->addDependency(DecayChannel_->parent()->radialSize());

    /// measured breakup momenta and four momenta never change, so no need to add them here
}

//-------------------------
std::complex<double> BlattWeisskopf::amplitude(DataPartition& d, std::shared_ptr<const ParticleCombination> pc) const
{
    /// \todo check
    unsigned symIndex = symmetrizationIndex(pc);

    if (Fq_r->calculationStatus(pc, symIndex, d.index()) == kUncalculated) {
        // nominal breakup momentum
        double m2_R = DecayChannel_->parent()->mass()->value();
        double m_a = initialStateParticle()->fourMomenta().m(d.dataPoint(), pc->daughters().at(0));
        double m_b = initialStateParticle()->fourMomenta().m(d.dataPoint(), pc->daughters().at(1));
        double q2 = MeasuredBreakupMomenta::calcQ2(m2_R, m_a, m_b);

        double R = DecayChannel_->parent()->radialSize()->value();
        double f = sqrt(F2(DecayChannel_->spinAmplitude()->twoL(), R*R, q2));
        Fq_r->setValue(f, d.dataPoint(), symIndex, d.index());

        DEBUG("Blatt-Weisskopf barrier factor Fq_r (L = " << spinToString(DecayChannel_->spinAmplitude()->twoL()) << ") = " << Fq_r->value(d.dataPoint(), symIndex));
    }

    if (Fq_ab->calculationStatus(pc, symIndex, d.index()) == kUncalculated) {
        // measured breakup momentum
        double q2 = initialStateParticle()->measuredBreakupMomenta().q2(d.dataPoint(), pc);

        double R = DecayChannel_->parent()->radialSize()->value();
        double f = sqrt(F2(DecayChannel_->spinAmplitude()->twoL(), R*R, q2));
        Fq_ab->setValue(f, d.dataPoint(), symIndex, d.index());

        DEBUG("Blatt-Weisskopf barrier factor Fq_ab (L = " << spinToString(DecayChannel_->spinAmplitude()->twoL()) << ") = " << Fq_ab->value(d.dataPoint(), symIndex));
    }

    double Fq_rOFq_ab = Fq_r->value(d.dataPoint(), symIndex) / Fq_ab->value(d.dataPoint(), symIndex);
    //DEBUG("Blatt-Weisskopf barrier factor ratio (L = " << spinToString(DecayChannel_->spinAmplitude()->twoL()) << ") = " << Fq_rOFq_ab);
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
            return 9. + 3.*z + z*z;
        default:
            LOG(ERROR) << "calculation of Blatt-Weisskopf barrier factor is not (yet) implemented for L = "
                       << spinToString(twoL) << ". returning 0." << std::endl;
            return 0;
    }
}

}


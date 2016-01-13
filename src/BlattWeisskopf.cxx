#include "BlattWeisskopf.h"

#include "CalculationStatus.h"
#include "Constants.h"
#include "DecayChannel.h"
#include "DecayingParticle.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "QuantumNumbers.h"
#include "Resonance.h"
#include "SpinAmplitude.h"

#include <stdexcept>

namespace yap {

//-------------------------
double BlattWeisskopf::F2(unsigned l, double z)
{
    switch (l) {
        case 0:
            return 1.;
        case 1:
            return 1. + z;
        case 3:
            return 9. + 3.*z + z * z;
        default:
            /// \todo put in generic formula for L > 2
            throw exceptions::Exception("BlattWeisskopf does not yet support l > 2", "BlattWeisskopf::F2");
    }
}

//-------------------------
BlattWeisskopf::BlattWeisskopf(unsigned L, DecayingParticle* dp) :
    DataAccessor(&ParticleCombination::equivDownByOrderlessContent),
    DecayingParticle_(dp),
    L_(L),
    Fq_r(new RealCachedDataValue(this)),
    Fq_ab(new RealCachedDataValue(this))
{
    if (!DecayingParticle_)
        throw exceptions::Exception("DecayingParticle unset", "BlattWeisskopf::BlattWeisskopf");

    if (!initialStateParticle())
        throw exceptions::Exception("InitialStateParticle unset", "BlattWeisskopf::BlattWeisskopf");

    Fq_r->addDependency(initialStateParticle()->fourMomenta().masses());
    Fq_r->addDependency(DecayingParticle_->mass());
    Fq_r->addDependency(DecayingParticle_->radialSize());

    Fq_ab->addDependency(initialStateParticle()->measuredBreakupMomenta().breakupMomenta());
    Fq_ab->addDependency(DecayingParticle_->radialSize());
}

//-------------------------
std::complex<double> BlattWeisskopf::amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc, unsigned dataPartitionIndex) const
{
    unsigned symIndex = symmetrizationIndex(pc);
    bool calc(false); // for debugging

    if (Fq_r->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated) {
        // nominal breakup momentum
        double m2_R = pow(DecayingParticle_->mass()->value(), 2);
        double m_a = initialStateParticle()->fourMomenta().m(d, pc->daughters().at(0));
        double m_b = initialStateParticle()->fourMomenta().m(d, pc->daughters().at(1));
        double q2 = MeasuredBreakupMomenta::calcQ2(m2_R, m_a, m_b);

        DEBUG(*(pc->daughters()[0]) << " " << * (pc->daughters()[1]));
        //DEBUG(initialStateParticle()->fourMomenta().symmetrizationIndex(pc->daughters().at(0)) <<
        //        " " << initialStateParticle()->fourMomenta().symmetrizationIndex(pc->daughters().at(1)));
        DEBUG(m_a << " " << m_b);

        double R = DecayingParticle_->radialSize()->value();
        double f = sqrt(F2(L_, R * R * q2));
        Fq_r->setValue(f, d, symIndex, dataPartitionIndex);

        calc = true;
        DEBUG("BlattWeisskopf::amplitude - calculated barrier factor Fq_r (L = " << L_ << ") = " << Fq_r->value(d, symIndex));
    }

    if (Fq_ab->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated) {
        // measured breakup momentum
        double q2 = initialStateParticle()->measuredBreakupMomenta().q2(d, pc);

        double R = DecayingParticle_->radialSize()->value();
        double f = sqrt(F2(L_, R * R * q2));
        Fq_ab->setValue(f, d, symIndex, dataPartitionIndex);

        calc = true;
        DEBUG("BlattWeisskopf::amplitude - calculated barrier factor Fq_ab (L = " << L_ << ") = " << Fq_ab->value(d, symIndex));
    }

    double Fq_rOFq_ab = Fq_r->value(d, symIndex) / Fq_ab->value(d, symIndex);

    if (calc) {
        DEBUG("BlattWeisskopf::amplitude - using calculated values to calculate Blatt-Weisskopf barrier factor ratio (L = " << L_ << ") = " << Fq_rOFq_ab);
    } else {
        DEBUG("BlattWeisskopf::amplitude - using cached values to calculate Blatt-Weisskopf barrier factor ratio (L = " << L_ << ") = " << Fq_rOFq_ab);
    }

    return std::complex<double>(Fq_rOFq_ab);
}

//-------------------------
bool BlattWeisskopf::consistent() const
{
    bool C = true;

    if (!Fq_r->dependsOn(initialStateParticle()->fourMomenta().masses())) {
        FLOG(ERROR) << "Fq_r doesn't have mass dependencies set";
        C &= false;
    }
    if (!Fq_ab->dependsOn(initialStateParticle()->measuredBreakupMomenta().breakupMomenta())) {
        FLOG(ERROR) << "Fq_ab doesn't have breakup-momenta dependencies set";
        C &= false;
    }

    return C;
}

//-------------------------
InitialStateParticle* BlattWeisskopf::initialStateParticle()
{
    return DecayingParticle_->initialStateParticle();
}

}


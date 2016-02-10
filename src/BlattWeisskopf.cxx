#include "BlattWeisskopf.h"

#include "CalculationStatus.h"
#include "Constants.h"
#include "DecayChannel.h"
#include "DecayingParticle.h"
#include "logging.h"
#include "Model.h"
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
        case 2:
            return 9. + 3.*z + z * z;
        default:
            /// \todo put in generic formula for L > 2
            throw exceptions::Exception("BlattWeisskopf does not yet support L = " + std::to_string(l) + " > 2",
                                        "BlattWeisskopf::F2");
    }
}

//-------------------------
BlattWeisskopf::BlattWeisskopf(unsigned L, DecayingParticle* dp) :
    DataAccessor(&ParticleCombination::equivDownByOrderlessContent),
    DecayingParticle_(dp),
    L_(L),
    Fq_r(RealCachedDataValue::create(this)),
    Fq_ab(RealCachedDataValue::create(this))
{
    if (!DecayingParticle_)
        throw exceptions::Exception("DecayingParticle unset", "BlattWeisskopf::BlattWeisskopf");

    if (!model())
        throw exceptions::Exception("Model unset", "BlattWeisskopf::BlattWeisskopf");

    Fq_r->addDependency(model()->fourMomenta().mass());
    Fq_r->addDependency(DecayingParticle_->mass());
    Fq_r->addDependency(DecayingParticle_->radialSize());

    Fq_ab->addDependency(model()->measuredBreakupMomenta().breakupMomenta());
    Fq_ab->addDependency(DecayingParticle_->radialSize());

    // register with model
    addToModel();
}

//-------------------------
double BlattWeisskopf::amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc, unsigned dataPartitionIndex) const
{
    unsigned symIndex = symmetrizationIndex(pc);

    if (Fq_r->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated) {
        // nominal breakup momentum
        double m2_R = pow(DecayingParticle_->mass()->value(), 2);
        double m_a = model()->fourMomenta().m(d, pc->daughters().at(0));
        double m_b = model()->fourMomenta().m(d, pc->daughters().at(1));
        double q2 = MeasuredBreakupMomenta::calcQ2(m2_R, m_a, m_b);

        double R = DecayingParticle_->radialSize()->value();
        double f = sqrt(F2(L_, R * R * q2));
        Fq_r->setValue(f, d, symIndex, dataPartitionIndex);

        // DEBUG("BlattWeisskopf::amplitude - calculated barrier factor Fq_r (L = " << L_ << ") = " << Fq_r->value(d, symIndex));
    }

    if (Fq_ab->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated) {
        // measured breakup momentum
        double q2 = model()->measuredBreakupMomenta().q2(d, pc);

        double R = DecayingParticle_->radialSize()->value();
        double f = sqrt(F2(L_, R * R * q2));
        Fq_ab->setValue(f, d, symIndex, dataPartitionIndex);

        // DEBUG("BlattWeisskopf::amplitude - calculated barrier factor Fq_ab (L = " << L_ << ") = " << Fq_ab->value(d, symIndex));
    }

    return Fq_r->value(d, symIndex) / Fq_ab->value(d, symIndex);
}

//-------------------------
bool BlattWeisskopf::consistent() const
{
    bool C = true;

    if (!Fq_r->dependsOn(model()->fourMomenta().mass())) {
        FLOG(ERROR) << "Fq_r doesn't have mass dependencies set";
        C &= false;
    }
    if (!Fq_ab->dependsOn(model()->measuredBreakupMomenta().breakupMomenta())) {
        FLOG(ERROR) << "Fq_ab doesn't have breakup-momenta dependencies set";
        C &= false;
    }

    return C;
}

//-------------------------
Model* BlattWeisskopf::model()
{
    return DecayingParticle_->model();
}

}


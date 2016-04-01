#include "BlattWeisskopf.h"

#include "CalculationStatus.h"
#include "DecayingParticle.h"
#include "Exceptions.h"
#include "FourMomenta.h"
#include "logging.h"
#include "MeasuredBreakupMomenta.h"
#include "Model.h"

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
    AmplitudeComponent(),
    DataAccessor(&ParticleCombination::equivDownByOrderlessContent),
    RequiresMeasuredBreakupMomenta(L > 0),
    DecayingParticle_(dp),
    L_(L)
{
    if (!DecayingParticle_)
        throw exceptions::Exception("DecayingParticle unset", "BlattWeisskopf::BlattWeisskopf");

    if (!model())
        throw exceptions::Exception("Model unset", "BlattWeisskopf::BlattWeisskopf");

    if (L_ > 2)
        throw exceptions::Exception("BlattWeisskopf does not yet support L = " + std::to_string(L_),
                                    "BlattWeisskopf::BlattWeisskopf");

    if (L_ > 0) {
        // register with model
        addToModel();

        Fq_r = RealCachedDataValue::create(this);
        Fq_r->addDependency(DaughterCachedDataValue(model()->fourMomenta()->mass(), 0));
        Fq_r->addDependency(DaughterCachedDataValue(model()->fourMomenta()->mass(), 1));
        Fq_r->addDependency(DecayingParticle_->mass());
        Fq_r->addDependency(DecayingParticle_->radialSize());

        Fq_ab = RealCachedDataValue::create(this);
        Fq_ab->addDependency(model()->measuredBreakupMomenta()->breakupMomenta());
        Fq_ab->addDependency(DecayingParticle_->radialSize());
    }

    // if L == 0, values are all always 1, no storage in DataPoint necessary

}

//-------------------------
double BlattWeisskopf::amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc, StatusManager& sm) const
{
    // spin 0 always has factors 1
    if (L() == 0)
        return 1.;

    unsigned symIndex = symmetrizationIndex(pc);

    if (sm.status(*Fq_r, symIndex) == CalculationStatus::uncalculated) {
        // nominal breakup momentum
        double m2_R = pow(DecayingParticle_->mass()->value(), 2);
        double m_a = model()->fourMomenta()->m(d, pc->daughters().at(0));
        double m_b = model()->fourMomenta()->m(d, pc->daughters().at(1));
        double q2 = MeasuredBreakupMomenta::calcQ2(m2_R, m_a, m_b);

        double R = DecayingParticle_->radialSize()->value();
        double f = sqrt(F2(L_, R * R * q2));
        Fq_r->setValue(f, d, symIndex, sm);

        DEBUG("BlattWeisskopf::amplitude :: calculated barrier factor Fq_r (L = " << L_ << ") = " << Fq_r->value(d, symIndex));
    }

    if (sm.status(*Fq_ab, symIndex) == CalculationStatus::uncalculated) {
        // measured breakup momentum
        double q2 = model()->measuredBreakupMomenta()->q2(d, pc);

        double R = DecayingParticle_->radialSize()->value();
        double f = sqrt(F2(L_, R * R * q2));
        Fq_ab->setValue(f, d, symIndex, sm);

        DEBUG("BlattWeisskopf::amplitude :: calculated barrier factor Fq_ab (L = " << L_ << ") = " << Fq_ab->value(d, symIndex));
    }

    return Fq_r->value(d, symIndex) / Fq_ab->value(d, symIndex);
}

//-------------------------
const Model* BlattWeisskopf::model() const
{
    return DecayingParticle_->model();
}

}


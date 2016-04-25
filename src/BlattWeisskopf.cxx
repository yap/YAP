#include "BlattWeisskopf.h"

#include "CalculationStatus.h"
#include "DataPartition.h"
#include "DecayingParticle.h"
#include "Exceptions.h"
#include "FourMomenta.h"
#include "logging.h"
#include "MeasuredBreakupMomenta.h"
#include "Model.h"

namespace yap {

//-------------------------
double f_inverse_square(unsigned l, double z)
{
    switch (l) {
        case 0:
            return 1.;
        case 1:
            return 1. + z;
        case 2:
            return 9. + 3. * z + z * z;
        default:
            /// \todo put in generic formula for L > 2
            throw exceptions::Exception("BlattWeisskopf does not yet support L = " + std::to_string(l) + " > 2",
                                        "f_inverse_square");
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

        BarrierFactor_ = RealCachedDataValue::create(this);
        BarrierFactor_->addDependency(DaughterCachedDataValue(model()->fourMomenta()->mass(), 0));
        BarrierFactor_->addDependency(DaughterCachedDataValue(model()->fourMomenta()->mass(), 1));
        BarrierFactor_->addDependency(DecayingParticle_->mass());
        BarrierFactor_->addDependency(DecayingParticle_->radialSize());
        BarrierFactor_->addDependency(model()->measuredBreakupMomenta()->breakupMomenta());
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

    if (sm.status(*BarrierFactor_, symIndex) == CalculationStatus::uncalculated) {

        double m2_R = pow(DecayingParticle_->mass()->value(), 2);
        double m_a = model()->fourMomenta()->m(d, pc->daughters().at(0));
        double m_b = model()->fourMomenta()->m(d, pc->daughters().at(1));

        // nominal breakup momentum
        double q2_nomi = MeasuredBreakupMomenta::calcQ2(m2_R, m_a, m_b);

        // measured breakup momentum
        double q2_meas = model()->measuredBreakupMomenta()->q2(d, pc);

        double r2 = pow(DecayingParticle_->radialSize()->value(), 2);
        double f2_nomi = f_inverse_square(L_, r2 * q2_nomi);
        double f2_meas = f_inverse_square(L_, r2 * q2_meas);

        double barrier_factor = sqrt(f2_nomi / f2_meas);

        BarrierFactor_->setValue(barrier_factor, d, symIndex, sm);

        return barrier_factor;
    }

    return BarrierFactor_->value(d, symIndex);
}

//-------------------------
double BlattWeisskopf::operator()(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
{
    return (L_ == 0) ? 1 : BarrierFactor_->value(d, symmetrizationIndex(pc));
}

//-------------------------
void BlattWeisskopf::calculate(DataPartition& D, const std::shared_ptr<ParticleCombination>& pc) const
{
    // spin 0 always has factors 1, no calculations necessary
    if (L_ == 0)
        return;

    unsigned symIndex = symmetrizationIndex(pc);

    if (D.status(*BarrierFactor_, symIndex) == CalculationStatus::uncalculated) {

        for (auto& d : D) {

            double m2_R = pow(DecayingParticle_->mass()->value(), 2);
            double m_a = model()->fourMomenta()->m(d, pc->daughters().at(0));
            double m_b = model()->fourMomenta()->m(d, pc->daughters().at(1));

            // nominal breakup momentum
            double q2_nomi = MeasuredBreakupMomenta::calcQ2(m2_R, m_a, m_b);

            // measured breakup momentum
            double q2_meas = model()->measuredBreakupMomenta()->q2(d, pc);

            double r2 = pow(DecayingParticle_->radialSize()->value(), 2);
            double f2_nomi = f_inverse_square(L_, r2 * q2_nomi);
            double f2_meas = f_inverse_square(L_, r2 * q2_meas);

            double barrier_factor = sqrt(f2_nomi / f2_meas);

            BarrierFactor_->setValue(barrier_factor, d, symIndex, D);
        }

        D.status(*BarrierFactor_, symIndex) = CalculationStatus::calculated;
    }
}

//-------------------------
CachedDataValueSet BlattWeisskopf::cachedDataValuesItDependsOn()
{
    /// \todo replace with actual dependencies of BarrierFactor_?
    CachedDataValueSet set;
    if (BarrierFactor_)
        set.insert(BarrierFactor_);
    return set;
}

//-------------------------
const Model* BlattWeisskopf::model() const
{
    return DecayingParticle_->model();
}

}


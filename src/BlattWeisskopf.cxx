#include "BlattWeisskopf.h"

#include "CalculationStatus.h"
#include "DataPartition.h"
#include "DecayingParticle.h"
#include "Exceptions.h"
#include "FourMomenta.h"
#include "logging.h"
#include "MeasuredBreakupMomenta.h"
#include "Model.h"
#include "Parameter.h"

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
    RecalculableDataAccessor(ParticleCombination::equivDownByOrderlessContent),
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

        addParameter(DecayingParticle_->mass());
        addParameter(DecayingParticle_->radialSize());

        BarrierFactor_ = RealCachedDataValue::create(this);
    }

    // if L == 0, values are all always 1, no storage in DataPoint necessary

}

//-------------------------
std::complex<double> BlattWeisskopf::value(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
{
    return (L_ == 0) ? 1 : BarrierFactor_->value(d, symmetrizationIndex(pc));
}

//-------------------------
void BlattWeisskopf::calculate(DataPartition& D) const
{
    // spin 0 always has factors 1, no calculations necessary
    if (L_ == 0)
        return;

    // loop over (ParticleCombination --> symmetrization index) map
    for (const auto& pc_symIndex : symmetrizationIndices()) {

        // check if barrier factor is uncalculated
        if (D.status(*BarrierFactor_, pc_symIndex.second) == CalculationStatus::uncalculated) {

            DEBUG("calculate BlattWeisskopf");

            // calculate on all data points in D
            for (auto& d : D) {

                double m2_R = pow(DecayingParticle_->mass()->value(), 2);
                double m_a = model()->fourMomenta()->m(d, pc_symIndex.first->daughters().at(0));
                double m_b = model()->fourMomenta()->m(d, pc_symIndex.first->daughters().at(1));

                // nominal breakup momentum
                double q2_nomi = MeasuredBreakupMomenta::calcQ2(m2_R, m_a, m_b);

                // measured breakup momentum
                double q2_meas = model()->measuredBreakupMomenta()->q2(d, pc_symIndex.first);

                double r2 = pow(DecayingParticle_->radialSize()->value(), 2);
                double f2_nomi = f_inverse_square(L_, r2 * q2_nomi);
                double f2_meas = f_inverse_square(L_, r2 * q2_meas);

                double barrier_factor = sqrt(f2_nomi / f2_meas);

                // store result in data point
                BarrierFactor_->setValue(barrier_factor, d, pc_symIndex.second, D);
            }

            // update status
            D.status(*BarrierFactor_, pc_symIndex.second) = CalculationStatus::calculated;
        }
    }
}

//-------------------------
void BlattWeisskopf::updateCalculationStatus(StatusManager& D) const
{
    if (variableStatus(*this) == VariableStatus::changed)
        D.set(*BarrierFactor_, CalculationStatus::uncalculated);
}

//-------------------------
const Model* BlattWeisskopf::model() const
{
    return DecayingParticle_->model();
}

}


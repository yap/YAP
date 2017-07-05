#include "NonrelativisticBreitWigner.h"

#include "BlattWeisskopf.h"
#include "CachedValue.h"
#include "CalculationStatus.h"
#include "DataPartition.h"
#include "DecayChannel.h"
#include "DecayingParticle.h"
#include "FourMomenta.h"
#include "logging.h"
#include "MathUtilities.h"
#include "MeasuredBreakupMomenta.h"
#include "Model.h"
#include "Parameter.h"

namespace yap {

//-------------------------
void NonrelativisticBreitWigner::calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const
{
    // if no calculation necessary, exit
    if (D.status(*T(), si) != CalculationStatus::uncalculated)
        return;

    /////////////////////////
    // common factors:

    // radial size squared
    double r2 = pow(owner()->radialSize()->value(), 2);

    // squared resonance mass
    double m2_R = pow(mass()->value(), 2);

    // i * nominal mass * nominal width
    auto mw_R = mass()->value() * width()->value();

    // J + 1/2
    unsigned twoLp1 = blattWeisskopf()->L() + 1;

    /////////////////////////

    // calculate BlattWeisskopf
    blattWeisskopf()->calculate(D);

    // T := 1 / (M^2 - m^2 - i * M * Gamma)
    for (auto& d : D) {

        // retrieve masses
        double m_ab = model()->fourMomenta()->m(d, pc);
        double m2_ab = m_ab * m_ab;
        double m_a  = model()->fourMomenta()->m(d, pc->daughters()[0]);
        double m_b  = model()->fourMomenta()->m(d, pc->daughters()[1]);

        // calculate nominal breakup momentum
        double q2_nomi = measured_breakup_momenta::q2(m2_R, m_a, m_b);

        // retrieve measured breakup momentum
        double q2_meas = measured_breakup_momenta::q2(m2_ab, m_a, m_b);

        double Q = sqrt(q2_meas / q2_nomi);
        
        auto w =  0.5 * mw_R / m_ab * pow(Q, twoLp1) * pow(blattWeisskopf()->value(d, pc), 2) / squared_barrier_factor(blattWeisskopf()->L(), q2_nomi * r2);

        T()->setValue(w / (mass()->value() - m_ab - 1_i * w), d, si, D);

    }

    D.status(*T(), si) = CalculationStatus::calculated;
}

}





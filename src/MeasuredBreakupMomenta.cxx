#include "MeasuredBreakupMomenta.h"

#include "FourMomenta.h"
#include "InitialStateParticle.h"
#include "logging.h"

#include <TLorentzVector.h>

#include <set>

namespace yap {

//-------------------------
MeasuredBreakupMomenta::MeasuredBreakupMomenta() :
    DataAccessor(&ParticleCombination::equivDownByOrderlessContent),
    Q2_(this)
{
}

//-------------------------
void MeasuredBreakupMomenta::calculate(DataPoint& d)
{
    Q2_.setCalculationStatus(kUncalculated);

    for (auto& kv : SymmetrizationIndices_) {

        // check if calculation necessary
        if (Q2_.calculationStatus(kv) == kCalculated)
            continue;

        if (kv.first->daughters().size() != 2) {
            LOG(ERROR) << "MeasuredBreakupMomenta::calculate - invalid number of daughters (" << kv.first->daughters().size() << " != 2)";
            return;
        }

        // Calculate
        double m2_R = initialStateParticle()->fourMomenta().m2(d, kv.first);
        double m_a  = initialStateParticle()->fourMomenta().m(d, kv.first->daughters()[0]);
        double m_b  = initialStateParticle()->fourMomenta().m(d, kv.first->daughters()[1]);

        double q2 = (m2_R - pow(m_a + m_b, 2)) * (m2_R - pow(m_a - m_b, 2)) / m2_R / 4.;

        Q2_.setValue(q2, d, kv.second);
    }
}


}

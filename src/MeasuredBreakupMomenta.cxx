#include "MeasuredBreakupMomenta.h"

#include "FourMomenta.h"
#include "InitialStateParticle.h"
#include "logging.h"

#include <set>

namespace yap {

//-------------------------
MeasuredBreakupMomenta::MeasuredBreakupMomenta() :
    DataAccessor(&ParticleCombination::equivDownByOrderlessContent)
{
}

//-------------------------
void MeasuredBreakupMomenta::calculate(DataPoint& d)
{
    std::set<unsigned> alreadyCalculated;

    for (auto& kv : SymmetrizationIndices_) {

        if (alreadyCalculated.find(kv.second) != alreadyCalculated.end())
            continue;

        if (kv.first->daughters().size() != 2) {
            LOG(ERROR) << "MeasuredBreakupMomenta::calculate - invalid number of daughters (" << kv.first->daughters().size() << " != 2)";
            return;
        }

        double m2_R = initialStateParticle()->fourMomenta().m2(d, kv.first);
        double m_a  = initialStateParticle()->fourMomenta().m(d, kv.first->daughters()[0]);
        double m_b  = initialStateParticle()->fourMomenta().m(d, kv.first->daughters()[1]);

        d.MeasuredBreakupMomenta_.at(kv.second) = (m2_R - pow(m_a + m_b, 2)) * (m2_R - pow(m_a - m_b, 2)) / m2_R / 4.0;
        alreadyCalculated.insert(kv.second);

        //DEBUG("breakup momentum for " << std::string(*kv.first) << " = " << d.MeasuredBreakupMomenta_.at(symmetrizationIndex(kv.first)));
    }
}


}

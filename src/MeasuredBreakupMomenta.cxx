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

    for (auto& pc : particleCombinations()) {
        unsigned index = symmetrizationIndex(pc);

        if (alreadyCalculated.find(index) != alreadyCalculated.end())
            continue;

        if (pc->daughters().size() != 2) {
            LOG(ERROR) << "MeasuredBreakupMomenta::calculate - invalid number of daughters (" << pc->daughters().size() << " != 2)";
            return;
        }

        double m2_R = initialStateParticle()->fourMomenta().m2(d, pc);
        double m_a = initialStateParticle()->fourMomenta().m(d, pc->daughters()[0]);
        double m_b = initialStateParticle()->fourMomenta().m(d, pc->daughters()[1]);

        d.MeasuredBreakupMomenta_.at(symmetrizationIndex(pc)) = (m2_R - (m_a + m_b) * (m_a + m_b)) * (m2_R - (m_a - m_b) * (m_a - m_b)) / m2_R / 4.0;
        alreadyCalculated.insert(index);

        //DEBUG("breakup momentum for " << std::string(*pc) << " = " << d.MeasuredBreakupMomenta_.at(symmetrizationIndex(pc)));
    }
}


}

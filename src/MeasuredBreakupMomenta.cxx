#include "MeasuredBreakupMomenta.h"

#include "FourMomenta.h"
#include "InitialStateParticle.h"
#include "logging.h"

#include <TLorentzVector.h>

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

        // Calculate
        double m2_R = initialStateParticle()->fourMomenta().m2(d, kv.first);
        double m_a  = initialStateParticle()->fourMomenta().m(d, kv.first->daughters()[0]);
        double m_b  = initialStateParticle()->fourMomenta().m(d, kv.first->daughters()[1]);

        d.MeasuredBreakupMomenta_.at(kv.second) = sqrt((m2_R - pow(m_a + m_b, 2)) * (m2_R - pow(m_a - m_b, 2)) / m2_R / 4.);
        alreadyCalculated.insert(kv.second);

        //DEBUG("breakup momentum for " << std::string(*kv.first) << " = " << d.MeasuredBreakupMomenta_.at(kv.second));

        // alternative way gives same results
        /*TLorentzVector R(0,0,0,0);
        for (ParticleIndex i : kv.first->indices())
            R += initialStateParticle()->fourMomenta().p(d, i);

        TLorentzVector A(0,0,0,0);
        for (ParticleIndex i : kv.first->daughters()[0]->indices())
            A += initialStateParticle()->fourMomenta().p(d, i);

        // boost
        A.Boost(-R.BoostVector());

        d.MeasuredBreakupMomenta_.at(kv.second) = A.Vect().Mag();
        alreadyCalculated.insert(kv.second);

        DEBUG("4vec breakup mom " << A.Vect().Mag());*/
    }
}


}

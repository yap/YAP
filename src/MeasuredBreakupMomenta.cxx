#include "MeasuredBreakupMomenta.h"

#include "Exceptions.h"
#include "FourMomenta.h"
#include "Model.h"

namespace yap {

//-------------------------
MeasuredBreakupMomenta::MeasuredBreakupMomenta(Model* m) :
    StaticDataAccessor(m, &ParticleCombination::equivDownByOrderlessContent),
    Q2_(RealCachedDataValue::create(this))
{
    if (!model()->fourMomenta())
        throw exceptions::Exception("Model's FourMomenta unset", "MeasuredBreakupMomenta::MeasuredBreakupMomenta");
    Q2_->addDependency(model()->fourMomenta()->mass());
    Q2_->addDependency(std::make_pair(model()->fourMomenta()->mass(), 0));
    Q2_->addDependency(std::make_pair(model()->fourMomenta()->mass(), 1));
}

//-------------------------
void MeasuredBreakupMomenta::calculate(DataPoint& d, unsigned dataPartitionIndex)
{
    Q2_->setCalculationStatus(kUncalculated, dataPartitionIndex);

    for (auto& kv : symmetrizationIndices()) {

        // check if calculation necessary
        if (Q2_->calculationStatus(kv.first, kv.second, dataPartitionIndex) == kCalculated)
            continue;

        if (kv.first->daughters().size() != 2)
            throw exceptions::Exception("invalid number of daughters (" + std::to_string(kv.first->daughters().size()) + ")",
                                        "MeasuredBreakupMomenta::calculate");

        // Calculate
        double m2_R = model()->fourMomenta()->m2(d, kv.first);
        double m_a  = model()->fourMomenta()->m(d, kv.first->daughters()[0]);
        double m_b  = model()->fourMomenta()->m(d, kv.first->daughters()[1]);

        Q2_->setValue(calcQ2(m2_R, m_a, m_b), d, kv.second, dataPartitionIndex);
    }
}

//-------------------------
double MeasuredBreakupMomenta::calcQ2(double m2_R, double m_a, double m_b)
{
    if (m_a == m_b)
        return m2_R / 4. - m_a * m_a;

    return (m2_R - pow(m_a + m_b, 2)) * (m2_R - pow(m_a - m_b, 2)) / m2_R / 4.;
}

//-------------------------
unsigned MeasuredBreakupMomenta::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    if (pc->isFinalStateParticle())
        throw exceptions::FinalStateParticleCombination("cannot calculate helicity angles for fsp", "MeasuredBreakupMomenta::addParticleCombination");
    return StaticDataAccessor::addParticleCombination(pc);
}


}

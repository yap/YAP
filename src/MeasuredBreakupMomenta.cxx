#include "MeasuredBreakupMomenta.h"

#include "CachedValue.h"
#include "CalculationStatus.h"
#include "Exceptions.h"
#include "FourMomenta.h"
#include "Model.h"
#include "ParticleCombination.h"
#include "StatusManager.h"

namespace yap {

//-------------------------
MeasuredBreakupMomenta::MeasuredBreakupMomenta(Model& m) :
    StaticDataAccessor(m, equal_down_by_orderless_content),
    Q2_(RealCachedValue::create(*this))
{
    registerWithModel();
}

//-------------------------
void MeasuredBreakupMomenta::addToStaticDataAccessors()
{
    // look for Model's FourMomenta_
    auto it_fm = std::find(staticDataAccessors().begin(), staticDataAccessors().end(), model()->fourMomenta().get());
    if (it_fm == staticDataAccessors().end())
        throw exceptions::Exception("MeasuredBreakupMomenta cannot be registered with the model before FourMomenta", "HelicityAngles::registerWithModel");
    // add this to just after FourMomenta_
    const_cast<StaticDataAccessorVector&>(staticDataAccessors()).insert(it_fm + 1, this);
}

//-------------------------
double MeasuredBreakupMomenta::q2(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
{
    return Q2_->value(d, symmetrizationIndex(pc));
}


//-------------------------
void MeasuredBreakupMomenta::calculate(DataPoint& d, StatusManager& sm) const
{
    // set Q2 uncalculated
    sm.set(*Q2_, CalculationStatus::uncalculated);

    for (auto& kv : symmetrizationIndices()) {

        // check if calculation unnecessary
        if (sm.status(*Q2_, kv.second) == CalculationStatus::calculated)
            continue;

        if (kv.first->daughters().size() != 2)
            throw exceptions::Exception("invalid number of daughters (" + std::to_string(kv.first->daughters().size()) + ")",
                                        "MeasuredBreakupMomenta::calculate");

        // Calculate
        double m2_R = model()->fourMomenta()->m2(d, kv.first);
        double m_a  = model()->fourMomenta()->m(d, kv.first->daughters()[0]);
        double m_b  = model()->fourMomenta()->m(d, kv.first->daughters()[1]);

        Q2_->setValue(squared_breakup_momentum(m2_R, m_a, m_b), d, kv.second, sm);
    }
}

//-------------------------
void MeasuredBreakupMomenta::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    if (pc->daughters().size() != 2)
        throw exceptions::NotTwoBodyParticleCombination("cannot calculate breakup momentum for "
                + std::to_string(pc->daughters().size()) + "-body decay",
                "MeasuredBreakupMomenta::addParticleCombination");
    StaticDataAccessor::addParticleCombination(pc);
}


}

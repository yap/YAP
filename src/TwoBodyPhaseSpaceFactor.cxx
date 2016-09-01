#include "TwoBodyPhaseSpaceFactor.h"

#include "DecayChannel.h"
#include "FourMomenta.h"
#include "MeasuredBreakupMomenta.h"
#include "Model.h"
#include "ParticleCombination.h"

namespace yap {

//-------------------------
TwoBodyPhaseSpaceFactor::TwoBodyPhaseSpaceFactor(const Model& m)
    : PhaseSpaceFactor(equal_down_by_orderless_content),
      RequiresMeasuredBreakupMomenta(true),
      Model_(&m)
{
}

//-------------------------
const std::complex<double> TwoBodyPhaseSpaceFactor::value(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
{
    return pow(4. * model()->measuredBreakupMomenta()->q2(d, pc) / model()->fourMomenta()->m2(d, pc), 0.25);
    // return sqrt(model()->measuredBreakupMomenta()->rho(d, pc));
}

//-------------------------
void TwoBodyPhaseSpaceFactor::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    if (pc->daughters().size() != 2)
        throw exceptions::NotTwoBodyParticleCombination("Not two-body decay", "PhaseSpaceFactor::addParticleCombination");
    PhaseSpaceFactor::addParticleCombination(pc);
}

//-------------------------
std::shared_ptr<PhaseSpaceFactor> TwoBodyPhaseSpaceFactorFactory::phaseSpaceFactor(const DecayChannel& dc, const SpinAmplitude& sa, std::shared_ptr<MassShape> ms)
{
    // for 3-or-more-body decay, phsp factor = 1 --> return nullptr
    if (dc.daughters().size() > 2)
        return nullptr;
    
    if (!dc.model())
        throw exceptions::Exception("Model is nullptr", "TwoBodyPhaseSpaceFactorFactory::phaseSpaceFactor");

    // search for existing phase-space factor for this model
    auto it = PHSPFactors_.find(dc.model());
    // if found, return it
    if (it != PHSPFactors_.end())
        return it->second;

    // else make it
    PHSPFactors_[dc.model()] = std::make_shared<TwoBodyPhaseSpaceFactor>(*dc.model());
    return PHSPFactors_[dc.model()];
}

}

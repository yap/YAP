#include "../include/ResonantDecayChannel.h"

#include "BlattWeisskopf.h"
#include "container_utils.h"
#include "CachedDataValue.h"
#include "CalculationStatus.h"
#include "DecayingParticle.h"
#include "Exceptions.h"
#include "FinalStateParticle.h"
#include "FreeAmplitude.h"
#include "logging.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleCombinationCache.h"
#include "spin.h"
#include "SpinAmplitude.h"
#include "SpinAmplitudeCache.h"
#include "StatusManager.h"

namespace yap {

//-------------------------
ResonantDecayChannel::ResonantDecayChannel(const ParticleVector& daughters) :
    DecayChannel(daughters)
{
    if (Daughters_.size() > 2)
        throw exceptions::Exception("More than two daughters", "ResonantDecayChannel::DecayChannel");
}

//-------------------------
void ResonantDecayChannel::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    DecayChannel::addParticleCombination(pc);

    // add to SpinAmplitude's
    for (auto& sa : SpinAmplitudes_)
        sa->addParticleCombination(pc);
}

//-------------------------
void ResonantDecayChannel::setDecayingParticle(DecayingParticle* dp)
{
    DecayChannel::setDecayingParticle(dp);

    // if SpinAmplitude's have already been added by hand, don't add automatically
    if (SpinAmplitudes_.empty()) {

        auto two_J = DecayingParticle_->quantumNumbers().twoJ();
        auto two_j1 = Daughters_[0]->quantumNumbers().twoJ();
        auto two_j2 = Daughters_[1]->quantumNumbers().twoJ();

        // create spin amplitudes
        // loop over possible S: |j1-j2| <= S <= (j1+j2)
        for (unsigned two_S = std::abs<int>(two_j1 - two_j2); two_S <= two_j1 + two_j2; two_S += 2) {
            // loop over possible L: |J-s| <= L <= (J+s)
            for (unsigned L = std::abs<int>(two_J - two_S) / 2; L <= (two_J + two_S) / 2; ++L) {
                // add SpinAmplitude retrieved from cache
                addSpinAmplitude(const_cast<Model*>(static_cast<const DecayChannel*>(this)->model())->spinAmplitudeCache()->spinAmplitude(two_J, two_j1, two_j2, L, two_S));
            }
        }
    } else {

        // check DecayingPartcle_'s quantum numbers against existing SpinAmplitude's
        if (DecayingParticle_->quantumNumbers().twoJ() != SpinAmplitudes_[0]->initialTwoJ())
            throw exceptions::Exception("Spins don't match ", "DecayChannel::setDecayingParticle");

    }

    // let DecayingParticle know to create a BlattWeisskopf objects for necessary orbital angular momenta
    for (auto& sa : SpinAmplitudes_)
        DecayingParticle_->storeBlattWeisskopf(sa->L());
}

//-------------------------
void ResonantDecayChannel::addSpinAmplitude(std::shared_ptr<SpinAmplitude> sa)
{
    // check number of daughters
    if (sa->finalTwoJ().size() != Daughters_.size())
        throw exceptions::Exception("Number of daughters doesn't match", "ResonantDecayChannel::addSpinAmplitude");

    /// \todo quantum numbers more completely?
    // check against daughter quantum numbers
    for (size_t i = 0; i < Daughters_.size(); ++i)
        if (Daughters_[i]->quantumNumbers().twoJ() != sa->finalTwoJ()[i])
            throw exceptions::Exception("Spins don't match daughter's", "ResonantDecayChannel::addSpinAmplitude");

    // check against DecayingParticle_ if set
    if (DecayingParticle_) {
        if (DecayingParticle_->quantumNumbers().twoJ() != sa->initialTwoJ())
            throw exceptions::Exception("Spins don't match DecayingParticle", "ResonantDecayChannel::addSpinAmplitude");
    } else {
        // else check against previously added SpinAmplitude's initial quantum numbers
        if (!SpinAmplitudes_.empty() and SpinAmplitudes_[0]->initialTwoJ() != sa->initialTwoJ())
            throw exceptions::Exception("Spins don't match previously added", "ResonantDecayChannel::addSpinAmplitude");
    }

    // add to SpinAmplitudes_
    SpinAmplitudes_.push_back(sa);

    // add this' ParticleCombination's to it
    for (auto& pc : particleCombinations())
        SpinAmplitudes_.back()->addParticleCombination(pc);

    sa->addToModel();
}

//-------------------------
bool ResonantDecayChannel::consistent() const
{
    bool C = DecayChannel::consistent();

    // loop over SpinAmplitude's
    for (const auto& sa : SpinAmplitudes_) {
        // check SpinAmplitude
        if (!sa) {
            FLOG(ERROR) << "A SpinAmplitude is empty";
            C &= false;
        } else
            C &= sa->consistent();
    }

    return C;
}

}

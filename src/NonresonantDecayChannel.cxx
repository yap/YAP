#include "../include/NonresonantDecayChannel.h"

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
NonresonantDecayChannel::NonresonantDecayChannel(const ParticleVector& daughters) :
    DecayChannel(daughters)
{
    if (Daughters_.size() < 3)
        throw exceptions::Exception("Less than three daughters", "NonresonantDecayChannel::NonresonantDecayChannel");

    for (auto& d : Daughters_)
        if (d->quantumNumbers().twoJ() != 0)
            throw exceptions::Exception("Spin != 0", "NonresonantDecayChannel::setDecayingParticle");
}

//-------------------------
void NonresonantDecayChannel::setDecayingParticle(DecayingParticle* dp)
{
    DecayChannel::setDecayingParticle(dp);

    // make sure that there is no relative angular momentum

    if (DecayingParticle_->quantumNumbers().twoJ() != 0)
        throw exceptions::Exception("Spin != 0", "NonresonantDecayChannel::setDecayingParticle");

    for (auto& d : Daughters_)
        if (d->quantumNumbers().twoJ() != 0)
            throw exceptions::Exception("Spin != 0", "NonresonantDecayChannel::setDecayingParticle");

}

}

#include "DecayChannel.h"

#include "DecayingParticle.h"
#include "FinalStateParticle.h"
#include "Resonance.h"
#include "logging.h"

namespace yap {

//-------------------------
DecayChannel::DecayChannel(Particle* daughterA, Particle* daughterB,
                           std::shared_ptr<SpinAmplitude> spinAmplitude) :
    Daughters_( {daughterA, daughterB}),
            BlattWeisskopf_(this),
            SpinAmplitude_(spinAmplitude),
            FreeAmplitude_(0)
{
  // set symmetrization indices
  std::vector<std::shared_ptr<ParticleCombination> > PCsA;
  if (dynamic_cast<DataAccessor*>(Daughters_[0]))
    PCsA = dynamic_cast<DataAccessor*>(Daughters_[0])->particleCombinations();
  else if (dynamic_cast<FinalStateParticle*>(Daughters_[0]))
    PCsA = dynamic_cast<FinalStateParticle*>(Daughters_[0])->particleCombinations();
  else
    LOG(ERROR) << "DecayChannel() - cannot get ParticleCombinations from daughterA";

  std::vector<std::shared_ptr<ParticleCombination> > PCsB;
  if (dynamic_cast<DataAccessor*>(Daughters_[1]))
    PCsB = dynamic_cast<DataAccessor*>(Daughters_[1])->particleCombinations();
  else if (dynamic_cast<FinalStateParticle*>(Daughters_[1]))
    PCsB = dynamic_cast<FinalStateParticle*>(Daughters_[1])->particleCombinations();
  else
    LOG(ERROR) << "DecayChannel() - cannot get ParticleCombinations from daughterB";

  for (std::shared_ptr<ParticleCombination> PCA : PCsA)
    for (std::shared_ptr<ParticleCombination> PCB : PCsB) {
      if (! PCA->shareIndices(PCB)) {
        std::shared_ptr<yap::ParticleCombination> a_b = yap::ParticleCombination::uniqueSharedPtr({PCA, PCB});
        // \todo what if daughterA == daughterB

        this->addSymmetrizationIndex(a_b);
        BlattWeisskopf_.addSymmetrizationIndex(a_b);
        SpinAmplitude_->addSymmetrizationIndex(a_b);
      }
    }
}

//-------------------------
Amp DecayChannel::amplitude(DataPoint& d)
{
    return Amp(1);
}

//-------------------------
bool DecayChannel::consistent() const
{
    bool consistent = true;

    consistent &= DataAccessor::consistent();

    // check daughters
    for (Daughters::const_iterator d = Daughters_.begin(); d != Daughters_.end(); ++d)
        consistent &= (*d)->consistent();
    if (! consistent) {
        LOG(ERROR) << "DecayChannel::consistent() - daughter inconsistent";
    }

    consistent &= BlattWeisskopf_.consistent();
    consistent &= spinAmplitude()->consistent();

    // check size of spin amplitude quantum numbers and size of daughters
    if (spinAmplitude()->finalQuantumNumbers().size() != Daughters_.size()) {
        LOG(ERROR) << "DecayChannel::consistent() - quantum numbers object and daughters object size mismatch";
        consistent = false;
    }

    // check if QuantumNumbers of SpinAmplitude objects match with Particles
    if (spinAmplitude()->initialQuantumNumbers() != parent()->quantumNumbers()) {
        LOG(ERROR) << "DecayChannel::consistent() - quantum numbers of parent and SpinResonance don't match.";
        consistent = false;
    }

    for (unsigned i = 0; i < Daughters_.size(); ++i) {
        if (spinAmplitude()->finalQuantumNumbers()[i] != Daughters_[i]->quantumNumbers()) {
            LOG(ERROR) << "DecayChannel::consistent() - quantum numbers of " << i << " and SpinResonance don't match.";
            consistent = false;
        }
    }

    // check if BlattWeisskopf points back to this DecayChannel
    if (this != BlattWeisskopf_.decayChannel()) {
        LOG(ERROR) << "DecayChannel::consistent() - BlattWeisskopf does not point back to this DecayChannel.";
        consistent =  false;
    }

    // check if INITIAL QuantumNumbers of SpinAmplitude objects match with this Resonance's QuantumNumbers
    if (spinAmplitude()->initialQuantumNumbers() != parent()->quantumNumbers()) {
        LOG(ERROR) << "DecayChannel::consistent() - quantum numbers of SpinAmplitude  and this channel's resonance don't match.";
        consistent =  false;
    }

    // check masses
    double finalMass = 0;
    for (Daughters::const_iterator d = Daughters_.begin(); d != Daughters_.end(); ++d)
        finalMass += (*d)->mass();
    if (finalMass > parent()->mass()) {
        LOG(ERROR) << "DecayChannel::consistent() - sum of daughter's masses is bigger than resonance mass.";
        consistent =  false;
    }

    if (! consistent) {
        LOG(ERROR) << "Channel is not consistent:  " << this->parent()->name() << " -> " << this->daughters()[0]->name() << " " << this->daughters()[1]->name() << "\n";
    }
    return consistent;
}

}

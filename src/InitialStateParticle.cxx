#include "InitialStateParticle.h"

#include "logging.h"

#include <TLorentzRotation.h>

#include <assert.h>

namespace yap {

//-------------------------
InitialStateParticle::InitialStateParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize) :
    DecayingParticle(this, q, mass, name, radialSize),
    FourMomenta_(this),
    HelicityAngles_(this)
{
}

//-------------------------
bool InitialStateParticle::consistent() const
{
    return DecayingParticle::consistent();
}

//-------------------------
void InitialStateParticle::setSymmetrizationIndexParents()
{
    std::cout << "InitialStateParticle::setSymmetrizationIndexParents()\n";

    ParticleCombination::makeParticleCombinationSetWithParents();

    ParticleCombination::printParticleCombinationSet();

    unsigned size = particleCombinations()[0]->indices().size();
    assert(size > 1);

    clearSymmetrizationIndices();

    for (auto& pc : ParticleCombination::particleCombinationSet()) {
        if (pc->indices().size() < size)
            continue;

        if (pc->daughters()[0]->parent() == nullptr)
            continue;

        std::cout << "  add " << std::string(*pc) << "\n";
        addSymmetrizationIndex(pc);
    }


    for (auto& ch : channels()) {
        std::vector<std::shared_ptr<ParticleCombination> > chPCs = ch->particleCombinations();
        ch->clearSymmetrizationIndices();

        // loop over channel's particle combinations
        for (auto& chPC : chPCs) {
            for (auto& pc : ParticleCombination::particleCombinationSet()) {
                if (ParticleCombination::equivDown(chPC, pc)) {

                    //std::cout << std::string(*chPC) << " == " << std::string(*pc) << "\n";
                    std::cout << "  add " << std::string(*pc);
                    if  (pc->parent()){
                      std::cout << " from decay " << std::string(*pc->parent());
                    }
                    std::cout << " to channel " << std::string(*ch) << "\n";

                    ch->addSymmetrizationIndex(pc);

                    // set PCs for channel's daughters
                    for (auto& pcDaughPC : pc->daughters()) {
                        for (const std::shared_ptr<Particle>& chDaugh : ch->daughters()) {
                            if (std::dynamic_pointer_cast<DecayingParticle>(chDaugh))
                                for (auto& chDaughPC : std::dynamic_pointer_cast<DecayingParticle>(chDaugh)->particleCombinations()) {
                                    if (ParticleCombination::equivDown(pcDaughPC, chDaughPC)) {
                                        std::dynamic_pointer_cast<DecayingParticle>(chDaugh)->addSymmetrizationIndex(pcDaughPC);
                                        std::cout << "   add " << std::string(*pcDaughPC);
                                        if (pcDaughPC->parent()) {
                                          std::cout << " from decay " << std::string(*pcDaughPC->parent());
                                        }
                                        std::cout << " to particle " << chDaugh->name() << "\n";
                                    }
                                }
                        }
                    }
                }
            }
        }

        std::cout << "  PCs in channel " << std::string(*ch) << ":\n";
        for (auto& meh : ch->particleCombinations())
            std::cout << "    " << std::string(*meh) << "\n";

        for (auto d : ch->daughters()) {
          if (std::dynamic_pointer_cast<DecayingParticle>(d)) {
            std::cout << "    PCs in particle " << d->name() << ":\n";
            for (auto& meh : std::dynamic_pointer_cast<DecayingParticle>(d)->particleCombinations()) {
                std::cout << "      " << std::string(*meh);
                if (meh->parent())
                  std::cout << " from decay " << std::string(*meh->parent());
                std::cout <<"\n";
            }
          }
        }
    }

    // next level
    for (auto& ch : channels()) {
        for (std::shared_ptr<Particle> d : ch->daughters()) {
            std::cout << "Now setSymmetrizationIndexParents for " << d->name() << "\n";
            d->setSymmetrizationIndexParents();
        }
    }
}

//-------------------------
/*void InitialStateParticle::addSymmetrizationIndex(std::shared_ptr<ParticleCombination> c)
{
  DecayingParticle::addSymmetrizationIndex(c);
  FourMomenta_.addSymmetrizationIndex(c);
  HelicityAngles_.addSymmetrizationIndex(c);
}*/

//-------------------------
bool InitialStateParticle::addDataPoint(DataPoint&& d)
{
    FourMomenta_.calculate(d);
    HelicityAngles_.calculate(d);
    if (!DataSet_.consistent(d))
        return false;
    return DataSet_.addDataPoint(d);;
}

//-------------------------
bool InitialStateParticle::addDataPoint(const DataPoint& d)
{
    return addDataPoint(DataPoint(d));
}


}

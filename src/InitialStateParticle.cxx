#include "InitialStateParticle.h"

#include "DataSet.h"
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
    addDataAccessor(this);
    removeDataAccessor(&HelicityAngles_);
}

//-------------------------
bool InitialStateParticle::consistent() const
{
    return DecayingParticle::consistent();
}

//-------------------------
bool InitialStateParticle::prepare()
{
    if (!consistent()) {
        LOG(ERROR) << "Cannot prepare InitialStateParticle, it is not consistent.";
        return false;
    }

    ParticleCombination::makeParticleCombinationSetWithParents();

    setSymmetrizationIndexParents();
    optimizeSpinAmplitudeSharing();

    // make sure that final state particles get the correct indices
    for (unsigned i = 0; i < particleCombinations()[0]->indices().size(); ++i) {
        std::shared_ptr<ParticleCombination> index = std::make_shared<ParticleCombination>(i);
        for (std::shared_ptr<ParticleCombination> pc : ParticleCombination::particleCombinationSet()) {
            if (ParticleCombination::equivByOrderlessContent(index, pc)) {
                FourMomenta_.addSymmetrizationIndex(pc);
                continue;
            }
        }
    }

    // add particle combinations to FourMomenta_ and HelicityAngles_
    for (auto& pc : ParticleCombination::particleCombinationSet()) {
        FourMomenta_.addSymmetrizationIndex(pc);

        if (pc->indices().size() > 1)
            HelicityAngles_.addSymmetrizationIndex(pc);
    }

    FourMomenta_.findInitialStateParticle();

    setDataAcessorIndices();

    if (!consistent()) {
        LOG(ERROR) << "Something went wrong while preparing InitialStateParticle, it is not consistent.";
        return false;
    }

    return true;
}

//-------------------------
void InitialStateParticle::setSymmetrizationIndexParents()
{
    //std::cout << "InitialStateParticle::setSymmetrizationIndexParents()\n";

    unsigned size = particleCombinations()[0]->indices().size();
    assert(size > 1);

    clearSymmetrizationIndices();

    // get initial state PCs from set and add them
    for (auto& pc : ParticleCombination::particleCombinationSet()) {
        if (pc->indices().size() < size)
            continue;

        if (pc->daughters()[0]->parent() == nullptr)
            continue;

        //std::cout << "  add " << std::string(*pc) << "\n";
        addSymmetrizationIndex(pc);
    }

    // next level
    for (auto& ch : channels())
        ch->setSymmetrizationIndexParents();

}

//-------------------------
bool InitialStateParticle::addDataPoint(DataPoint&& d)
{
    d.allocateStorage(FourMomenta_, HelicityAngles_, DataAccessors_);

    FourMomenta_.calculate(d);
    HelicityAngles_.calculate(d);

    /*for (auto& pc : FourMomenta_.particleCombinations()) {
      std::cout << FourMomenta_.symmetrizationIndex(pc) << "  \t" << std::string(*pc) << "  \t";
      FourMomenta_.p(d, FourMomenta_.symmetrizationIndex(pc)).Print();
    }

    for (auto& pc : HelicityAngles_.particleCombinations()) {
      std::cout << HelicityAngles_.symmetrizationIndex(pc) << "  \t" << std::string(*pc) << "  \t";
      for (double a : HelicityAngles_.helicityAngles(d, HelicityAngles_.symmetrizationIndex(pc)))
          std::cout << a << "  \t";
      std::cout << "\n";
    }*/

    if (!DataSet_.consistent(d))
        return false;

    return DataSet_.addDataPoint(d);;
}

//-------------------------
bool InitialStateParticle::addDataPoint(const DataPoint& d)
{
    return addDataPoint(DataPoint(d));
}

//-------------------------
void InitialStateParticle::printDataAccessors()
{
    std::cout << "DataAccessors of " << name() << "\n"
              << "index \tnSymIndices \tname  \t\tparticleCombinations\n";
    for (DataAccessor* d : DataAccessors_) {
        std::cout << d->index() << "  \t" << d->maxSymmetrizationIndex() + 1
                  << "  \t(" << typeid(*d).name() << ")  \t";
        for (auto& pc : d->particleCombinations())
            std::cout << std::string(*pc) << ";  ";
        std::cout << "\n";
    }
    std::cout << std::endl;
}

//-------------------------
void InitialStateParticle::setDataAcessorIndices()
{
    unsigned i(0);
    for (DataAccessor* d : DataAccessors_)
        d->setIndex(i++);
}


}

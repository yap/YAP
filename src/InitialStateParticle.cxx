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
    //std::cout << "InitialStateParticle::setSymmetrizationIndexParents()\n";

    ParticleCombination::makeParticleCombinationSetWithParents();

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

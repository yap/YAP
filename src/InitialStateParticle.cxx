#include "InitialStateParticle.h"

#include "logging.h"

#include <TLorentzRotation.h>

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
    for (std::shared_ptr<yap::ParticleCombination>& pc : particleCombinations()) {
        pc->setParents();
    }
}

//-------------------------
bool InitialStateParticle::addDataPoint(DataPoint&& d)
{
    FourMomenta_.calculate(d);
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

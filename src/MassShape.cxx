#include "MassShape.h"

#include "logging.h"
#include "ParticleCombination.h"

namespace yap {

//-------------------------
MassShape::MassShape(std::initializer_list<double> pars) :
    DataAccessor(&ParticleCombination::equivByOrderlessContent),
    Parameters_(pars)
{
}

//-------------------------
bool MassShape::setParameters(const ParticleTableEntry& entry)
{
    if (entry.MassShapeParameters_.size() < Parameters_.size()) {
        LOG(ERROR) << "MassShape::setParameters - entry.MassShapeParameters_ is too small for MassShape object.";
        return false;
    }
    for (unsigned i = 0; i < entry.MassShapeParameters_.size(); ++i) {
        Parameters_.at(i)->setValue(entry.MassShapeParameters_[i]);
    }
    return true;
}

}

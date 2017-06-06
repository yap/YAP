#include "MassShape.h"

#include "CachedValue.h"
#include "CalculationStatus.h"
#include "DataPartition.h"
#include "DecayingParticle.h"
#include "Exceptions.h"
#include "logging.h"
#include "Parameter.h"
#include "ParticleCombination.h"
#include "VariableStatus.h"

namespace yap {

//-------------------------
MassShape::MassShape() :
    RecalculableAmplitudeComponent(equal_by_orderless_content),
    Model_(nullptr)
{}

//-------------------------
void MassShape::calculate(DataPartition& D) const
{
    // loop over (ParticleCombination --> symmetrization index) map
    for (const auto& pc_si : symmetrizationIndices())
        calculate(D, pc_si.first, pc_si.second);
}

//-------------------------
bool MassShape::consistent() const
{
    bool C = DataAccessor::consistent();

    // check if owner is set
    if (!Model_) {
        FLOG(ERROR) << "Model isn't set";
        C &= false;
    }

    return C;
}

//-------------------------
void MassShape::setModel(const Model& m)
{
    Model_ = &m;
}

}

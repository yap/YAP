#include "PhaseSpaceFactor.h"

#include "TwoBodyPhaseSpaceFactor.h"

namespace yap {

//-------------------------
std::shared_ptr<PhaseSpaceFactorFactory> DefaultPHSPFactory = TwoBodyPhaseSpaceFactorFactory::instance();

}

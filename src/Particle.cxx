#include "Particle.h"
#include "logging.h"

namespace yap {

//-------------------------
bool Particle::consistent() const
{

    if (Mass_ <= 0.) {
        LOG(ERROR) << "Particle::consistent() - mass not positive.";
        return false;
    }

    return true;
}

}

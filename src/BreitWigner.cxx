#include "BreitWigner.h"

#include "Constants.h"
#include "logging.h"

namespace yap {

//-------------------------
BreitWigner::BreitWigner(double mass, double width) :
    MassShape(2)
{
    setMass(mass);
    setWidth(width);
}

//-------------------------
Amp BreitWigner::amplitude(DataPoint& d)
{
    return Complex_0;
}

//-------------------------
Amp BreitWigner::amplitude(double s)
{
    return 1. / Amp(squaredmass(s) - s, mass() * width(s));
}

//-------------------------
bool BreitWigner::consistent() const
{
    bool consistent = true;
    if (mass() <= 0) {
        LOG(ERROR) << "BreitWigner::consistent() - mass <= 0";
        consistent = false;
    }
    if (width() <= 0) {
        LOG(ERROR) << "BreitWigner::consistent() - width <= 0";
        consistent = false;
    }

    return consistent;
}

}





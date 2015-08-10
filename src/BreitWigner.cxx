#include "BreitWigner.h"

#include "Constants.h"

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
    return (mass() > 0 && width() > 0);
}

}




